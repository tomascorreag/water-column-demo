// Water Column Console — demo v3 (data-first rework).
//
// Encoding rules:
//   colour      = the selected physical variable (temp / salinity / sigma0),
//                 or the support proxy when "Support" is chosen. One legend.
//   solid/holes = support of the interpolated field: solid where kernel
//                 support >= 0.5, screen-door stipple 0.08..0.5, absent below.
//   geometry    = where measurements are (tubes), where a density surface
//                 sits (isopycnal sheets), where the field is cut (slice/section).
//   accents     = pinned profiles only (orange / pink / green); never data.
// Casts are drawn from binned measurements; slices, sections and isopycnal
// sheets are drawn from the interpolated grid and say so in the readout.

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// ---------- showcase build ----------
// A reduced build for sharing: casts, the main thermocline, time, and nothing
// else. Same code and same data as the console — the difference is what is
// exposed. build_demo_page.py --showcase stamps the class on <body>;
// ?showcase=1 does the same when serving demo/ directly.
const SHOWCASE = document.body.classList.contains('showcase')
  || new URLSearchParams(location.search).has('showcase');
if (SHOWCASE) document.body.classList.add('showcase', 'explain');
if (SHOWCASE) {
  // every panel starts hidden; the introduction is the only thing that reveals
  // them, one step at a time
  for (const id of ['title', 'exag-badge', 'guide', 'tour', 'row-casts', 'row-thermo',
    'row-bathy', 'insets', 'timeline', 'legend', 'compass'])
    document.getElementById(id).classList.add('rv');
  // where each panel sits is style.css's business (one grid, grid-area per
  // mode and breakpoint): the guide has the top of the left column, the
  // toggles, the colour legend and the plots stack at its foot, and the
  // timeline has the rest of the bottom edge
  document.querySelector('#title h1').textContent = 'Water column walkthrough';
  document.querySelector('#title .sub').textContent =
    '83 Argo profiles · north-west Atlantic · June 2023 · rotate with the mouse, zoom with the wheel';
}

// ---------- small screen ----------
// One definition of "small", shared by the layout (style.css keys on
// body.mobile and nothing there keys on a width), the input model (a touch
// screen has no hover, so nothing may be told only by hovering) and the
// notice. The short-viewport half is qualified by a coarse pointer so that a
// desktop window dragged flat does not become a phone. ?mobile=1 forces it for
// testing on a desktop, ?mobile=0 suppresses it so a scripted capture at a
// small window still gets the full console.
const MOBILE_MQ = matchMedia('(max-width: 820px), ((max-height: 460px) and (pointer: coarse))');
const MOBILE_FORCE = new URLSearchParams(location.search).get('mobile');
let MOBILE = MOBILE_FORCE === '1' ? true : MOBILE_FORCE === '0' ? false : MOBILE_MQ.matches;
const MOBILE_AT_LOAD = MOBILE;
document.body.classList.toggle('mobile', MOBILE);

// ---------- data decoding ----------

function decodeU16(b64, lo, hi) {
  const bin = atob(b64);
  const n = bin.length >> 1;
  const out = new Float32Array(n);
  const span = (hi - lo) / 65000;
  for (let i = 0; i < n; i++) {
    const u = bin.charCodeAt(2 * i) | (bin.charCodeAt(2 * i + 1) << 8);
    out[i] = lo + u * span;
  }
  return out;
}

const G = OCEAN.grid;
const R = OCEAN.rng;
const casts = OCEAN.casts.map(c => ({
  ...c,
  pres: decodeU16(c.pres, ...R.pres),
  temp: decodeU16(c.temp, ...R.temp),
  psal: decodeU16(c.psal, ...R.psal),
  sigma0: decodeU16(c.sigma0, ...R.sigma0),
  // stratification, shipped as the buoyancy frequency N (TEOS-10, smoothed
  // over 50 dbar in build_demo_data.py) and squared back here: N^2 is what the
  // pycnocline detector reads.
  n2: decodeU16(c.nfreq, ...R.nfreq).map(v => v * v),
  // per-grid-level kernel sums: sv = weight this cast puts on level k,
  // vb* = the value it carries there. Enough to re-assemble the whole field
  // under any time weighting without shipping a 4-D grid.
  sv: decodeU16(c.sv, 0, OCEAN.time.svMax),
  vb: {
    temp: decodeU16(c.vbTemp, ...R.temp),
    psal: decodeU16(c.vbPsal, ...R.psal),
    sigma0: decodeU16(c.vbSigma0, ...R.sigma0),
  },
}));
const NX = G.nx, NY = G.ny, NZ = G.nz;
const gzAt = k => k * G.presMax / (NZ - 1);

const NCELL = NX * NY * NZ;
const ACC = { temp: new Float32Array(NCELL), psal: new Float32Array(NCELL), sigma0: new Float32Array(NCELL), den: new Float32Array(NCELL) };

// per-cast horizontal footprint: the Gaussian is cut at 3 L_h (weight 0.011),
// which turns the field rebuild from every node into a small box per cast
const HCUT = 3 * G.lh;
const dxCell = (G.lon1 - G.lon0) / (NX - 1) * G.kmLon;
const dyCell = (G.lat1 - G.lat0) / (NY - 1) * G.kmLat;
const foot = casts.map(c => {
  const cx = (c.lon - G.lon0) / (G.lon1 - G.lon0) * (NX - 1);
  const cy = (c.lat - G.lat0) / (G.lat1 - G.lat0) * (NY - 1);
  const x0 = Math.max(0, Math.ceil(cx - HCUT / dxCell)), x1 = Math.min(NX - 1, Math.floor(cx + HCUT / dxCell));
  const y0 = Math.max(0, Math.ceil(cy - HCUT / dyCell)), y1 = Math.min(NY - 1, Math.floor(cy + HCUT / dyCell));
  const w = new Float32Array(Math.max(0, (x1 - x0 + 1) * (y1 - y0 + 1)));
  let i = 0;
  for (let y = y0; y <= y1; y++) {
    const dky = (y - cy) * dyCell;
    for (let x = x0; x <= x1; x++) {
      const dkx = (x - cx) * dxCell;
      w[i++] = Math.exp(-(dkx * dkx + dky * dky) / (2 * G.lh * G.lh));
    }
  }
  return { x0, x1, y0, y1, w };
});

// Assemble the field from the per-cast sums: node (k, y, x) takes
// sum_c a_c[k] sv_c[k] wh_c[y, x] in the denominator and the same times
// vbar_c[k] in the numerator. weightAt(c, a) fills a[k] with cast c's weight
// per level and returns false if the cast contributes nowhere; awMin is a
// speed-only cutoff on a[k] sv[k].
function accumulateField(weightAt, awMin) {
  ACC.temp.fill(0); ACC.psal.fill(0); ACC.sigma0.fill(0); ACC.den.fill(0);
  const a = new Float32Array(NZ);
  for (let ci = 0; ci < casts.length; ci++) {
    const c = casts[ci];
    if (!weightAt(c, a)) continue;
    const fp = foot[ci];
    for (let k = 0; k < NZ; k++) {
      const aw = c.sv[k] * a[k];
      if (aw <= awMin) continue;
      const at = aw * c.vb.temp[k], ap = aw * c.vb.psal[k], as = aw * c.vb.sigma0[k];
      const base = k * NY * NX;
      for (let y = fp.y0, wi = 0; y <= fp.y1; y++) {
        const row = base + y * NX;
        for (let x = fp.x0; x <= fp.x1; x++, wi++) {
          const w = fp.w[wi];
          const o = row + x;
          ACC.den[o] += aw * w;
          ACC.temp[o] += at * w; ACC.psal[o] += ap * w; ACC.sigma0[o] += as * w;
        }
      }
    }
  }
  const ref = G.ref, mean = G.mean;
  for (let i = 0; i < NCELL; i++) {
    const d = ACC.den[i];
    F.conf[i] = d > 0 ? Math.min(d / ref, 1) : 0;
    if (d > 0) {
      F.temp[i] = ACC.temp[i] / d; F.psal[i] = ACC.psal[i] / d; F.sigma0[i] = ACC.sigma0[i] / d;
    } else {
      F.temp[i] = mean.temp; F.psal[i] = mean.psal; F.sigma0[i] = mean.sigma0;
    }
  }
}

// The pooled-month field is not shipped: it is the timeline's own sum with
// every cast at weight 1, so it is built here once (~15 ms) rather than
// carried as 2 MB of base64. build_demo_data.py still forms the same grid,
// without the 3 L_h cut, to fix `ref` and `mean`; the two agree to ~1 mK
// except at nodes only a far tail reaches (support ~0.02-0.05), where this
// one, the cut kernel, is the documented method and the one the timeline uses.
const F = { temp: new Float32Array(NCELL), psal: new Float32Array(NCELL), sigma0: new Float32Array(NCELL), conf: new Float32Array(NCELL) };
accumulateField((c, a) => { a.fill(1); return true; }, 0);

// ---------- coordinates (km, y up; world.scale.y = vertical exaggeration) ----------

const lonMid = (G.lon0 + G.lon1) / 2, latMid = (G.lat0 + G.lat1) / 2;
const xOf = lon => (lon - lonMid) * G.kmLon;
const zOf = lat => -(lat - latMid) * G.kmLat;
const yOf = pres => -pres / 1000;                 // dbar ~ m -> km
const lonOfX = x => lonMid + x / G.kmLon;
const latOfZ = z => latMid - z / G.kmLat;
const DEEP = G.presMax / 1000;

let RC = 0;                                        // circular domain, circumscribes casts
for (const c of casts) RC = Math.max(RC, Math.hypot(xOf(c.lon), zOf(c.lat)));
RC += 60;
const clamp01 = v => Math.max(0, Math.min(1, v));

// ---------- time state (declared early: cast colouring reads it at build) ----------
// The heavy lifting is in "time: the field as it was known at an instant" below.

const T = OCEAN.time;
const TIME = {
  on: false,
  t: T.tMax,                                       // playhead, days since T.t0
  scale: 1,                                        // multiplier on L_t (expert)
  playing: false,
  speed: 3.0,                                      // days per second
};
const LT = Float32Array.from(T.lt);                // L_t per grid level, days
const ltAt = pres => {
  const f = clamp01(pres / G.presMax) * (NZ - 1);
  const k = Math.min(Math.floor(f), NZ - 2);
  return (LT[k] + (LT[k + 1] - LT[k]) * (f - k)) * TIME.scale;
};
// How much a cast says about the water at the playhead, at one pressure.
// Two-sided: the record is already collected, so the estimate at t is a
// smoother, not a filter, and a cast taken just after t constrains t as well
// as one taken just before. Without that, a cast switched on at full weight
// the instant it was taken and the field jumped between frames.
//
// The uncertainty weight deliberately does NOT appear here. This kernel is
// what makes the field an estimate *of an instant* — drop it and the playhead
// stops meaning anything, which is what an earlier version did at UNC = 0.
// What UNC suppresses is the *display* of the resulting support: at 0 the
// field still moves with the playhead, it just no longer admits that a cast
// three weeks away is a weak claim on it.
const timeWeight = (dt, L) => Math.exp(-dt * dt / (2 * L * L));
function ageWeight(c, pres) {
  if (!TIME.on) return 1;
  return timeWeight(TIME.t - c.t, ltAt(pres));
}
const castExists = () => true;
const dayLabel = d => {                            // days since t0 -> "Jun 14 06:00Z"
  const dt = new Date(Date.parse(T.t0.replace('Z', ':00Z')) + d * 86400000);
  const mon = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'][dt.getUTCMonth()];
  return `${mon} ${dt.getUTCDate()} ${String(dt.getUTCHours()).padStart(2, '0')}:00Z`;
};

// ---------- colour: variable -> colormap ----------

// matplotlib LUTs, 32 stops. Inferno is sampled from 0.12 so the cold end is
// dark purple, not black, and stays visible against the mid-grey ground.
// Bright ends are clipped (~0.9) so the warmest/saltiest water does not vanish
// into the paper-grey ground.
const LUT = {
  inferno: [[0.082, 0.043, 0.215], [0.136, 0.047, 0.3], [0.19, 0.039, 0.361], [0.245, 0.037, 0.4], [0.297, 0.047, 0.42], [0.354, 0.067, 0.431], [0.404, 0.086, 0.433], [0.454, 0.104, 0.43], [0.503, 0.122, 0.423], [0.56, 0.141, 0.41], [0.609, 0.159, 0.394], [0.658, 0.179, 0.373], [0.706, 0.201, 0.348], [0.758, 0.229, 0.315], [0.802, 0.259, 0.283], [0.842, 0.293, 0.249], [0.878, 0.332, 0.212], [0.913, 0.382, 0.17], [0.939, 0.43, 0.13], [0.959, 0.482, 0.089], [0.974, 0.537, 0.048], [0.985, 0.601, 0.024], [0.988, 0.66, 0.052], [0.986, 0.721, 0.112], [0.977, 0.782, 0.186], [0.963, 0.851, 0.286]],
  viridis: [[0.267, 0.005, 0.329], [0.277, 0.05, 0.376], [0.282, 0.095, 0.417], [0.283, 0.136, 0.453], [0.278, 0.18, 0.487], [0.269, 0.219, 0.51], [0.257, 0.256, 0.527], [0.243, 0.292, 0.539], [0.226, 0.331, 0.547], [0.211, 0.364, 0.552], [0.196, 0.395, 0.555], [0.182, 0.426, 0.557], [0.168, 0.46, 0.558], [0.156, 0.49, 0.558], [0.145, 0.519, 0.557], [0.134, 0.549, 0.554], [0.123, 0.582, 0.547], [0.119, 0.611, 0.539], [0.125, 0.64, 0.527], [0.143, 0.669, 0.511], [0.181, 0.701, 0.488], [0.226, 0.729, 0.463], [0.281, 0.755, 0.433], [0.344, 0.78, 0.397], [0.422, 0.806, 0.352], [0.497, 0.826, 0.306], [0.576, 0.845, 0.256], [0.658, 0.86, 0.203], [0.752, 0.875, 0.143], [0.835, 0.886, 0.103]],
  // YlGnBu reversed again = light (buoyant) -> dark blue (dense)
  ylgnbu: [[0.924, 0.97, 0.695], [0.887, 0.956, 0.697], [0.849, 0.941, 0.7], [0.812, 0.926, 0.703], [0.76, 0.906, 0.708], [0.69, 0.878, 0.715], [0.619, 0.851, 0.722], [0.548, 0.823, 0.728], [0.472, 0.794, 0.737], [0.411, 0.772, 0.746], [0.35, 0.749, 0.755], [0.289, 0.726, 0.764], [0.235, 0.693, 0.766], [0.2, 0.657, 0.762], [0.164, 0.62, 0.759], [0.129, 0.584, 0.755], [0.117, 0.533, 0.736], [0.122, 0.483, 0.713], [0.127, 0.433, 0.689], [0.132, 0.383, 0.665], [0.136, 0.334, 0.642], [0.139, 0.292, 0.623], [0.142, 0.251, 0.603], [0.145, 0.21, 0.583], [0.117, 0.182, 0.522], [0.088, 0.159, 0.463], [0.06, 0.136, 0.404], [0.031, 0.114, 0.345]],
  // support: pale -> dark slate; single hue, "more = darker" on a light ground
  support: [[0.86, 0.90, 0.93], [0.66, 0.78, 0.85], [0.45, 0.62, 0.71], [0.33, 0.46, 0.55], [0.24, 0.31, 0.38], [0.14, 0.17, 0.21]],
  // depth of a mapped surface: pale violet -> deep indigo, "deeper = darker".
  // Its own hue on purpose — inferno, viridis, YlGnBu and slate already mean
  // temperature, salinity, density and support, and a mapped depth is none of
  // those. It keeps the indigo the pycnocline volume has always been drawn in.
  // rose (shallow) -> indigo (deep). Both ends have to be dark enough to show
  // through a sheet drawn at a third opacity on a light ground, so this ramp
  // carries depth mostly in hue and only partly in lightness; a pale shallow
  // end simply disappeared. It stops short of orange so it cannot be read as
  // the temperature ramp.
  depth: [[0.851, 0.353, 0.545], [0.729, 0.286, 0.549], [0.608, 0.243, 0.545], [0.482, 0.216, 0.522], [0.365, 0.196, 0.478], [0.255, 0.176, 0.404], [0.149, 0.129, 0.318]],
  // thickness of a layer: thin -> thick, pale sage to deep green. It can start
  // pale because the sheet it colours is drawn opaque, so nothing shows through
  // it; and it is the one ramp here with no other meaning attached.
  thick: [[0.906, 0.929, 0.882], [0.788, 0.851, 0.769], [0.643, 0.769, 0.639], [0.478, 0.667, 0.522], [0.325, 0.549, 0.427], [0.192, 0.416, 0.333], [0.086, 0.271, 0.235]],
};
function lutColor(lut, u) {
  u = clamp01(u);
  const f = u * (lut.length - 1);
  const i = Math.min(Math.floor(f), lut.length - 2), w = f - i;
  return [0, 1, 2].map(k => lut[i][k] + (lut[i + 1][k] - lut[i][k]) * w);
}
function niceTicks(lo, hi, n = 5) {
  const raw = (hi - lo) / n, p = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * p).find(s => (hi - lo) / s <= n + 1);
  const t = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) t.push(+v.toFixed(6));
  return { ticks: t, step };
}
const VARS = {
  temp: { name: 'Temperature (in situ)', unit: '°C', lut: LUT.inferno, lo: R.temp[0], hi: R.temp[1], contour: 2, dec: 1 },
  psal: { name: 'Practical salinity', unit: 'PSS-78', lut: LUT.viridis, lo: R.psal[0], hi: R.psal[1], contour: 0.25, dec: 2 },
  sigma0: { name: 'Potential density anomaly σ₀ (TEOS-10)', unit: 'kg m⁻³', lut: LUT.ylgnbu, lo: R.sigma0[0], hi: R.sigma0[1], contour: 0.25, dec: 2 },
  conf: { name: 'Field support (kernel-weight proxy)', unit: '0–1', lut: LUT.support, lo: 0, hi: 1, contour: 0.25, dec: 2 },
};
let curVar = 'temp';
const V = () => VARS[curVar];
const varColor = (v, key = curVar) => {
  const d = VARS[key];
  return lutColor(d.lut, (v - d.lo) / (d.hi - d.lo));
};
const css = rgb => `rgb(${rgb.map(v => Math.round(v * 255)).join(',')})`;

// support thresholds for the interpolated field
let CUT = 0.0;                                     // not drawn below (default: draw everything)
const SOLID = 0.5;                                 // solid at/above
let stippleOn = false;
let contoursOn = true;

// ---------- uncertainty weight: how much of the doubt is drawn ----------
//
// UNC = 1 is the honest render: support drives the VSUP tree, the opacity and
// the temporal decay. UNC = 0 is the render an ordinary pipeline would give
// you — every point that any cast reaches is drawn at full colour and full
// opacity, and every cast counts equally whatever its age. The slider is
// there so the two can be compared on the same data rather than described;
// it suppresses the *display* of uncertainty, it does not reduce it, and
// anything below 1 raises a badge saying so.
//
//   support seen by the encoding : s' = 1 - UNC (1 - s)
//
// It acts on the encoding only. The spatial and temporal kernels are left
// alone, so at UNC = 0 the field still varies across the box and still moves
// with the playhead — it is drawn everywhere the kernel reaches, at full
// colour and full opacity, as an ordinary pipeline would draw it.
let UNC = 1;
const supEff = s => 1 - UNC * (1 - s);

// VSUP (Correll, Moritz & Heer 2018): a quantization tree over (value,
// uncertainty). Layer 0 (support >= SOLID) has 2^(L-1) value bins in the
// full ramp; each layer up halves the bins and moves toward a neutral grey,
// so low-support field is drawn with less value precision, not just flagged.
const VSUP_LAYERS = 5;
const VSUP_NEUTRAL = [0.48, 0.51, 0.55];             // #7a828c, darker than the ground
function vsupLayer(sup) {
  sup = supEff(sup);
  if (sup >= SOLID) return 0;
  return 1 + Math.min(VSUP_LAYERS - 2, Math.floor(clamp01((SOLID - sup) / (SOLID - CUT)) * (VSUP_LAYERS - 1)));
}
// the same tree on an arbitrary ramp and an already-normalised value, for the
// layer surfaces (whose value is a thickness, not one of the field variables)
// The tree has 31 leaves per ramp (16+8+4+2+1 bins), so every colour an
// interpolated layer can take is one of 31 arrays per ramp, built on first
// use and shared after that. UNC is folded in by vsupLayer before the lookup,
// so (ramp, layer, bin) is the whole key. The array returned is shared: read
// it, never write to it. The slice alone asks 40k times per rebuild, and
// allocating two arrays per ask was most of what it cost.
const vsupMemo = new Map();                         // lut -> [layer * 16 + bin] -> [r, g, b]
function vsupLeaf(lut, layer, b) {
  let tab = vsupMemo.get(lut);
  if (!tab) { tab = new Array(VSUP_LAYERS * 16).fill(null); vsupMemo.set(lut, tab); }
  const slot = layer * 16 + b;
  let rgb = tab[slot];
  if (!rgb) {
    const bins = 1 << (VSUP_LAYERS - 1 - layer);
    const m = layer / (VSUP_LAYERS - 1);              // 0 = full colour, 1 = neutral
    rgb = tab[slot] = lutColor(lut, (b + 0.5) / bins).map((c, i) => c + (VSUP_NEUTRAL[i] - c) * m);
  }
  return rgb;
}
function vsupOn(lut, u, sup) {
  const layer = vsupLayer(sup);
  const bins = 1 << (VSUP_LAYERS - 1 - layer);
  return vsupLeaf(lut, layer, Math.min(bins - 1, Math.floor(clamp01(u) * bins)));
}
function vsupColor(v, sup, key = curVar) {
  const d = VARS[key];
  const layer = vsupLayer(sup);
  const bins = 1 << (VSUP_LAYERS - 1 - layer);
  const u = clamp01((v - d.lo) / (d.hi - d.lo));
  return vsupLeaf(d.lut, layer, Math.min(bins - 1, Math.floor(u * bins)));
}
// opacity: interpolated geometry fades out entirely as support -> 0. Support
// exactly 0 means no cast is within the kernel cutoff at all, so nothing is
// drawn there at any uncertainty weight — UNC removes the graded fade, not
// the distinction between "thin" and "nothing".
//
// EDGE is the last sliver of that: the horizontal kernel is cut at 3 L_h, so
// support does not decay to zero, it stops. At UNC = 1 the fade hides that,
// but at 0 the cut would be drawn at full opacity as a cell-aligned staircase
// — the shape of the grid, not of anything measured. Ramping the alpha over
// the last 0.02 of support puts the edge on the support contour instead.
const EDGE = 0.02;
const supportAlpha = sup => (sup <= 0 ? 0
  : Math.pow(clamp01(supEff(sup) / SOLID), 0.75) * clamp01(sup / EDGE));
// support range covered by each layer, for the legend
function vsupRange(layer) {
  if (layer === 0) return [SOLID, 1];
  const step = (SOLID - CUT) / (VSUP_LAYERS - 1);
  return [SOLID - layer * step, SOLID - (layer - 1) * step];
}

const castSigma = c => c.terr ?? ({ D: 0.002, A: 0.005, R: 0.02 }[c.mode] ?? 0.02);
const modeName = m => ({ D: 'delayed-mode (calibrated)', A: 'adjusted', R: 'real-time (uncalibrated)' }[m] ?? m);

// ---------- scene ----------

THREE.ColorManagement.enabled = false;             // hex/float colours are display values as written
const container = document.getElementById('scene');
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, MOBILE ? 1.5 : 2));
renderer.setSize(innerWidth, innerHeight);
renderer.outputColorSpace = THREE.LinearSRGBColorSpace;   // colours are display values
renderer.setClearColor(0xe4e7eb, 1);
container.appendChild(renderer.domElement);

const scene = new THREE.Scene();
const persp = new THREE.PerspectiveCamera(42, innerWidth / innerHeight, 1, 60000);
const ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, -60000, 60000);
let camera = persp;                                // swapped to ortho on gizmo snaps
// The showcase HUD owns the bottom and the left, so a scene centred in the
// viewport sits behind the panels. setViewOffset shifts the rendered window
// rather than the camera, so orbiting still turns about the domain centre and
// the projected DOM labels follow without a second correction.
// The constants were tuned at 1920x1080, where the left column takes 0.33 of
// the width and the bottom band 0.24 of the height. On any other viewport the
// shift scales with the fraction the HUD actually takes, read off the laid-out
// panels: a .rv panel is laid out at opacity 0, so this is valid before
// anything is revealed.
const SC_SHIFT_X = 0.045, SC_SHIFT_Y = 0.10;        // fractions of the viewport: right, up
const SC_REF = { l: 0.33, b: 0.24 };
const scBase = { x: SC_SHIFT_X, y: SC_SHIFT_Y };    // the working view's shift on this viewport
const scShift = { x: SC_SHIFT_X, y: SC_SHIFT_Y };   // live: the opening frames the chart elsewhere
function computeViewShift() {
  if (!SHOWCASE) return;
  // on a phone the HUD is one column across the bottom, so there is no left
  // column to shift out of; the plots are inside a sheet that is translated
  // off-screen while closed and would report a meaningless rect anyway
  if (MOBILE) { scBase.x = 0; scBase.y = 0.12; return; }
  const rect = id => {
    const el = document.getElementById(id);
    return el && getComputedStyle(el).display !== 'none' ? el.getBoundingClientRect() : null;
  };
  let l = 0, top = innerHeight, right = innerWidth;
  for (const id of ['tour', 'insets']) { const r = rect(id); if (r && r.width) l = Math.max(l, r.right); }
  for (const id of ['insets', 'timeline']) { const r = rect(id); if (r && r.height) top = Math.min(top, r.top); }
  // the walkthrough card is in the right column, so it takes room off the
  // other side: the shift is what the left column occupies *net of* it, and it
  // can reach zero -- with both edges covered the scene belongs in the middle
  { const r = rect('guide'); if (r && r.width) right = Math.min(right, r.left); }
  const lf = (l - (innerWidth - right)) / innerWidth, bf = (innerHeight - top) / innerHeight;
  scBase.x = Math.min(0.12, Math.max(0, SC_SHIFT_X * lf / SC_REF.l));
  scBase.y = Math.min(0.16, Math.max(0.05, SC_SHIFT_Y * bf / SC_REF.b));
}
function applyViewShift() {
  if (!SHOWCASE) return;
  const w = innerWidth, h = innerHeight;
  for (const c of [persp, ortho]) c.setViewOffset(w, h, -w * scShift.x, h * scShift.y, w, h);
}
applyViewShift();
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.enablePan = !SHOWCASE;         // showcase navigation is rotate, zoom and time
function setView(name) {
  const views = {
    overview: [[RC * 0.7, RC * 1.0, RC * 1.9], [0, -120, 0]],
    map: [[0, RC * 2.9, 1], [0, -50, 0]],
    section: [[RC * 1.05, RC * 0.42, RC * 0.6], [0, -110, 0]],
    low: [[RC * 1.6, RC * 0.35, RC * 1.6], [0, -100, 0]],
    density: [[RC * 1.35, RC * 0.7, RC * 1.45], [0, -110, 0]],
  }[name];
  if (camera === ortho) goPersp();
  camera.position.set(...views[0]);
  controls.target.set(...views[1]);
  if (name === 'section' && typeof sect !== 'undefined') {   // centre on the plane
    const dx = sect.mesh.position.x, dz = sect.mesh.position.z;
    camera.position.x += dx; camera.position.z += dz; controls.target.x += dx; controls.target.z += dz;
  }
  controls.update();
}
setView('overview');

const world = new THREE.Group();                  // world.scale.y = exaggeration
scene.add(world);

// ---------- labels projected from 3D to HTML ----------

const labelRoot = document.getElementById('labels');
const labels = [];                                 // { el, pos: Vector3 (pre-exag), fixedY? }
function addLabel(text, x, y, z, cls = '') {
  const el = document.createElement('div');
  el.className = 'lb ' + cls;
  el.textContent = text;
  labelRoot.appendChild(el);
  const L = { el, pos: new THREE.Vector3(x, y, z), visible: true };
  labels.push(L);
  return L;
}
const _v = new THREE.Vector3();
function updateLabels() {
  const w = innerWidth, h = innerHeight;
  for (const L of labels) {
    if (!L.visible) { L.el.style.display = 'none'; continue; }
    _v.set(L.pos.x, L.pos.y * exag, L.pos.z).project(camera);
    if (_v.z > 1 || Math.abs(_v.x) > 1.05 || Math.abs(_v.y) > 1.05) { L.el.style.display = 'none'; continue; }
    L.el.style.display = 'block';
    L.el.style.left = ((_v.x + 1) / 2 * w) + 'px';
    L.el.style.top = ((1 - _v.y) / 2 * h) + 'px';
  }
}

// ---------- reference frame: graticule, depth axis, scale bar, north ----------
// grouped so ?grid=off can drop the whole frame for concept stills that want
// the surface alone (the DOM labels go with it via body.nogrid)
const frameGroup = new THREE.Group();
world.add(frameGroup);
// the surface graticule and its lon/lat labels, kept apart from the rest of the
// frame so the showcase can drop just those (see setGraticule)
let gratMesh = null;
const gratLabels = [];
// seen from straight above, the depth axis collapses into a stack of labels on
// top of each other and says nothing: there is no depth in a map view yet
const depthLabels = [];
function setDepthAxis(on) { for (const L of depthLabels) L.visible = on; }
{
  const gridMat = new THREE.LineBasicMaterial({ color: 0xb3bac3 });
  const rimMat = new THREE.LineBasicMaterial({ color: 0x6b737d });
  const pts = [];
  const inCircle = (x, z) => Math.hypot(x, z) <= RC;
  // graticule every 2 deg on the surface, clipped to the circle
  const lonA = Math.ceil(lonOfX(-RC) / 2) * 2, lonB = Math.floor(lonOfX(RC) / 2) * 2;
  const latA = Math.ceil(latOfZ(RC) / 2) * 2, latB = Math.floor(latOfZ(-RC) / 2) * 2;
  for (let lon = lonA; lon <= lonB; lon += 2) {
    const x = xOf(lon);
    if (Math.abs(x) >= RC) continue;
    const zz = Math.sqrt(RC * RC - x * x);
    pts.push(new THREE.Vector3(x, 0, -zz), new THREE.Vector3(x, 0, zz));
    gratLabels.push(addLabel(`${Math.abs(lon)}°W`, x, 0, zz + 25));
  }
  for (let lat = latA; lat <= latB; lat += 2) {
    const z = zOf(lat);
    if (Math.abs(z) >= RC) continue;
    const xx = Math.sqrt(RC * RC - z * z);
    pts.push(new THREE.Vector3(-xx, 0, z), new THREE.Vector3(xx, 0, z));
    gratLabels.push(addLabel(`${lat}°N`, xx + 30, 0, z));
  }
  gratMesh = new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(pts), gridMat);
  frameGroup.add(gratMesh);
  // rim rings at surface and floor
  const ring = (r, y, mat) => {
    const p = [];
    for (let i = 0; i <= 128; i++) {
      const a = i / 128 * Math.PI * 2;
      p.push(new THREE.Vector3(r * Math.cos(a), y, r * Math.sin(a)));
    }
    frameGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(p), mat));
  };
  ring(RC, 0, rimMat); ring(RC, -DEEP, gridMat);
  // depth axis at the west rim: ticks every 500 dbar, labelled
  const ax = -RC, az = 0;
  const axPts = [new THREE.Vector3(ax, 0, az), new THREE.Vector3(ax, -DEEP, az)];
  for (let p = 0; p <= G.presMax; p += 500) {
    axPts.push(new THREE.Vector3(ax, yOf(p), az), new THREE.Vector3(ax - 25, yOf(p), az));
    depthLabels.push(addLabel(`${p} dbar`, ax - 95, yOf(p), az, 'axis'));
  }
  // three more posts, unlabelled, so the volume reads as a cylinder
  for (const a of [Math.PI / 2, 0, 3 * Math.PI / 2])
    axPts.push(new THREE.Vector3(RC * Math.cos(a), 0, RC * Math.sin(a)),
      new THREE.Vector3(RC * Math.cos(a), -DEEP, RC * Math.sin(a)));
  frameGroup.add(new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(axPts), rimMat));
  // 200 km scale bar on the surface, south side
  const sb = [new THREE.Vector3(-100, 0, RC + 60), new THREE.Vector3(100, 0, RC + 60),
  new THREE.Vector3(-100, 0, RC + 50), new THREE.Vector3(-100, 0, RC + 70),
  new THREE.Vector3(100, 0, RC + 50), new THREE.Vector3(100, 0, RC + 70)];
  frameGroup.add(new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(sb),
    new THREE.LineBasicMaterial({ color: 0x1c2128 })));
  addLabel('200 km', 0, 0, RC + 95, 'axis');
  addLabel('N ↑', 0, 0, -RC - 60, 'axis');
}

// ---------- sea floor (GEBCO_2021, Seabed 2030) ----------
// Context, not a measurement being interpolated: it is a published grid, drawn
// as it is. It exists to say how much of this water column the profiles do not
// reach — they stop at 2000 dbar over a floor that is 4.5-5.5 km down.
// Elevation is metres and the casts are dbar; yOf treats them as the same, an
// error of about 1.5 % at 5 km that no reading of this surface depends on.
// The mesh, not the grid, was the limit: 150 nodes across 1800 km is 12 km, and
// the New England Seamounts are a few km across at the summit, so each came out
// as a single spike and read as noise. 400 nodes is 4.5 km, against 5.5 km data.
const BATHY_N = 400;                               // mesh nodes across the disc
const bathyGroup = new THREE.Group();
world.add(bathyGroup);
// transparent so the toggle can fade it, depthWrite kept because it is a solid
// heightfield and the water column has to sit in front of it
const bathyMat = new THREE.MeshBasicMaterial({
  vertexColors: true, side: THREE.DoubleSide, transparent: true, opacity: 0, depthWrite: true
});
const bathyZ = decodeU16(BATHY.z, ...BATHY.rng);
function bathyAt(lon, lat) {
  const fx = (lon - BATHY.lon0) / (BATHY.lon1 - BATHY.lon0) * (BATHY.nlon - 1);
  const fy = (lat - BATHY.lat0) / (BATHY.lat1 - BATHY.lat0) * (BATHY.nlat - 1);
  const a = Math.max(0, Math.min(BATHY.nlon - 2, Math.floor(fx)));
  const b = Math.max(0, Math.min(BATHY.nlat - 2, Math.floor(fy)));
  const u = Math.max(0, Math.min(1, fx - a)), v = Math.max(0, Math.min(1, fy - b));
  const o = b * BATHY.nlon + a;
  return (bathyZ[o] * (1 - u) + bathyZ[o + 1] * u) * (1 - v)
    + (bathyZ[o + BATHY.nlon] * (1 - u) + bathyZ[o + BATHY.nlon + 1] * u) * v;
}
// Flat colour on a MeshBasicMaterial would hide every ridge, and adding lights
// for one mesh is not worth it: shade from the slope instead, one dot product
// against a fixed sun. The ramp is deliberately outside the magma scale — this
// is not a measured field and must not read as one.
function bathyBuild() {
  for (const ch of [...bathyGroup.children]) { ch.geometry.dispose(); bathyGroup.remove(ch); }
  const n = BATHY_N, step = 2 * RC / (n - 1);
  const pos = new Float32Array(n * n * 3), col = new Float32Array(n * n * 3);
  const el = new Float32Array(n * n), raw = new Float32Array(n * n);
  for (let b = 0; b < n; b++) for (let a = 0; a < n; a++) {
    const i = b * n + a, x = -RC + a * step, z = -RC + b * step;
    const e = bathyAt(lonOfX(x), latOfZ(z));
    raw[i] = e;
    el[i] = Math.min(0, e);                          // land sits at the surface plane
    pos[3 * i] = x; pos[3 * i + 1] = el[i] / 1000; pos[3 * i + 2] = z;
  }
  const SUN = [0.45, 0.78, 0.44];
  for (let b = 0; b < n; b++) for (let a = 0; a < n; a++) {
    const i = b * n + a;
    const ea = el[b * n + Math.min(n - 1, a + 1)] - el[b * n + Math.max(0, a - 1)];
    const eb = el[Math.min(n - 1, b + 1) * n + a] - el[Math.max(0, b - 1) * n + a];
    // normal of the true-scale surface: metres of relief over metres of ground
    const nx = -ea / (2 * step * 1000), nz = -eb / (2 * step * 1000);
    const inv = 1 / Math.hypot(nx, 1, nz);
    const lam = Math.max(0, (nx * SUN[0] + SUN[1] + nz * SUN[2]) * inv);
    const sh = 0.62 + 0.38 * lam;
    // land is a different substance from sea floor and must not shade into the
    // shelf: the NW rim of this disc reaches the American coast
    let r, g, bl;
    if (raw[i] > 0) { r = 0.80; g = 0.75; bl = 0.64; }
    else {
      const t = Math.max(0, Math.min(1, -raw[i] / 6000));   // 0 shore, 1 abyss
      r = 0.72 - 0.44 * t; g = 0.78 - 0.46 * t; bl = 0.82 - 0.42 * t;
    }
    col[3 * i] = r * sh; col[3 * i + 1] = g * sh; col[3 * i + 2] = bl * sh;
  }
  const idx = [];
  const inside = i => Math.hypot(pos[3 * i], pos[3 * i + 2]) <= RC;
  for (let b = 0; b < n - 1; b++) for (let a = 0; a < n - 1; a++) {
    const i = b * n + a, j = i + 1, k = i + n, l = k + 1;
    if (!(inside(i) && inside(j) && inside(k) && inside(l))) continue;
    idx.push(i, k, j, j, k, l);
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  geo.setAttribute('color', new THREE.BufferAttribute(col, 3));
  geo.setIndex(idx);
  bathyGroup.add(new THREE.Mesh(geo, bathyMat));
}
function bathyStats() {
  let lo = 0, sum = 0, n = 0;
  for (let b = 0; b < BATHY.nlat; b++) for (let a = 0; a < BATHY.nlon; a++) {
    const x = (a / (BATHY.nlon - 1)) * (BATHY.lon1 - BATHY.lon0) + BATHY.lon0;
    const y = (b / (BATHY.nlat - 1)) * (BATHY.lat1 - BATHY.lat0) + BATHY.lat0;
    if (Math.hypot(xOf(x), zOf(y)) > RC) continue;
    const e = bathyZ[b * BATHY.nlon + a];
    if (e >= 0) continue;
    lo = Math.min(lo, e); sum += e; n++;
  }
  return { deepest: -lo, mean: -sum / Math.max(1, n) };
}

// The surface graticule and the layer's fishnet are two grids over one scene
// and read as a single confused mesh when both are up. The showcase runs one at
// a time: coordinates while the water column is bare, the fishnet once there is
// a surface whose shape is the thing to read. The rim, depth axis and scale bar
// stay either way — those are what keep the box legible.
function setGraticule(on) {
  if (gratMesh) gratMesh.visible = on;
  for (const L of gratLabels) L.visible = on;
}

// ---------- casts: constant-radius tubes, per-vertex colour = variable ----------

const RINGS = 96, RADIAL = 10;
let TUBE_R = 5.5;                                  // km
const castGroup = new THREE.Group();
world.add(castGroup);
const castMeshes = [];

function tubeGeometry(c) {
  const x = xOf(c.lon), z = zOf(c.lat);
  const p0 = c.pres[0], p1 = c.pres[c.pres.length - 1];
  const n = (RINGS + 1) * (RADIAL + 1);
  const pos = new Float32Array(n * 3), col = new Float32Array(n * 4);
  const val = new Float32Array(RINGS + 1);        // per-ring sample index
  let bin = 0, vi = 0;
  for (let i = 0; i <= RINGS; i++) {
    const p = p0 + (p1 - p0) * i / RINGS;
    while (bin < c.pres.length - 1 && c.pres[bin + 1] < p) bin++;
    val[i] = bin;
    const y = yOf(p);
    for (let j = 0; j <= RADIAL; j++) {
      const a = j / RADIAL * Math.PI * 2;
      pos[vi++] = x + TUBE_R * Math.cos(a); pos[vi++] = y; pos[vi++] = z + TUBE_R * Math.sin(a);
    }
  }
  const idx = [];
  for (let i = 0; i < RINGS; i++) for (let j = 0; j < RADIAL; j++) {
    const a = i * (RADIAL + 1) + j, b = a + RADIAL + 1;
    idx.push(a, b, a + 1, a + 1, b, b + 1);
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  geo.setAttribute('color', new THREE.BufferAttribute(col, 4));
  geo.setIndex(idx);
  geo.userData.bins = val;
  return geo;
}
function colorTube(mesh) {
  const c = casts[mesh.userData.ci];
  const col = mesh.geometry.getAttribute('color');
  const bins = mesh.geometry.userData.bins;
  const key = curVar === 'conf' ? null : curVar;
  let vi = 0;
  const w = TIME.on ? ageWeight(c, 0) : 1;
  const a = TIME.on ? TUBE_AFLOOR + (1 - TUBE_AFLOOR) * Math.pow(supEff(w), 0.75) : 1;
  const hi = isHi(mesh.userData.ci);
  mesh.material = a >= 0.995 ? (hi ? tubeMatHi : tubeMat) : (hi ? tubeMatFadeHi : tubeMatFade);
  for (let i = 0; i <= RINGS; i++) {
    // With the timeline off a measurement is support 1 by definition and gets
    // the continuous ramp. With it on, the tube carries its own temporal
    // support: still a perfect measurement of when it was taken, but a
    // weakening claim on the water at the playhead — so it climbs the same
    // VSUP tree as the interpolated field.
    //
    // ONE weight for the whole tube, at the depth where L_t is resolved.
    // Earlier this was per-ring, ageWeight(c, p), so a cast thinned from the
    // top down. That looked like a measured depth dependence and mostly was
    // not: L_t is fitted only above 150 dbar (6.4 and 6.8 d in the two bands,
    // i.e. no depth dependence within the part that IS resolved), and below
    // that it is a log-ramp to the record length that this file elsewhere
    // calls assumed. So the gradient down the tube was drawing the assumption,
    // not the data. Fading the whole cast at the one timescale the record
    // actually resolves says less and claims only what was measured. It also
    // matches the timeline's per-cast ticks, which have always used this same
    // weight, so a tick and its tube now agree.
    const rgb = key
      ? (TIME.on ? vsupColor(c[key][bins[i]], w) : varColor(c[key][bins[i]]))
      : lutColor(LUT.support, w);
    for (let j = 0; j <= RADIAL; j++) {
      col.array[vi++] = rgb[0]; col.array[vi++] = rgb[1]; col.array[vi++] = rgb[2]; col.array[vi++] = a;
    }
  }
  col.needsUpdate = true;
}
// 0, not a floor. An earlier version held a floor here on the argument that a
// cast stays a perfect measurement of the water at the moment it was taken, so
// it should never disappear. That argument is about the OBSERVATION; what the
// tube colour is drawing is the observation's claim on the water at the
// playhead, and that claim does decay with distance from the cast — in time
// here, since a tube sits at zero distance in space from itself. The "a cast is
// here" statement is carried by the surface marker, which does not fade.
const TUBE_AFLOOR = 0;
// depthWrite stays on: these are solid cylinders, and 83 of them sorted only by
// object order would blend through each other. Rejecting the far ones by depth
// costs the correct blend where two tubes overlap and buys back the geometry.
let hoveredCast = null;
const pins = [];                                   // cast indices, max 3
// Tube against sheet is the one place this scene needs per-pixel transparency
// ordering. Three sorts transparent objects by the depth of the object's
// origin, and every tube and every sheet here has its geometry in world
// coordinates with the mesh at the origin, so they all tie and the order is
// creation order: tubes, then sheets. A tube that writes depth is then
// correctly in front of the sheet where it is nearer; a tube that does not is
// painted over by the sheet everywhere. Two materials switched at a weight of
// 0.995 gave the first to opaque tubes and the second to faded ones, and the
// switch was a visible snap: a tube at 0.99 lost the sheet's occlusion in one
// frame.
//
// With MSAA on, alpha-to-coverage does it per sample instead: the tube covers
// a fraction alpha of each pixel's samples and writes depth for those only, so
// a sheet drawn afterwards fills the rest where the tube is in front and
// blends over all of them where it is behind. Both cases come out as the exact
// blend, for every weight, with one material and no threshold. The cost is
// that alpha is quantised to the sample count (dithered on most GPUs), so a
// faded tube reads as a fine stipple rather than a smooth wash. The material
// is not `transparent`: alpha drives coverage and must not also blend, and the
// opaque list is drawn before every sheet, which is the order required.
// Without MSAA, alpha-to-coverage degrades to a 0.5 cutoff, so the two-material
// scheme is kept as the fallback there.
//
// The trap in that, and it bit: three always asks for a canvas WITH an alpha
// channel, and a material that is not `transparent` has blending switched off,
// so the tube wrote its own alpha into the framebuffer at every sample it
// covered. That is a hole in the canvas's opacity, and the sheet — depth-
// rejected at exactly those samples — could not fill it back in. The browser
// then composited the hole over the page and the page colour came through,
// bleaching a tube-shaped band out of everything already drawn behind it.
// Measured: colour 0.599 at alpha 0.31 over the 0.894 ground resolved to 0.973,
// brighter than the ground it should have darkened, so the pycnocline read as
// missing behind a faded tube and came back only once the tube reached zero and
// stopped covering anything. CustomBlending fixes it without moving the tube
// out of the opaque queue (three forces blending off only for NormalBlending):
// ONE/ZERO on RGB is what no blending did, since coverage has already done the
// mixing, and ZERO/ONE on alpha leaves the canvas opaque. Anything else that
// draws with blending off must keep its alpha at 1 for the same reason.
const MSAA = (() => { const gl = renderer.getContext(); return gl.getParameter(gl.SAMPLES) > 0; })();
const tubeMat = MSAA
  ? new THREE.MeshBasicMaterial({
    vertexColors: true, alphaToCoverage: true, depthWrite: true,
    blending: THREE.CustomBlending,
    blendSrc: THREE.OneFactor, blendDst: THREE.ZeroFactor,          // replace RGB
    blendSrcAlpha: THREE.ZeroFactor, blendDstAlpha: THREE.OneFactor, // keep the canvas opaque
  })
  : new THREE.MeshBasicMaterial({ vertexColors: true, transparent: true, depthWrite: true });
// fallback only: a faded tube is a decayed claim, not a wall, so below full
// weight it stops writing depth and takes the transparent sort
const tubeMatFade = MSAA ? tubeMat
  : new THREE.MeshBasicMaterial({ vertexColors: true, transparent: true, depthWrite: false });
// Hover and pin brighten the cast rather than recolouring it: the ramp is the
// data and must not move. material.color multiplies the vertex colour, and
// THREE.Color does not clamp, so a value above 1 lifts the whole tube toward
// the top of the ramp without changing which end of it a depth sits at.
const tubeMatHi = tubeMat.clone();
const tubeMatFadeHi = MSAA ? tubeMatHi : tubeMatFade.clone();
tubeMatHi.color.setRGB(1.45, 1.45, 1.45);
tubeMatFadeHi.color.setRGB(1.45, 1.45, 1.45);
const isHi = ci => ci === hoveredCast || pins.includes(ci);
// only the material changes, so this never touches the colour attribute
function applyHilite(ci) {
  const m = castMeshes[ci];
  if (!m) return;
  const faded = m.material === tubeMatFade || m.material === tubeMatFadeHi;
  m.material = faded ? (isHi(ci) ? tubeMatFadeHi : tubeMatFade)
    : (isHi(ci) ? tubeMatHi : tubeMat);
}
// surface markers: filled disc = delayed-mode/adjusted, ring = real-time
const markGroup = new THREE.Group();
world.add(markGroup);
const discGeo = new THREE.CircleGeometry(7, 20);
const ringGeo = new THREE.RingGeometry(4.5, 7.5, 20);
const markMat = new THREE.MeshBasicMaterial({ color: 0x1c2128, side: THREE.DoubleSide });
// The lift is divided by the exaggeration, so it stays a fixed distance off the
// surface plane on screen. A constant pre-exag offset is 30 % of the whole
// 2-unit column (DEEP = 2), so the markers floated above the water and climbed
// as the exaggeration went up.
const MARK_LIFT = 0.6;
for (let ci = 0; ci < casts.length; ci++) {
  const c = casts[ci];
  const mesh = new THREE.Mesh(tubeGeometry(c), tubeMat);
  mesh.userData.ci = ci;
  colorTube(mesh);
  castGroup.add(mesh);
  castMeshes.push(mesh);
  const mk = new THREE.Mesh(c.mode === 'R' ? ringGeo : discGeo, markMat);
  mk.rotation.x = -Math.PI / 2;
  mk.position.set(xOf(c.lon), MARK_LIFT, zOf(c.lat));   // setExag divides this by the exaggeration
  mk.userData.ci = ci;
  markGroup.add(mk);
}
function rebuildCasts() {
  for (const m of castMeshes) {
    m.geometry.dispose();
    m.geometry = tubeGeometry(casts[m.userData.ci]);
    colorTube(m);
  }
}
function recolorCasts() { for (const m of castMeshes) colorTube(m); }

// ---------- iso-surfaces: one surface, dragged by value ----------
//
// One surface at a time, on whichever variable is being coloured: an
// isopycnal for sigma0, an isotherm for temp, an isohaline for psal. The
// level is continuous (slider), and the sheet is cut from the live field, so
// it moves as the playhead moves and greys and fades with support like every
// other interpolated layer. Levels are no longer a fixed menu because there
// is nothing special about 27.0 — the useful question is what the surface
// does as you sweep it through the water column.

const ISO = {
  on: true,
  // remembered per variable, so switching colour-by and back keeps the level
  level: { temp: 15, psal: 36.0, sigma0: 27.0 },
};
// support is a property of the sampling, not of the water: an iso-surface of
// it would be circular, so that mode keeps the density surface
const isoVar = () => (curVar === 'conf' ? 'sigma0' : curVar);
const isoLevel = () => ISO.level[isoVar()];
const ISO_NAME = { temp: ['Isotherm', '°C', 'T'], psal: ['Isohaline', 'PSS-78', 'S'], sigma0: ['Isopycnal', 'kg m⁻³', 'σ₀'] };
const isoDec = () => (isoVar() === 'temp' ? 1 : 2);
const isoText = () => `${ISO_NAME[isoVar()][2]} ${isoLevel().toFixed(isoDec())}`;

// Mesh for the iso-surface and the layer: two nodes per grid cell by default.
// ?fine=k multiplies that — for a still, where the silhouette of a sheet is
// read at print size and the frame rate does not matter. Read here rather than
// in the URL block at the bottom because every array below is sized from it.
const FINE = Math.max(1, Math.min(4, +(new URLSearchParams(location.search).get('fine') || 1)));
const FX = FINE * 2 * (NX - 1) + 1, FY = FINE * 2 * (NY - 1) + 1;
const isoGroup = new THREE.Group();
world.add(isoGroup);
const isoLabels = [];

// deepest crossing of `level` in one grid column, scanned bottom-up: the
// stable position, where a top-down scan picks up spurious mixed-layer
// crossings next to outcrop holes. Read off the grid nodes and interpolated
// afterwards - the field is already smooth on the 25 km kernel, so sampling
// it at half-cells before the crossing only costs time.
function isoSampleNode(ix, iy, level, arr) {
  const o = iy * NX + ix, plane = NY * NX;
  let s1 = arr[(NZ - 1) * plane + o];
  for (let k = NZ - 2; k >= 0; k--) {
    const s0 = arr[k * plane + o];
    if ((s0 - level) * (s1 - level) <= 0 && s0 !== s1) {
      const f = (level - s0) / (s1 - s0);
      const c0 = F.conf[k * plane + o], c1 = F.conf[(k + 1) * plane + o];
      return { pres: gzAt(k) + f * (gzAt(k + 1) - gzAt(k)), infl: c0 + f * (c1 - c0) };
    }
    s1 = s0;
  }
  return null;
}
const isoMat = new THREE.ShaderMaterial({
  vertexShader: `
    attribute vec3 col; attribute float alpha; varying vec3 vC; varying vec3 vN; varying float vA;
    void main(){ vC = col; vA = alpha; vN = normalize(normalMatrix * normal);
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position,1.0); }`,
  fragmentShader: `
    varying vec3 vC; varying vec3 vN; varying float vA;
    void main(){ float l = 0.7 + 0.4 * abs(vN.z);       // headlight lambert, two-sided
      gl_FragColor = vec4(vC * l, vA); }`,
  transparent: true, depthWrite: true,
  side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1,
});
// the slab has two faces the eye must see through, so it does not write depth
const slabMat = isoMat.clone();
slabMat.depthWrite = false;
// the slab's vertical walls are wound outward and drawn single-sided: the far
// wall of a closed volume is culled instead of being laid over the near one,
// which is most of the fog you get from a two-sided box
const wallMat = slabMat.clone();
wallMat.side = THREE.FrontSide;
const ISO_NEUTRAL = [0.66, 0.69, 0.73];
// sheet colour: the level's own colour on the same ramp as the field, pulled
// toward neutral by the same VSUP layer rule as the slices, so "greyer = less
// support" is one rule everywhere
// five colours per (variable, level), shared arrays like vsupLeaf's
const isoColMemo = { key: null, base: null, tab: [] };
const isoVertexColor = (sup) => {
  const key = `${isoVar()}|${isoLevel()}`;
  if (isoColMemo.key !== key) { isoColMemo.key = key; isoColMemo.base = varColor(isoLevel(), isoVar()); isoColMemo.tab.length = 0; }
  const layer = vsupLayer(sup), m = layer / (VSUP_LAYERS - 1);
  return isoColMemo.tab[layer] ??= isoColMemo.base.map((c, i) => c + (ISO_NEUTRAL[i] - c) * m);
};
const isoLineMat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.6 });
// A fishnet over a heightSurface: same vertices, every `step` row and column,
// dropping any segment that would bridge a slope the surface does not have.
// Colour blends the ground toward ink by the vertex's own alpha, so the net
// carries the support encoding rather than overriding it.
// `gain` rescales that alpha before it is used as the blend: a sheet painted at
// a low peak opacity (the pycnocline tops out at 0.34, so its net came out a
// pale grey barely off the ground colour) still wants a net that reaches full
// ink where support is full. Gain is 1/peak-alpha, so the blend runs 0..1 over
// the sheet's actual range instead of over 0..0.34 of it.
// scratch for fishnet: a node has at most one segment to the right and one
// down, so 2 * FX * FY segments of two vertices bounds it
const FNET = { p: new Float32Array(FX * FY * 12), c: new Float32Array(FX * FY * 12) };
function fishnet(S, step, slopeMax, { gain = 1, ink = 0.11, mat = isoLineMat } = {}) {
  const { pos, alp, vid } = S;
  const lp = FNET.p, lc = FNET.c;
  const [g0, g1, g2] = GROUND_RGB;
  let n = 0;
  const put = a => {
    const f = clamp01(alp[a] * gain), o = 3 * n++;
    lp[o] = pos[3 * a]; lp[o + 1] = pos[3 * a + 1]; lp[o + 2] = pos[3 * a + 2];
    lc[o] = g0 + (ink - g0) * f; lc[o + 1] = g1 + (ink - g1) * f; lc[o + 2] = g2 + (ink - g2) * f;
  };
  const sm2 = slopeMax * slopeMax;
  const seg = (a, b) => {
    if (a < 0 || b < 0) return;
    const dx = pos[3 * a] - pos[3 * b], dy = pos[3 * a + 1] - pos[3 * b + 1], dz = pos[3 * a + 2] - pos[3 * b + 2];
    if (dy * dy > sm2 * (dx * dx + dz * dz)) return;
    put(a); put(b);
  };
  for (let fy = 0; fy < FY; fy += step) for (let fx = 0; fx < FX - 1; fx++) seg(vid[fy * FX + fx], vid[fy * FX + fx + 1]);
  for (let fx = 0; fx < FX; fx += step) for (let fy = 0; fy < FY - 1; fy++) seg(vid[fy * FX + fx], vid[(fy + 1) * FX + fx]);
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(lp.slice(0, 3 * n), 3));
  g.setAttribute('color', new THREE.BufferAttribute(lc.slice(0, 3 * n), 3));
  return new THREE.LineSegments(g, mat);
}
const GROUND_RGB = [0.894, 0.906, 0.922];
const ISO_STEP = 3;                                  // fishnet spacing, fine cells
const SLOPE_MAX = 0.008;                             // km/km; real fronts here are <0.006
const px = fx => xOf(G.lon0 + (G.lon1 - G.lon0) * fx / (FX - 1));
const pz = fy => zOf(G.lat0 + (G.lat1 - G.lat0) * fy / (FY - 1));

// 3x3 median (>=5 valid neighbours) kills single-cell spikes and ragged rims
function median3x3(src) {
  const out = new Float32Array(FX * FY).fill(NaN);
  const nb = new Float32Array(9);                    // insertion sort in place: no per-node allocation
  for (let fy = 0; fy < FY; fy++) for (let fx = 0; fx < FX; fx++) {
    let n = 0;
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      const yy = fy + dy, xx = fx + dx;
      if (yy < 0 || yy >= FY || xx < 0 || xx >= FX) continue;
      const v = src[yy * FX + xx];
      if (Number.isNaN(v)) continue;
      let k = n++;
      while (k > 0 && nb[k - 1] > v) { nb[k] = nb[k - 1]; k--; }
      nb[k] = v;
    }
    if (n >= 5) out[fy * FX + fx] = nb[n >> 1];
  }
  return out;
}

// coarse grid node -> fine mesh node. FX-1 = 2(NX-1), so this is a plain
// half-cell bilinear; a node is undefined unless every corner it needs is.
// The iso-surface can be expanded this way because its depth is a continuous
// function of the field; a gradient layer's outline is not (see below).
function expandFine(c) {
  const out = new Float32Array(FX * FY).fill(NaN);
  for (let fy = 0; fy < FY; fy++) {
    const gy = fy / 2, iy = Math.min(fy >> 1, NY - 2), wy = gy - iy;
    for (let fx = 0; fx < FX; fx++) {
      const gx = fx / 2, ix = Math.min(fx >> 1, NX - 2), wx = gx - ix;
      const w = [(1 - wx) * (1 - wy), wx * (1 - wy), (1 - wx) * wy, wx * wy];
      const v = [c[iy * NX + ix], c[iy * NX + ix + 1], c[(iy + 1) * NX + ix], c[(iy + 1) * NX + ix + 1]];
      let acc = 0, bad = false;
      for (let q = 0; q < 4; q++) {
        if (w[q] <= 0) continue;
        if (Number.isNaN(v[q])) { bad = true; break; }
        acc += v[q] * w[q];
      }
      if (!bad) out[fy * FX + fx] = acc;
    }
  }
  return out;
}

// a heightfield over the fine grid -> triangles, per-vertex colour and alpha.
// Returns the vertex ids so a second surface can be stitched to this one.
// opt.xs / opt.zs override the regular node positions (the slab moves its rim
// nodes onto the contour where the layer ends); opt.soft marks those moved
// nodes, whose edges skip the slope guard — a rim that tapers to nothing is
// steep by construction and the guard would delete exactly it.
// scratch for heightSurface, sized for the whole fine grid; the geometry gets
// a copy of the used prefix, so two surfaces built in a row (the layer's top
// and base) do not share a buffer
const HS = {
  pos: new Float32Array(FX * FY * 3), col: new Float32Array(FX * FY * 3), alp: new Float32Array(FX * FY),
  soft: new Uint8Array(FX * FY), idx: new Uint32Array((FX - 1) * (FY - 1) * 6), vid: new Int32Array(FX * FY),
};
function heightSurface(depth, sup, colFn, alpFn, slopeMax, opt = {}) {
  const { xs = null, zs = null, soft = null, box = null } = opt;
  const x0 = box ? box.x0 : 0, x1 = box ? box.x1 : FX - 1;
  const y0 = box ? box.y0 : 0, y1 = box ? box.y1 : FY - 1;
  const { pos, col, alp, idx, vid } = HS, vsoft = HS.soft;
  vid.fill(-1);
  let nv = 0, ni = 0;
  for (let fy = y0; fy <= y1; fy++) for (let fx = x0; fx <= x1; fx++) {
    const i = fy * FX + fx;
    const d = depth[i];
    if (Number.isNaN(d)) continue;
    vid[i] = nv;
    const s = sup[i], o = 3 * nv;
    pos[o] = xs ? xs[i] : px(fx); pos[o + 1] = yOf(d); pos[o + 2] = zs ? zs[i] : pz(fy);
    const rgb = colFn(s, d, i);                      // d, i: for a sheet coloured by its own depth or by a per-node field
    col[o] = rgb[0]; col[o + 1] = rgb[1]; col[o + 2] = rgb[2];
    alp[nv] = alpFn(s);
    vsoft[nv] = soft ? soft[i] : 0;
    nv++;
  }
  const sm2 = slopeMax * slopeMax;
  const okTri = (a, b, c) => {
    const s = (i, j) => {
      if (vsoft[i] || vsoft[j]) return true;
      const dx = pos[3 * i] - pos[3 * j], dy = pos[3 * i + 1] - pos[3 * j + 1], dz = pos[3 * i + 2] - pos[3 * j + 2];
      return dy * dy <= sm2 * (dx * dx + dz * dz);
    };
    return s(a, b) && s(b, c) && s(a, c);
  };
  for (let fy = y0; fy < y1; fy++) for (let fx = x0; fx < x1; fx++) {
    const a = vid[fy * FX + fx], b = vid[fy * FX + fx + 1],
      c = vid[(fy + 1) * FX + fx], d = vid[(fy + 1) * FX + fx + 1];
    if (a >= 0 && b >= 0 && c >= 0 && okTri(a, b, c)) { idx[ni++] = a; idx[ni++] = c; idx[ni++] = b; }
    if (b >= 0 && c >= 0 && d >= 0 && okTri(b, c, d)) { idx[ni++] = b; idx[ni++] = c; idx[ni++] = d; }
  }
  const P = pos.slice(0, 3 * nv), A = alp.slice(0, nv), I = idx.slice(0, ni);
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(P, 3));
  geo.setAttribute('col', new THREE.BufferAttribute(col.slice(0, 3 * nv), 3));
  geo.setAttribute('alpha', new THREE.BufferAttribute(A, 1));
  geo.setAttribute('normal', new THREE.BufferAttribute(vertexNormals(P, I, nv), 3));
  geo.setIndex(new THREE.BufferAttribute(I, 1));
  return { geo, pos: P, alp: A, vid: vid.slice() };
}
// Area-weighted vertex normals over the flat arrays: the same sum of face
// normals computeVertexNormals forms, without a Vector3 per corner. Three's
// version was half of heightSurface's time across the three sheets a tick
// builds. Face normal (C - B) x (A - B), as three winds it.
const VN = new Float32Array(FX * FY * 3);
function vertexNormals(P, I, nv) {
  const n = VN.subarray(0, 3 * nv);
  n.fill(0);
  for (let k = 0; k < I.length; k += 3) {
    const a = 3 * I[k], b = 3 * I[k + 1], c = 3 * I[k + 2];
    const cbx = P[c] - P[b], cby = P[c + 1] - P[b + 1], cbz = P[c + 2] - P[b + 2];
    const abx = P[a] - P[b], aby = P[a + 1] - P[b + 1], abz = P[a + 2] - P[b + 2];
    const nx = cby * abz - cbz * aby, ny = cbz * abx - cbx * abz, nz = cbx * aby - cby * abx;
    n[a] += nx; n[a + 1] += ny; n[a + 2] += nz;
    n[b] += nx; n[b + 1] += ny; n[b + 2] += nz;
    n[c] += nx; n[c + 1] += ny; n[c + 2] += nz;
  }
  for (let v = 0; v < 3 * nv; v += 3) {
    const l = Math.hypot(n[v], n[v + 1], n[v + 2]) || 1;
    n[v] /= l; n[v + 1] /= l; n[v + 2] /= l;
  }
  return n.slice();
}

// measured crossings: a short collar on each tube where the binned profile
// crosses the level (first crossing from the top). These are data; the sheet
// between them is interpolation.
const isoRingGroup = new THREE.Group();
world.add(isoRingGroup);
// one geometry holding a copy of `proto` at every (x, y, z) in `places`: the
// collars and the iso rings are 83-162 identical open cylinders each, and as
// separate meshes they were 83-162 draw calls for a few hundred triangles.
// Nothing picks them (pickAt reads the tubes and the planes), so nothing
// needs them apart. Position only: their materials are unlit.
function mergedCylinders(proto, places) {
  const pp = proto.getAttribute('position').array, pi = proto.getIndex().array;
  const nv = pp.length / 3, n = places.length / 3;
  const pos = new Float32Array(pp.length * n), idx = new Uint32Array(pi.length * n);
  for (let m = 0; m < n; m++) {
    const x = places[3 * m], y = places[3 * m + 1], z = places[3 * m + 2], vo = m * nv;
    for (let v = 0; v < nv; v++) {
      pos[(vo + v) * 3] = pp[3 * v] + x; pos[(vo + v) * 3 + 1] = pp[3 * v + 1] + y; pos[(vo + v) * 3 + 2] = pp[3 * v + 2] + z;
    }
    for (let k = 0; k < pi.length; k++) idx[m * pi.length + k] = pi[k] + vo;
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setIndex(new THREE.BufferAttribute(idx, 1));
  return g;
}
function buildIsoRings() {
  for (const ch of [...isoRingGroup.children]) { ch.geometry.dispose(); isoRingGroup.remove(ch); }
  if (!ISO.on) return;
  const key = isoVar(), level = isoLevel();
  const mat = new THREE.MeshBasicMaterial({ color: new THREE.Color(...varColor(level, key)) });
  const geo = new THREE.CylinderGeometry(TUBE_R * 2.2, TUBE_R * 2.2, 0.03, 14, 1, true);
  // castVisible is constant (castExists is () => true); if it ever becomes a
  // filter, syncIsoRings has to rebuild here rather than toggle
  const places = [];
  casts.forEach((c) => {
    if (!castVisible(c)) return;
    for (let i = 1; i < c.pres.length; i++) {
      if ((c[key][i - 1] - level) * (c[key][i] - level) > 0) continue;
      const f = (level - c[key][i - 1]) / (c[key][i] - c[key][i - 1] || 1);
      places.push(xOf(c.lon), yOf(c.pres[i - 1] + f * (c.pres[i] - c.pres[i - 1])), zOf(c.lat));
      break;
    }
  });
  if (places.length) isoRingGroup.add(new THREE.Mesh(mergedCylinders(geo, places), mat));
  geo.dispose();
  syncIsoRings();
}
function syncIsoRings() {
  isoRingGroup.visible = ISO.on && $('tg-casts').checked;
}
function buildIsoSheets() {
  for (const ch of [...isoGroup.children]) { ch.geometry.dispose(); isoGroup.remove(ch); }
  for (const L of isoLabels) { L.el.remove(); labels.splice(labels.indexOf(L), 1); }
  isoLabels.length = 0;
  if (!ISO.on) return;
  const arr = F[isoVar()], level = isoLevel();
  const cDepth = new Float32Array(NX * NY).fill(NaN);
  const cInfl = new Float32Array(NX * NY).fill(NaN);
  for (let iy = 0; iy < NY; iy++) for (let ix = 0; ix < NX; ix++) {
    if (Math.hypot(px(2 * ix), pz(2 * iy)) > RC) continue;
    const s = isoSampleNode(ix, iy, level, arr);
    if (!s) continue;
    cDepth[iy * NX + ix] = s.pres; cInfl[iy * NX + ix] = s.infl;
  }
  const inflE = expandFine(cInfl);
  const infl = new Float32Array(FX * FY);
  for (let i = 0; i < infl.length; i++) infl[i] = Number.isNaN(inflE[i]) ? 0 : inflE[i];
  const { geo, pos, alp, vid } = heightSurface(median3x3(expandFine(cDepth)), infl, isoVertexColor, supportAlpha, SLOPE_MAX);
  isoGroup.add(new THREE.Mesh(geo, isoMat));
  isoGroup.add(fishnet({ pos, alp, vid }, ISO_STEP, SLOPE_MAX));
  // the sheets all sit at the origin, so among themselves three's transparent
  // sort falls back to object id, i.e. whichever was rebuilt last. Pinned to
  // the order the per-tick rebuild used to give: iso-surface, then the layer.
  for (const m of isoGroup.children) m.renderOrder = 1;
  // label the shallowest well-supported point on the surface
  let best = null;
  for (let i = 0; i < FX * FY; i++) {
    const v = vid[i];
    if (v < 0 || infl[i] <= 0.6) continue;
    const y = pos[3 * v + 1];
    if (!best || y > best.y) best = { y, x: pos[3 * v], z: pos[3 * v + 2] };
  }
  if (best) {
    const L = addLabel(isoText(), best.x, best.y + 0.02, best.z, 'iso');
    L.el.style.color = css(varColor(level, isoVar()));
    isoLabels.push(L);
  }
}

// ---------- the permanent pycnocline: measured at each cast, mapped between them ----------
//
// WHICH thermocline. These profiles have two. The seasonal one is shallow and
// sharp — top ~8 dbar, base ~98 dbar in June, peak |dT/dp| 0.076 degC/dbar —
// and it carries 22 % of the temperature drop between the surface and 1000
// dbar. Below it, in 91 of 100 casts, sits the permanent (main) thermocline:
// 528 dbar to 978 dbar, a third of the peak gradient but five times as thick,
// and 59 % of that drop. An earlier version of this layer took the strongest
// gradient run in the upper 300 dbar, so it could only ever return the
// seasonal one; this one goes after the permanent pycnocline, which is the
// structure the profile plots are dominated by.
//
// HOW it is found. Not by a gradient threshold — the two features would then
// have to be separated by a number chosen to separate them, and the deep one
// is the weaker of the two per dbar. Instead this follows OAC-P, the objective
// algorithm of Feucher et al. (JTECH 2016; JGR-Oceans 2019, doi
// 10.1029/2018JC014526), which works on stratification and on the ORDER of its
// extrema:
//   - the seasonal pycnocline is the N^2 maximum in the top 200 dbar;
//   - below it, above 500 dbar, the subtropical mode water is an N^2 MINIMUM
//     (a thick, weakly stratified lens);
//   - the permanent pycnocline is the N^2 maximum below that minimum.
// No threshold decides which is which — their order does, so the detector
// cannot slide from one feature to the other. N^2 comes from TEOS-10 and is
// smoothed over 50 dbar, the scale Feucher et al. use for the North Atlantic
// (they reduce it to 20 m elsewhere; the North Atlantic pycnocline is the
// deepest and thickest of the five subtropical gyres).
//
// The one free parameter is where the layer ENDS, and that is a choice, hence
// the slider: the top and base are where the smoothed N^2 falls to a fraction
// of its peak, walking up and down from the core. Half-max is the default. The
// upward walk is stopped at the mode-water minimum, so the seasonal pycnocline
// can never be absorbed into the layer no matter how low the fraction goes.
//
// Sanity check against the literature: over these 83 casts the core sits at
// 820 dbar with N^2 = 2.5e-5 s^-2, and Feucher et al. (2019, their Table 1)
// give 722 m and 2.3e-5 s^-2 for the western North Atlantic subtropical gyre.
// At the default fraction the layer is 625 dbar thick against their 439 m,
// which is a difference in where the edge is declared, not in what was found —
// their half-Gaussian fit and this half-max walk agree at a fraction of 0.75.
//
// WHY it is not cut from the interpolated field, like the slice and the
// section are. That field is a normalised kernel average with the Gaussian cut
// at 3 L_h = 75 km, and half of the grid nodes that have any data at all
// (748 of 1524) are reached by exactly ONE cast. Where one observation
// dominates, a normalised weighted mean *is* that observation: at those nodes
// the field equals that cast's own profile to within 0.35 degC in 94% of
// cases. So the field is a flat disc of radius 75 km around every isolated
// cast, with all the variation crushed into the narrow overlaps — plateau,
// cliff, plateau. Any layer read off it inherits that shape.
//
// So the layer is taken where it is actually measured — in each cast's own
// 5-dbar profile — and the two depths are then mapped between casts by
// objective analysis (Gauss-Markov / optimal interpolation, the standard
// method for scattered hydrographic data). Unlike a kernel average, OI has no
// compact support and no plateau: away from a cast the estimate relaxes toward
// the mean rather than holding that cast's value. It also returns the thing
// the kernel-weight proxy only gestures at — a real normalised error variance,
// used here as the layer's support.
//
// WHAT THE MAP THEN SHOWS, and it is the opposite of the seasonal layer's
// result. Hold one cast out, map the other 80, predict it: the base is missed
// by 109 dbar against 157 for guessing the basin mean and the top by 67
// against 96 — about half the variance of each, against 12 % for the seasonal
// layer under the same test. So this volume is drawn bowed, and it is entitled
// to be: down past 1200 dbar in the southwest, up toward 800 in the northeast.
// The two surfaces still get their own covariances and their own error
// variances — the base is coherent over ~350 km and the top over ~220 — and
// each sheet is shaded by its own.

// kind: which layer is drawn — 'permanent' is the main pycnocline found on N^2
// (see below), 'seasonal' is the steep thermocline in the top few hundred dbar,
// found on the temperature gradient itself. render: 'slab' is the two-sheet
// volume, 'sheet' one opaque surface through the layer's middle coloured by its
// thickness, 'cloud' the field inside the layer as points.
const THERMO = {
  on: false, frac: 0.5, lo: 0, hi: 1, kind: 'permanent', render: 'slab',
  mask: true, colorBy: 'depth', walls: true
};
let LTH = new Float32Array(0);                       // layer thickness per fine node, for colorBy 'thick'
const THERMO_RGB = [0.26, 0.28, 0.47];
const LAY = {
  seasMax: 200,                                      // dbar, band the seasonal N^2 peak is sought in
  mwMax: 500,                                        // dbar, deepest a mode-water N^2 minimum may sit
  rise: 4,                                           // bins of sustained N^2 rise that end the search for it
  floor: 0.02,                                       // skill below which it is not drawn
  eps: 0.002,                                        // and none at all below this
  nugMin: 0.02,                                      // the nugget at uncertainty weight 0
  alphaFloor: 0,                                     // no floor: it fades to nothing (see below)
  slope: 0.02,                                       // km/km, mesh slope guard
  net: 4,                                            // fishnet spacing, fine cells
  aPeak: 0.34,                                       // peak sheet opacity; the net's gain
  gn: 25,                                            // OI evaluation grid (gn x gn over the box)
  dtFit: 5,                                          // days; pairs this close are "synoptic" for the fit
};
// darker and more opaque than isoLineMat: this net is drawn over a sheet that
// is itself only a third opaque, sitting on a light ground, so at the iso
// material's 0.6 it disappeared into the background.
const netMat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.85 });
const thermoGroup = new THREE.Group();
world.add(thermoGroup);
const thermoLabels = [];

// Sheet colour is the mapped quantity itself — the depth of the surface — on
// the depth ramp, climbing the same VSUP tree as everything else: fewer depth
// bins and closer to neutral as the analysis's own skill falls. It was one flat
// indigo for both sheets, which drew the layer's existence but not its shape;
// the two sheets are 400-800 dbar apart and the base bows ~300 dbar across the
// box, and none of that was in the colour. The ramp spans the drawn depths of
// both sheets together (THERMO.lo/hi, set in buildThermo), so top and base sit
// at opposite ends of it and the thickness reads as the gap between them.
// colorBy 'thick' is the layer's own thickness at that column, on the magma
// ramp and the slate neutral — the same tree, and the same palette, as the VSUP
// figure the proposal opens with. 'depth' is each sheet's own depth, which is
// what the two-sheet pycnocline volume used before there was a thickness to
// draw. Either way the value climbs the tree: fewer bins, then neutral, as the
// cross-validated skill falls, and the opacity fades on top of that.
function thermoColor(sup, d, i) {
  // 'temp': temperature as a monotone function of depth — the magma ramp read
  // the way an oceanographer will read it anyway, warm and bright at the top of
  // the layer, cold and dark at its base. In the synthetic scene this is as
  // invented as the surfaces; on measured data it is a depth proxy, not the
  // field's own temperature, and a caption has to say which it is.
  if (THERMO.colorBy === 'temp')
    return vsupOn(LUT.inferno, 1 - (d - THERMO.lo) / Math.max(1e-6, THERMO.hi - THERMO.lo), sup);
  const u = THERMO.colorBy === 'thick'
    ? (LTH[i] - THERMO.lo) / Math.max(1e-6, THERMO.hi - THERMO.lo)
    : (d - THERMO.lo) / Math.max(1e-6, THERMO.hi - THERMO.lo);
  return THERMO.colorBy === 'thick'
    ? vsupOn(LUT.inferno, u, sup)
    : vsupOn(LUT.depth, u, sup);
}

// --- 1. the layer in one cast, from that cast's own binned N^2 profile ---
// seasonal peak -> mode-water minimum -> permanent pycnocline peak, in that
// order, then walk out to f x peak on each side. Ordering does the
// identification; f only sets the edge. Returns null if the profile has no
// mode-water minimum above 500 dbar or nothing below it, which no cast in this
// month does.
function castLayer(c, f) {
  const p = c.pres, n2 = c.n2, n = p.length;
  let seas = -1;
  for (let i = 0; i < n && p[i] < LAY.seasMax; i++) if (seas < 0 || n2[i] > n2[seas]) seas = i;
  if (seas < 0) return null;
  // The mode-water minimum is the FIRST minimum below the seasonal peak, not
  // the smallest N^2 in some window: in six of these casts the pycnocline core
  // itself sits above 500 dbar, and a window argmin walks straight past the
  // minimum and down the pycnocline's far flank to the edge of the window.
  // "First" needs a persistence guard or noise ends the search early —
  // LAY.rise bins of sustained rise, 20 dbar, well inside the 50 dbar the
  // profile was already smoothed over. The window is then only an outer bound.
  let mw = seas, run = 0;
  for (let i = seas + 1; i < n && p[i] < LAY.mwMax; i++) {
    if (n2[i] < n2[mw]) { mw = i; run = 0; } else if (++run >= LAY.rise) break;
  }
  let pk = -1;
  for (let i = mw + 1; i < n; i++) if (pk < 0 || n2[i] > n2[pk]) pk = i;
  // No minimum and no maximum below it: N^2 just falls off monotonically. That
  // is not a detector failure, it is a profile with no subtropical mode water
  // and so no permanent pycnocline to find — two casts here, both north of
  // 39 degN on the cold side of the Gulf Stream. They are dropped, not forced.
  if (pk < 0 || pk - mw < 2) return null;
  const th = f * n2[pk];
  let a = pk; while (a > mw + 1 && n2[a - 1] > th) a--;
  let b = pk; while (b < n - 1 && n2[b + 1] > th) b++;
  // the base can run past the end of the profile: these floats stop near 2000
  // dbar, and a cast whose N^2 is still above the edge there has had its layer
  // cut by the instrument, not by the ocean. Flagged so the caption can say so.
  return { top: p[a], bot: p[b], core: p[pk], pk: n2[pk], mw: p[mw], trunc: b === n - 1 };
}

// The SEASONAL thermocline, for when that is what is being asked for: the steep
// layer in the top few hundred dbar that caps the mixed layer, and the one with
// biological consequences — light above it, nutrients below, and a gradient
// sharp enough to sit across an organism's tolerance. It is found on dT/dp and
// not on N^2: the word is thermocline, and temperature is what it names. The
// gradient is smoothed over SEAS.smooth bins either side (5-dbar bins, so ~25
// dbar), the core is the largest |dT/dp| above SEAS.maxP, and the edges are
// where the smoothed gradient falls to the same fraction of that peak the other
// detector uses.
//
// WHAT THIS SURFACE IS AND IS NOT. The permanent pycnocline's depth is a field
// this array can map — hold a cast out and the map predicts it to about half
// its variance. The seasonal layer is NOT: the same test on it resolves about
// 12 %, so its shape between casts is very largely the covariance model's
// assumption and not a measurement. That is why the shading still runs and why
// the drawn region is still cut by cross-validated skill. Any caption for a
// figure of this surface has to say so.
const SEAS = { maxP: 300, smooth: 2 };
function castThermo(c, f) {
  const p = c.pres, T = c.temp, n = p.length;
  if (n < 8) return null;
  const g = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const a = Math.max(0, i - SEAS.smooth), b = Math.min(n - 1, i + SEAS.smooth);
    g[i] = p[b] > p[a] ? Math.abs((T[b] - T[a]) / (p[b] - p[a])) : 0;
  }
  let pk = -1;
  for (let i = 0; i < n && p[i] <= SEAS.maxP; i++) if (pk < 0 || g[i] > g[pk]) pk = i;
  if (pk < 0 || !(g[pk] > 0)) return null;
  const edge = f * g[pk];
  let a = pk; while (a > 0 && g[a - 1] >= edge) a--;
  let b = pk; while (b < n - 1 && g[b + 1] >= edge) b++;
  if (!(p[b] > p[a])) return null;
  return { top: p[a], bot: p[b], core: p[pk], pk: g[pk], trunc: b === n - 1 };
}
const castLayerOf = (c, f) => (THERMO.kind === 'seasonal' ? castThermo(c, f) : castLayer(c, f));

// --- 2. covariance chosen by cross-validation over the casts ---
// Not by fitting a variogram. That was tried, and it is a trap here: the pair
// counts run ~30:1 in favour of the far separation bins, which carry no
// information about the nugget or the correlation length, and the resulting
// covariance (L = 240 km, nugget 0.41) predicts a held-out cast WORSE than the
// basin mean does — it explained -5 % of the variance out of sample while
// looking smooth and confident over 93 % of the box. So the two parameters are
// picked by leave-one-out cross-validation instead, using Dubrule's shortcut:
// the LOO residual is (A^-1 d)_i / (A^-1)_ii, so each candidate costs one
// factorisation rather than N of them. (The shortcut holds the mean at its
// full-sample value, so it differs from brute-force LOO by <0.5 % here — fine
// for choosing between models.)
// The base is correlated over hundreds of km, the top over ~100, so the grid
// has to span both. The nugget grid stops at 0.05 rather than running down to
// LAY.nugMin: 0 is the interpolant endpoint the uncertainty slider drives to,
// and the cross-validation surface is a flat ridge down there anyway (the base
// scores 63.6 % at L = 350, nugget 0.02 against 62.5 % at L = 220, nugget 0.1),
// so letting it pick the endpoint would cost nothing in skill and would leave
// the slider with no geometry to move.
const COV_L = [45, 70, 100, 150, 220, 350, 500];
const COV_NUG = [0.05, 0.1, 0.2, 0.3, 0.45, 0.6];
function looRms(X, Z, d, N, L, nug) {
  const rho = (dx, dz) => (1 - nug) * Math.exp(-(dx * dx + dz * dz) / (2 * L * L));
  const A = [];
  for (let i = 0; i < N; i++) {
    A.push(new Float64Array(N));
    for (let j = 0; j < N; j++) A[i][j] = rho(X[i] - X[j], Z[i] - Z[j]) + (i === j ? nug : 0);
  }
  const M = cholesky(A, N);
  if (!M) return Infinity;
  const mu = d.reduce((a2, b2) => a2 + b2, 0) / N;
  const w = cholSolve(M, Float64Array.from(d.map(v => v - mu)), N);
  const e = new Float64Array(N);
  let acc = 0, vacc = 0;
  for (let i = 0; i < N; i++) {
    e.fill(0); e[i] = 1;
    const col = cholSolve(M, e, N);
    const r = w[i] / col[i];
    acc += r * r;
    // The matching half of the same identity: the model's own leave-one-out
    // error variance at cast i is 1/(A^-1)_ii, in the normalised units A is
    // written in (A_ii = (1-nug) + nug = 1). Free here — col is already solved.
    // Comparing it with the residual above is what calibrates the error map.
    vacc += 1 / col[i];
  }
  return { rms: Math.sqrt(acc / N), predVar: vacc / N };
}
// The two surfaces have different covariances and must not share one. The base
// is coherent over ~220 km and the top over ~100 with a much larger nugget;
// forcing the top onto the base's covariance costs it about ten points of
// cross-validated skill and, worse, draws it as if it were as well determined.
// So each is fitted, mapped and shaded on its own.
function fitLayerCov(sel) {
  const idx = sel.map((r, i) => (r ? i : -1)).filter(i => i >= 0);
  const N = idx.length;
  const X = idx.map(i => xOf(casts[i].lon)), Z = idx.map(i => zOf(casts[i].lat));
  const one = key => {
    const d = idx.map(i => sel[i][key]);
    const mu = d.reduce((a, b) => a + b, 0) / N;
    const sd = Math.sqrt(d.reduce((a, b) => a + (b - mu) * (b - mu), 0) / N);
    let best = { L: 220, nug: 0.2, loo: Infinity, predVar: 1 };
    for (const L of COV_L) for (const nug of COV_NUG) {
      const r = looRms(X, Z, d, N, L, nug);
      if (r.rms < best.loo) best = { L, nug, loo: r.rms, predVar: r.predVar };
    }
    // how much of the cast-to-cast variance the array actually resolves. The
    // reference is the basin mean: predicting that alone has RMS = sd.
    const r2 = 1 - (best.loo / sd) * (best.loo / sd);
    // CALIBRATION. c'A^-1 c answers "how tightly does this covariance model pin
    // the field here", which is not "how wrong is this map". Cross-validation
    // picks (L, nugget) by prediction error but nothing then makes the error
    // FIELD reproduce that error: at the casts the base's map claims 41 dbar
    // and actually misses 109, and the top claims 37 and misses 67. Both halves
    // of that comparison are already computed above, in the same units, so the
    // ratio scales the error variance everywhere: cal = measured LOO variance /
    // the model's own. 7.2 on the base, 3.3 on the top. It assumes the SHAPE of
    // the error field is right and only its scale is wrong, and it is anchored
    // at the casts, which is where support is highest — so if anything it still
    // flatters the gaps between them.
    const cal = (best.loo * best.loo) / (sd * sd) / best.predVar;
    return {
      L: best.L, nug: best.nug, sd, loo: best.loo, r2, mu, n: N,
      cal, predSd: Math.sqrt(best.predVar) * sd
    };
  };
  return { top: one('top'), bot: one('bot') };
}

// --- 3. objective analysis over the cast positions ---
function cholesky(A, n) {
  const M = [];
  for (let i = 0; i < n; i++) M.push(new Float64Array(n));
  for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
    let s = A[i][j];
    for (let k = 0; k < j; k++) s -= M[i][k] * M[j][k];
    if (i === j) { if (s <= 1e-12) return null; M[i][i] = Math.sqrt(s); }
    else M[i][j] = s / M[j][j];
  }
  return M;
}
function cholSolve(M, b, n) {
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= M[i][k] * y[k]; y[i] = s / M[i][i]; }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= M[k][i] * x[k]; x[i] = s / M[i][i]; }
  return x;
}

let layerSel = [], layerCov = null, layerFrac = null;
function layerCasts() {                              // re-extract only when the edge fraction or the layer moves
  const key = `${THERMO.frac}|${THERMO.kind}`;
  if (layerFrac === key) return layerSel;
  layerFrac = key;
  layerSel = casts.map(c => castLayerOf(c, THERMO.frac));
  layerCov = fitLayerCov(layerSel);
  return layerSel;
}
// the map: top, bot and skill on a gn x gn grid over the box. All three are
// smooth, so the mesh samples this bilinearly — a depth survives interpolation,
// which is exactly what a threshold does not.
const MAP = { top: null, bot: null, skTop: null, skBot: null, n: 0, ok: false, skKey: null, skT: 0 };
// The covariance between casts depends on their own separations in space and
// time, never on where the playhead is — only the vector c(x) does. So the
// factorisation, the weights and the Cholesky factor are cached and the
// playhead only re-runs the evaluation loop. Two systems now, one per surface.
let OI = null, OIkey = '';
function oiSystem(idx, X, Z, N, key, cov, sel) {
  // The uncertainty weight moves the nugget, so it changes the *geometry* and
  // not only the shading. At 1 the nugget is the cross-validated one, which
  // damps the analysis toward the basin mean by however much the array cannot
  // resolve — a lot for the top, much less for the base. At 0 the nugget goes
  // to LAY.nugMin and the same equations become an interpolant that passes
  // through every measured depth, with relief that nothing cross-validates.
  const nug = cov.nug * UNC + LAY.nugMin * (1 - UNC);
  const L = cov.L;
  // with the timeline on the covariance is space-time: two casts far apart in
  // time tell you less about each other, exactly as they do about the field.
  // At the depth this layer lives, though, L_t is pinned at the record length
  // — 30 days resolves no decorrelation down there — so this factor is close
  // to 1 for every pair, and the deep layer barely moves with the playhead.
  // That is a measured absence, not a claim that the deep ocean is steady.
  const Lt = ltAt(800);
  const rho = (dx, dz, dt) => (1 - nug) * Math.exp(-(dx * dx + dz * dz) / (2 * L * L))
    * (TIME.on ? Math.exp(-dt * dt / (2 * Lt * Lt)) : 1);
  const A = [];
  for (let i = 0; i < N; i++) {
    A.push(new Float64Array(N));
    for (let j = 0; j < N; j++)
      A[i][j] = rho(X[i] - X[j], Z[i] - Z[j], casts[idx[i]].t - casts[idx[j]].t) + (i === j ? nug : 0);
  }
  const M = cholesky(A, N);
  if (!M) return null;
  const d = idx.map(i => sel[i][key]);
  const mu = d.reduce((x, y) => x + y, 0) / N;
  const w = cholSolve(M, Float64Array.from(d.map(v => v - mu)), N);
  // The error variance needs c'A^-1 c at every grid node. A = LL', so that is
  // |L^-1 c|^2 — one forward substitution, half the work of a full solve and no
  // inverse to form. M is flattened for the inner loop.
  const Mf = new Float64Array(N * N);
  for (let i = 0; i < N; i++) for (let j = 0; j <= i; j++) Mf[i * N + j] = M[i][j];
  // The calibration rides on the uncertainty weight, because it is the same
  // statement: at 1 the error map is scaled to the error it is measured to
  // make, at 0 it reports the raw c'A^-1 c an ordinary pipeline would print —
  // which is the over-confident render the slider exists to contrast against,
  // and which is also what keeps the volume closed over the whole disc at 0.
  const cal = 1 + (cov.cal - 1) * UNC;
  return { rho, mu, w, Mf, nug, cal, L, Lt };
}
function oiSetup() {
  const key = `${THERMO.frac}|${THERMO.kind}|${TIME.on}|${TIME.scale}|${UNC}`;
  if (OI && OIkey === key) return OI;
  const sel = layerCasts();
  const idx = sel.map((r, i) => (r ? i : -1)).filter(i => i >= 0);
  const N = idx.length;
  OIkey = key;
  if (N < 4) { OI = null; return null; }
  const X = idx.map(i => xOf(casts[i].lon)), Z = idx.map(i => zOf(casts[i].lat));
  const top = oiSystem(idx, X, Z, N, 'top', layerCov.top, sel);
  const bot = oiSystem(idx, X, Z, N, 'bot', layerCov.bot, sel);
  if (!top || !bot) { OI = null; return null; }
  // The covariance between a grid node and a cast separates: c_j(node, t) =
  // Cs_j(node) * f_j(t), the spatial part times the temporal factor. Cs is
  // the same at every playhead position, so it is formed here, once per
  // system, and a tick pays N exps for the f's instead of gn * gn * N. Nodes
  // more than two cells outside the domain disc are never sampled by the mesh
  // (bilinear reads the cell a point is in, and every point is inside the
  // disc), so they are marked and skipped.
  const gn = LAY.gn, cell = 2 * RC / (gn - 1);
  const rOK = new Uint8Array(gn * gn);
  for (const sy of [top, bot]) {
    sy.Cs = new Float32Array(gn * gn * N);
    for (let b = 0; b < gn; b++) for (let a = 0; a < gn; a++) {
      const pxx = (a / (gn - 1) * 2 - 1) * RC, pz = (b / (gn - 1) * 2 - 1) * RC, o = b * gn + a;
      rOK[o] = Math.hypot(pxx, pz) <= RC + 2 * cell ? 1 : 0;
      for (let j = 0; j < N; j++)
        sy.Cs[o * N + j] = (1 - sy.nug) * Math.exp(-((pxx - X[j]) ** 2 + (pz - Z[j]) ** 2) / (2 * sy.L * sy.L));
    }
  }
  OI = { idx, N, X, Z, top, bot, sel, rOK };
  return OI;
}
function solveLayer() {
  MAP.ok = false;
  const oi = oiSetup();
  if (!oi) return;
  const { idx, N, X, Z } = oi;
  const gn = LAY.gn;
  MAP.n = gn;
  if (!MAP.top || MAP.top.length !== gn * gn) {
    MAP.top = new Float32Array(gn * gn); MAP.bot = new Float32Array(gn * gn);
    MAP.skTop = new Float32Array(gn * gn); MAP.skBot = new Float32Array(gn * gn);
  }
  // The estimate costs N per node; the error variance costs N^2/2 (the
  // forward substitution) and is what a tick spends its time on. Its only
  // dependence on the playhead is through the temporal factors, whose scale
  // is L_t(800) = the record length, so over 0.5 d of playhead they move by
  // under 1 %: the variance is re-solved every 0.5 d (or when the system
  // changes), the estimate every call.
  const doSk = MAP.skKey !== OIkey || !TIME.on || Math.abs(TIME.t - MAP.skT) >= 0.5;
  if (doSk) { MAP.skKey = OIkey; MAP.skT = TIME.t; }
  const c = new Float64Array(N), u = new Float64Array(N), f = new Float64Array(N);
  const surf = [[oi.top, MAP.top, MAP.skTop], [oi.bot, MAP.bot, MAP.skBot]];
  for (const [sy, val, sk] of surf) {
    // the temporal factor per cast at this playhead (1 with the timeline off);
    // the spatial part per node is in sy.Cs, see oiSetup
    for (let j = 0; j < N; j++) f[j] = TIME.on ? timeWeight(TIME.t - casts[idx[j]].t, sy.Lt) : 1;
    const Cs = sy.Cs, Mf = sy.Mf, w = sy.w;
    for (let o = 0; o < gn * gn; o++) {
      if (!oi.rOK[o]) { val[o] = clamp01(sy.mu / G.presMax) * G.presMax; sk[o] = 0; continue; }
      let v = sy.mu;
      const co = o * N;
      for (let j = 0; j < N; j++) {
        c[j] = Cs[co + j] * f[j];
        v += w[j] * c[j];
      }
      val[o] = clamp01(v / G.presMax) * G.presMax;
      if (!doSk) continue;
      let s2 = 0;                                    // |L^-1 c|^2
      for (let i = 0; i < N; i++) {
        let t = c[i];
        const row = i * N;
        for (let k = 0; k < i; k++) t -= Mf[row + k] * u[k];
        t /= Mf[row + i];
        u[i] = t; s2 += t * t;
      }
      sk[o] = clamp01(1 - sy.cal * (1 - Math.min(1, s2)));
    }
  }
  // a layer cannot be inverted by the map
  for (let i = 0; i < gn * gn; i++) if (MAP.bot[i] < MAP.top[i]) MAP.bot[i] = MAP.top[i];
  MAP.ok = true;
}
// ?synth=1 — a SYNTHETIC layer, for concept figures only. Every depth,
// thickness and skill value below is invented: no cast informs it, and it
// exists because the real June array maps the seasonal thermocline at ~12 % of
// its variance, which is honest and visually flat. Any figure rendered from it
// must be captioned as synthetic / illustrative. It replaces only the layer
// maps: the casts, the field and the cloud still read the measured data.
// The on-screen default is set for a console being read close up next to the
// slices; alone on a light ground at high exaggeration the layer nearly
// vanishes. aPeak rescales the whole alpha ramp, so the fade still encodes the
// same skill — this is exposure, not encoding (same knob as ?apeak).
if (SHOWCASE) {
  LAY.aPeak = 0.6;
  // no skirt between the two sheets: it closes the volume off at the rim, where
  // the analysis has least to stand on, and reads as a hard edge to a surface
  // whose whole point is that it fades out
  THERMO.walls = false;
}

let SYNTH = false;
function synthLayer() {
  const gn = LAY.gn;
  MAP.n = gn;
  MAP.top = new Float32Array(gn * gn); MAP.bot = new Float32Array(gn * gn);
  MAP.skTop = new Float32Array(gn * gn); MAP.skBot = new Float32Array(gn * gn);
  const sig = t => 1 / (1 + Math.exp(-t));
  for (let b = 0; b < gn; b++) {
    const z = b / (gn - 1) * 2 - 1;                  // -1..1 across the disc
    for (let a = 0; a < gn; a++) {
      const x = a / (gn - 1) * 2 - 1, o = b * gn + a;
      // a meandering front: a shallow, thin thermocline on one side dropping
      // to a deep, thick one across a sigmoid whose axis wanders
      const along = (x - z) * 0.7071, across = (x + z) * 0.7071;
      const meander = 0.28 * Math.sin(2.4 * along + 0.7) + 0.10 * Math.sin(5.1 * along - 1.9);
      const f = sig((across + meander) / 0.16);
      // an eddy doming the deep side back up
      const eddy = Math.exp(-((x - 0.38) ** 2 + (z + 0.30) ** 2) / (2 * 0.16 ** 2));
      const core = 34 + 96 * f - 42 * eddy + 6 * Math.sin(3.3 * x + 1.2) * Math.sin(2.9 * z);
      const th = 16 + 46 * f + 20 * eddy;
      MAP.top[o] = Math.max(4, core - th / 2);
      MAP.bot[o] = core + th / 2;
      // skill: a survey track crossing the front, three stations off it, and a
      // weak broad background — so the VSUP grading and the fade have real
      // structure: vivid along the track, ghostly between, gone at the rim
      const dTrack = Math.abs(across - 0.35 * Math.sin(1.8 * along));
      let sk = 0.96 * Math.exp(-dTrack * dTrack / (2 * 0.13 ** 2));
      for (const [px2, pz2, w] of [[-0.55, 0.5, 0.8], [0.5, 0.55, 0.85], [-0.35, -0.62, 0.75]])
        sk = Math.max(sk, w * Math.exp(-((x - px2) ** 2 + (z - pz2) ** 2) / (2 * 0.18 ** 2)));
      sk = Math.max(sk, 0.3 * Math.exp(-(x * x + z * z) / (2 * 0.55 ** 2)));
      MAP.skBot[o] = sk;
      MAP.skTop[o] = sk * 0.9;
    }
  }
  MAP.ok = true;
}
const mapAt = (arr, x, z) => {                       // bilinear sample of the OI grid
  const gn = MAP.n;
  const fx = clamp01((x / RC + 1) / 2) * (gn - 1), fy = clamp01((z / RC + 1) / 2) * (gn - 1);
  const i = Math.min(Math.floor(fx), gn - 2), j = Math.min(Math.floor(fy), gn - 2);
  const wx = fx - i, wy = fy - j, o = j * gn + i;
  return arr[o] * (1 - wx) * (1 - wy) + arr[o + 1] * wx * (1 - wy)
    + arr[o + gn] * (1 - wx) * wy + arr[o + gn + 1] * wx * wy;
};

// --- 4. mesh: one volume, cut only where we cut it ---
const FN = FX * FY;
const LF = {
  top: new Float32Array(FN), bot: new Float32Array(FN), sk: new Float32Array(FN),
  m: new Float32Array(FN), inc: new Uint8Array(FN), soft: new Uint8Array(FN),
  x: new Float32Array(FN), z: new Float32Array(FN), ok: new Uint8Array(FN),
  topF: new Float32Array(FN), botF: new Float32Array(FN), skF: new Float32Array(FN),
  skT: new Float32Array(FN), skTF: new Float32Array(FN),
};
// This layer gets its OWN fine grid, spanning the domain circle rather than the
// lon/lat box the other layers use. Nothing about the analysis stops at the box
// — the OI is defined on cast positions and MAP is already evaluated over
// [-RC, RC] — so a mesh clipped to the box was drawing a rectangle where the
// domain is a disc, and it showed as soon as the skill stopped cutting the
// shape for us. The cells are anisotropic (the domain is square in world units
// but FX != FY); the circle is still a circle.
const lx = fx => (fx / (FX - 1) * 2 - 1) * RC;
const lz = fy => (fy / (FY - 1) * 2 - 1) * RC;
const LCELL = Math.max(2 * RC / (FX - 1), 2 * RC / (FY - 1));
for (let fy = 0; fy < FY; fy++) for (let fx = 0; fx < FX; fx++) {
  const i = fy * FX + fx;
  LF.x[i] = lx(fx); LF.z[i] = lz(fy); LF.ok[i] = 1;
}
// how far inside the drawn region a node is: skill above the floor, a hard stop
// where the analysis knows nothing at all, and the edge of the domain disc. The
// radial term is scaled by one cell so the rim interpolation below lands on
// r = RC rather than on a staircase of cell corners. At uncertainty weight 0
// the first two terms are positive everywhere, so the surface closes over the
// whole disc — and the disc is the only thing still bounding it.
// THERMO.mask = false drops both skill terms and lets the analysis run to the
// edge of the disc. It is for illustrating what a surface looks like, and it
// draws the field where the covariance has extrapolated it back to the basin
// mean — any figure using it has to say so.
const layerMargin = (sk, r) => (THERMO.mask
  ? Math.min((supEff(sk) - LAY.floor) / LAY.floor, (sk - LAY.eps) / LAY.eps, (RC - r) / LCELL)
  : (RC - r) / LCELL);

function layerField() {
  LF.inc.fill(0); LF.soft.fill(0);
  for (let i = 0; i < FN; i++) {
    LF.x[i] = lx(i % FX); LF.z[i] = lz((i - i % FX) / FX);
    // the drawn region is decided by the BASE's skill — that is the surface
    // the array resolves and the one that gives the object its shape. The top
    // carries its own, much lower skill into its own shading.
    const sk = mapAt(MAP.skBot, LF.x[i], LF.z[i]);
    LF.skT[i] = mapAt(MAP.skTop, LF.x[i], LF.z[i]);
    LF.top[i] = mapAt(MAP.top, LF.x[i], LF.z[i]);
    LF.bot[i] = mapAt(MAP.bot, LF.x[i], LF.z[i]);
    LF.sk[i] = sk;
    LF.m[i] = layerMargin(sk, Math.hypot(LF.x[i], LF.z[i]));
    if (LF.m[i] >= 0) LF.inc[i] = 1;
  }
  // rim: move each outside node that touches the surface onto the contour where
  // the margin crosses zero, so the cut follows the skill contour and not the mesh
  const NB = [-1, 1, -FX, FX];
  for (let i = 0; i < FN; i++) {
    if (!LF.ok[i] || LF.inc[i]) continue;
    const fx = i % FX;
    let n = 0, X = 0, Z = 0, tp = 0, bt = 0, sk = 0, skt = 0;
    for (const d of NB) {
      const j = i + d;
      if (j < 0 || j >= FN || LF.inc[j] !== 1) continue;
      if ((d === -1 && fx === 0) || (d === 1 && fx === FX - 1)) continue;
      const s = LF.m[j] / (LF.m[j] - LF.m[i]);       // in (0,1]: 0 at j, 1 at i
      X += LF.x[j] + s * (lx(fx) - LF.x[j]);
      Z += LF.z[j] + s * (lz((i - fx) / FX) - LF.z[j]);
      tp += LF.top[j] + s * (LF.top[i] - LF.top[j]);
      bt += LF.bot[j] + s * (LF.bot[i] - LF.bot[j]);
      sk += LF.sk[j] + s * (LF.sk[i] - LF.sk[j]);
      skt += LF.skT[j] + s * (LF.skT[i] - LF.skT[j]);
      n++;
    }
    if (!n) continue;
    LF.x[i] = X / n; LF.z[i] = Z / n;
    LF.top[i] = tp / n; LF.bot[i] = bt / n; LF.sk[i] = sk / n; LF.skT[i] = skt / n;
    LF.inc[i] = 2; LF.soft[i] = 1;                   // 2 = rim, never snapped to
  }
  LF.topF.fill(NaN); LF.botF.fill(NaN); LF.skF.fill(0); LF.skTF.fill(0);
  for (let i = 0; i < FN; i++) {
    if (!LF.inc[i]) continue;
    LF.topF[i] = LF.top[i]; LF.botF[i] = LF.bot[i]; LF.skF[i] = LF.sk[i]; LF.skTF[i] = LF.skT[i];
  }
}

// --- render 'sheet': ONE opaque surface through the middle of the layer ---
//
// Geometry is the layer's own centre depth — the shape already says where it
// is, so colour is free to say something else, and what it says is how THICK
// the layer is there. Nothing is encoded in opacity: the sheet is solid, which
// is a deliberate simplification for a still. Uncertainty is still in the
// colour (the same VSUP tree: fewer thickness bins, then neutral, as the
// cross-validated skill falls) and still in the cut — the sheet stops where
// holding a cast out shows the analysis has no skill left.
const sheetMat = isoMat.clone();
sheetMat.transparent = false;
sheetMat.depthWrite = true;
function buildSheet(topF, botF, skF, skTF) {
  const mid = new Float32Array(topF.length), sk = new Float32Array(topF.length);
  for (let i = 0; i < topF.length; i++) {
    if (Number.isNaN(topF[i]) || Number.isNaN(botF[i])) { mid[i] = NaN; continue; }
    mid[i] = (topF[i] + botF[i]) / 2;
    // the sheet is one object made of two analyses, so it can be no better
    // determined than the weaker of them
    sk[i] = Math.min(skF[i], skTF[i]);
  }
  const S1 = heightSurface(mid, sk, thermoColor, () => 1, LAY.slope,
    { xs: LF.x, zs: LF.z, soft: LF.soft });
  thermoGroup.add(new THREE.Mesh(S1.geo, sheetMat));
}

// --- render 'cloud': the field inside the layer, as points ---
//
// DENSITY is the certainty channel: candidate points on a jittered lattice are
// kept with probability support^gamma, so the cloud is dense where the casts
// constrain the field, thins as they stop, and is empty where nothing was
// measured. Colour is still temperature through the VSUP tree (greying to
// slate with the same support), so the two encodings agree. The jitter and the
// keep-test use a seeded PRNG, so a rebuild reproduces the same cloud.
const cloudMat = new THREE.ShaderMaterial({
  vertexShader: `
    attribute vec3 col; attribute float alpha; uniform float size;
    varying vec3 vC; varying float vA;
    void main(){ vC = col; vA = alpha;
      vec4 mv = modelViewMatrix * vec4(position, 1.0);
      gl_PointSize = size * (3000.0 / -mv.z);        // scene units are km; the box is ~1300 across
      gl_Position = projectionMatrix * mv; }`,
  fragmentShader: `
    varying vec3 vC; varying float vA;
    void main(){ vec2 d = gl_PointCoord - vec2(0.5);
      if (dot(d, d) > 0.25) discard;                 // round points, not squares
      gl_FragColor = vec4(vC, vA); }`,
  uniforms: { size: { value: 7 } },
  transparent: true, depthWrite: false,
});
const CLOUD = { res: 1, gamma: 1.5, dens: 1 };       // candidates per grid cell per axis; keep-prob = dens * conf^gamma
function buildCloud() {
  const t0 = performance.now();
  const pos = [], col = [], alp = [];
  // mulberry32, fixed seed: the cloud is stochastic but a capture must reproduce
  let rs = 0x9e3779b9;
  const rnd = () => {
    rs |= 0; rs = rs + 0x6D2B79F5 | 0;
    let t = Math.imul(rs ^ rs >>> 15, 1 | rs);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
  // CLOUD.res > 1 places candidates between the field's grid nodes, through the
  // same trilinear sampler the readout and the 2D insets use. It is
  // interpolation of an already-smoothed field (25 km, 30 dbar kernel), not new
  // information, so a caption has to say the point spacing is finer than the
  // grid — the density gradient is the information, not the point positions.
  const r = CLOUD.res, sx = Math.round((NX - 1) * r), sy = Math.round((NY - 1) * r),
    sk = Math.round((NZ - 1) * r);
  const dlon = (G.lon1 - G.lon0) / sx, dlat = (G.lat1 - G.lat0) / sy, dp = G.presMax / sk;
  for (let b = 0; b <= sy; b++) for (let a = 0; a <= sx; a++) {
    for (let c = 0; c <= sk; c++) {
      // jitter inside the cell: with density carrying the encoding, a regular
      // lattice with holes would read as moiré, not as thinning support
      const lon = G.lon0 + a * dlon + (rnd() - 0.5) * dlon,
        lat = G.lat0 + b * dlat + (rnd() - 0.5) * dlat,
        p = c * dp + (rnd() - 0.5) * dp;
      const X = xOf(lon), Z = zOf(lat);
      if (Math.hypot(X, Z) > RC) continue;
      const tp = mapAt(MAP.top, X, Z), bt = mapAt(MAP.bot, X, Z);
      if (!(bt > tp) || p < tp || p > bt) continue;
      const f = fieldAt(lon, lat, p);
      if (!(f.conf > 0.01)) continue;                // nothing measured here
      if (rnd() >= Math.min(1, CLOUD.dens * Math.pow(f.conf, CLOUD.gamma))) continue;
      pos.push(X, yOf(p), Z);
      // same tree and same palette as the surfaces and as the proposal's VSUP
      // figure: value on the magma ramp, greying to slate as support falls
      col.push(...vsupColor(f.temp, f.conf, 'temp'));
      alp.push(0.85);                                // constant: density carries the certainty now
    }
  }
  console.log(`cloud: ${pos.length / 3} points, res ${r}, ${(performance.now() - t0).toFixed(0)} ms`);
  if (!pos.length) return;
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  geo.setAttribute('col', new THREE.Float32BufferAttribute(col, 3));
  geo.setAttribute('alpha', new THREE.Float32BufferAttribute(alp, 1));
  thermoGroup.add(new THREE.Points(geo, cloudMat));
}

function buildThermo() {
  for (const ch of [...thermoGroup.children]) { ch.geometry.dispose(); thermoGroup.remove(ch); }
  for (const L of thermoLabels) { L.el.remove(); labels.splice(labels.indexOf(L), 1); }
  thermoLabels.length = 0;
  if (!THERMO.on) return;
  if (SYNTH) synthLayer(); else solveLayer();
  if (!MAP.ok) return;
  layerField();
  const topF = LF.topF, botF = LF.botF, skF = LF.skF, skTF = LF.skTF;
  // thickness per node and the range the colour ramp spans, for either sheet
  // render. Computed here so the slab, the single sheet and the walls all read
  // the same array and the same scale.
  if (LTH.length !== topF.length) LTH = new Float32Array(topF.length);
  {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < topF.length; i++) {
      LTH[i] = botF[i] - topF[i];
      if (!LF.inc[i] || Number.isNaN(LTH[i])) continue;
      if (THERMO.colorBy === 'thick') {
        if (LTH[i] < lo) lo = LTH[i];
        if (LTH[i] > hi) hi = LTH[i];
      } else {
        for (const d of [topF[i], botF[i]]) {
          if (d < lo) lo = d;
          if (d > hi) hi = d;
        }
      }
    }
    if (lo < hi) { THERMO.lo = lo; THERMO.hi = hi; }
  }
  if (THERMO.render === 'sheet') { buildSheet(topF, botF, skF, skTF); return; }
  if (THERMO.render === 'cloud') { buildCloud(); return; }
  // Support here is the CALIBRATED error variance (see fitLayerCov), so it is
  // on the same footing as the kernel proxy the other layers use and gets the
  // same treatment: full indigo and full alpha only where the analysis has been
  // shown to be right, greying up the VSUP tree and fading to nothing as it has
  // not. Earlier versions kept a floor of opacity here because the raw skill
  // was near-saturated and the surface would otherwise have been drawn solid
  // everywhere OR invisible everywhere, neither of which was true; with the
  // calibration in, the honest fade is also the readable one. Over the disc the
  // base now spans all five VSUP layers and 26 % of it is below the floor and
  // not drawn at all — the outer ring, where there are no casts.
  const aFloor = LAY.alphaFloor, aSpan = 1 - aFloor;
  // Peak opacity is back at 0.34 now that the floor is gone. It was dropped to
  // 0.22 when this layer got five times thicker than the seasonal one it
  // replaced and the floor was still holding the whole sheet visible; with the
  // floor at 0 the thin end takes care of itself, and the well-supported end
  // has to read as solid or the fade encodes nothing.
  const alpha = s => (aFloor + aSpan * Math.pow(supportAlpha(s), 1.6)) * LAY.aPeak;
  const opt = { xs: LF.x, zs: LF.z, soft: LF.soft };
  // each sheet shaded by its OWN error variance: the base is resolved and
  // draws near full indigo, the top is not and stays close to neutral. Same
  // object, two very different levels of evidence, and that shows.
  const T1 = heightSurface(topF, skTF, thermoColor, alpha, LAY.slope, opt);
  const B1 = heightSurface(botF, skF, thermoColor, alpha, LAY.slope, opt);
  thermoGroup.add(new THREE.Mesh(T1.geo, slabMat), new THREE.Mesh(B1.geo, slabMat));
  // fishnet over each sheet, the same one the iso-surface carries: without it a
  // near-flat translucent sheet has no cue to read its shape from except its
  // silhouette, and the relief here (the base bows ~300 dbar across the box) is
  // exactly the thing worth seeing. Same vertices, every LAY.net'th row and
  // column, same slope guard, and the line takes its own vertex's alpha — so
  // the net thins out with the sheet instead of drawing a confident grid over
  // ground the analysis does not stand on.
  for (const S1 of [T1, B1]) thermoGroup.add(fishnet(S1, LAY.net, LAY.slope,
    { gain: 1 / LAY.aPeak, ink: 0.05, mat: netMat }));
  // wall on every node edge where the meshed region ends. Single-sided and
  // wound outward, so the far side of the volume is culled instead of laid over
  // the near side. Every wall here is a cut we made — the analysis does not
  // stop on its own — so they are all honest edges of the drawn region.
  const okNode = i => LF.inc[i] !== 0;
  const okCell = (cx, cy) => cx >= 0 && cy >= 0 && cx < FX - 1 && cy < FY - 1 &&
    okNode(cy * FX + cx) && okNode(cy * FX + cx + 1) && okNode((cy + 1) * FX + cx) && okNode((cy + 1) * FX + cx + 1);
  const sp = [], sc = [], sa = [], si = [];
  const cellMid = (cx, cy) => [(LF.x[cy * FX + cx] + LF.x[cy * FX + cx + 1]) / 2,
  (LF.z[cy * FX + cx] + LF.z[(cy + 1) * FX + cx]) / 2];
  const wall = (i0, i1, inX, inY) => {
    if (!okNode(i0) || !okNode(i1)) return;
    if (Math.abs(botF[i0] - topF[i0]) < 1 && Math.abs(botF[i1] - topF[i1]) < 1) return;
    const n = sp.length / 3;
    for (const [i, d] of [[i0, topF[i0]], [i1, topF[i1]], [i1, botF[i1]], [i0, botF[i0]]]) {
      sp.push(LF.x[i], yOf(d), LF.z[i]);
      sc.push(...thermoColor(skF[i], d, i));         // the wall shades from top depth to base depth
      sa.push((aFloor + aSpan * Math.pow(supportAlpha(skF[i]), 1.6)) * 0.46);
    }
    const ex = LF.x[i1] - LF.x[i0], ez = LF.z[i1] - LF.z[i0];
    const nx = ez, nz = -ex;                         // face normal of the quad as wound below
    const [mx, mz] = cellMid(inX, inY);
    const out = nx * (LF.x[i0] - mx) + nz * (LF.z[i0] - mz) > 0;
    if (out) si.push(n, n + 1, n + 2, n, n + 2, n + 3);
    else si.push(n, n + 2, n + 1, n, n + 3, n + 2);
  };
  for (let fy = 0; fy < FY - 1; fy++) for (let fx = 0; fx < FX; fx++) {
    const l = okCell(fx - 1, fy), r = okCell(fx, fy);
    if (l !== r) wall(fy * FX + fx, (fy + 1) * FX + fx, l ? fx - 1 : fx, fy);
  }
  for (let fx = 0; fx < FX - 1; fx++) for (let fy = 0; fy < FY; fy++) {
    const u = okCell(fx, fy - 1), d = okCell(fx, fy);
    if (u !== d) wall(fy * FX + fx, fy * FX + fx + 1, fx, u ? fy - 1 : fy);
  }
  // ?walls=off leaves the two sheets unjoined — a concept still may want the
  // surfaces alone; interactively the walls mark every cut as a cut, so they
  // stay the default
  if (si.length && THERMO.walls) {
    const sgeo = new THREE.BufferGeometry();
    sgeo.setAttribute('position', new THREE.Float32BufferAttribute(sp, 3));
    sgeo.setAttribute('col', new THREE.Float32BufferAttribute(sc, 3));
    sgeo.setAttribute('alpha', new THREE.Float32BufferAttribute(sa, 1));
    sgeo.setIndex(si);
    sgeo.computeVertexNormals();
    thermoGroup.add(new THREE.Mesh(sgeo, wallMat));
  }
  // collar on each tube: where that cast actually measured the layer. These are
  // data; the sheets between them are the analysis.
  const sel = layerCasts();
  const collar = new THREE.CylinderGeometry(TUBE_R * 1.9, TUBE_R * 1.9, 0.004, 12, 1, true);
  const cmat = new THREE.MeshBasicMaterial({ color: new THREE.Color(...THERMO_RGB) });
  const places = [];
  casts.forEach((c, ci) => {
    const r = sel[ci];
    if (!r || !$('tg-casts').checked || !castVisible(c)) return;
    places.push(xOf(c.lon), yOf(r.top), zOf(c.lat), xOf(c.lon), yOf(r.bot), zOf(c.lat));
  });
  if (places.length) thermoGroup.add(new THREE.Mesh(mergedCylinders(collar, places), cmat));
  collar.dispose();
  for (const m of thermoGroup.children) m.renderOrder = 2;   // after the iso-surface, see buildIsoSheets
  let best = null;
  for (let i = 0; i < FN; i++) {
    if (LF.inc[i] !== 1) continue;
    if (!best || skF[i] > best.sk) best = { sk: skF[i], th: botF[i] - topF[i], i, mid: (topF[i] + botF[i]) / 2 };
  }
  // no floating label on this layer in the showcase: the caption already
  // names it, and a tooltip pinned to the deepest node is one more thing
  // over a scene that is meant to be looked at
  if (best && !SHOWCASE) {
    const L = addLabel(`${SHOWCASE ? 'main thermocline' : 'permanent pycnocline'} ${best.th.toFixed(0)} dbar thick`,
      LF.x[best.i], yOf(best.mid), LF.z[best.i], 'iso');
    L.el.style.color = css(THERMO_RGB);
    thermoLabels.push(L);
  }
}
// layer statistics for the caption: measured at the casts, plus how much of the
// box the analysis actually resolves
function thermoStats() {
  const sel = layerCasts();
  if (!MAP.ok) solveLayer();
  const q = (a, f) => (a.length ? a[Math.floor(f * (a.length - 1))] : NaN);
  const have = sel.filter(Boolean);
  const th = have.map(r => r.bot - r.top).sort((a, b) => a - b);
  const tops = have.map(r => r.top).sort((a, b) => a - b);
  const bots = have.map(r => r.bot).sort((a, b) => a - b);
  let drawn = 0, n = 0;
  if (MAP.ok) for (let i = 0; i < MAP.n * MAP.n; i++) {
    const a = i % MAP.n, b = (i - a) / MAP.n;
    const X = (a / (MAP.n - 1) * 2 - 1) * RC, Z = (b / (MAP.n - 1) * 2 - 1) * RC;
    if (Math.hypot(X, Z) > RC) continue;
    n++; if (MAP.skBot[i] >= LAY.floor) drawn++;
  }
  const cores = have.map(r => r.core).sort((a, b) => a - b);
  const pks = have.map(r => r.pk).sort((a, b) => a - b);
  const cv = layerCov || { top: {}, bot: {} };
  return {
    nCast: have.length, nTrunc: have.filter(r => r.trunc).length,
    med: q(th, 0.5), thin: q(th, 0.05), thick: q(th, 0.95),
    medTop: q(tops, 0.5), medBot: q(bots, 0.5), topP90: q(tops, 0.9),
    medCore: q(cores, 0.5), medPk: q(pks, 0.5), botLo: q(bots, 0.05), botHi: q(bots, 0.95),
    topLo: q(tops, 0.05), topHi: q(tops, 0.95),
    L: cv.bot.L, nug: cv.bot.nug, sd: cv.bot.sd, loo: cv.bot.loo, r2: cv.bot.r2,
    tL: cv.top.L, tNug: cv.top.nug, tSd: cv.top.sd, tLoo: cv.top.loo, tR2: cv.top.r2,
    cal: cv.bot.cal, predSd: cv.bot.predSd, tCal: cv.top.cal, tPredSd: cv.top.predSd,
    nugEff: OI ? OI.bot.nug : NaN, cover: n ? drawn / n : 0
  };
}

// ---------- interpolated slice + section, colour = variable, holes = low support ----------

function trilin(arr, gx, gy, kf) {
  const x0 = Math.min(Math.max(Math.floor(gx), 0), NX - 2), wx = clamp01(gx - x0);
  const y0 = Math.min(Math.max(Math.floor(gy), 0), NY - 2), wy = clamp01(gy - y0);
  const k0 = Math.min(Math.max(Math.floor(kf), 0), NZ - 2), wk = clamp01(kf - k0);
  const bl = (k) => {
    const base = k * NY * NX + y0 * NX + x0;
    return arr[base] * (1 - wx) * (1 - wy) + arr[base + 1] * wx * (1 - wy)
      + arr[base + NX] * (1 - wx) * wy + arr[base + NX + 1] * wx * wy;
  };
  return bl(k0) * (1 - wk) + bl(k0 + 1) * wk;
}
// one variable and its support at (lon, lat, pres), into a shared object.
// The slice and the section sample 40k and 24k points per rebuild; fieldAt's
// four trilinears and a fresh object per sample were most of what they cost.
const SMP = { v: 0, conf: 0 };
function fieldSample(lon, lat, pres) {
  const gx = (lon - G.lon0) / (G.lon1 - G.lon0) * (NX - 1);
  const gy = (lat - G.lat0) / (G.lat1 - G.lat0) * (NY - 1);
  const kf = clamp01(pres / G.presMax) * (NZ - 1);
  const inBox = gx >= 0 && gx <= NX - 1 && gy >= 0 && gy <= NY - 1;
  const cx = Math.min(Math.max(gx, 0), NX - 1), cy = Math.min(Math.max(gy, 0), NY - 1);
  const x0 = Math.min(Math.max(Math.floor(cx), 0), NX - 2), wx = clamp01(cx - x0);
  const y0 = Math.min(Math.max(Math.floor(cy), 0), NY - 2), wy = clamp01(cy - y0);
  const k0 = Math.min(Math.max(Math.floor(kf), 0), NZ - 2), wk = clamp01(kf - k0);
  const b0 = k0 * NY * NX + y0 * NX + x0, b1 = b0 + NY * NX;
  const w00 = (1 - wx) * (1 - wy), w10 = wx * (1 - wy), w01 = (1 - wx) * wy, w11 = wx * wy;
  const bl = (arr, o) => arr[o] * w00 + arr[o + 1] * w10 + arr[o + NX] * w01 + arr[o + NX + 1] * w11;
  const conf = inBox ? bl(F.conf, b0) * (1 - wk) + bl(F.conf, b1) * wk : 0;
  SMP.conf = conf;
  SMP.v = curVar === 'conf' ? conf : bl(F[curVar], b0) * (1 - wk) + bl(F[curVar], b1) * wk;
  return SMP;
}
// sample the field at (lon, lat, pres); null outside the gridded box
function fieldAt(lon, lat, pres) {
  const gx = (lon - G.lon0) / (G.lon1 - G.lon0) * (NX - 1);
  const gy = (lat - G.lat0) / (G.lat1 - G.lat0) * (NY - 1);
  const kf = clamp01(pres / G.presMax) * (NZ - 1);
  const inBox = gx >= 0 && gx <= NX - 1 && gy >= 0 && gy <= NY - 1;
  const cx = Math.min(Math.max(gx, 0), NX - 1), cy = Math.min(Math.max(gy, 0), NY - 1);
  return {
    temp: trilin(F.temp, cx, cy, kf), psal: trilin(F.psal, cx, cy, kf),
    sigma0: trilin(F.sigma0, cx, cy, kf), conf: inBox ? trilin(F.conf, cx, cy, kf) : 0, inBox
  };
}

// marching squares on a W x H scalar array (NaN = undefined) -> canvas lines
// marching-squares case table: edge pairs per corner code (edges 0 top, 1
// right, 2 bottom, 3 left), flat so the inner loop allocates nothing
const MS_SEGS = [null, [3, 0], [0, 1], [3, 1], [1, 2], [3, 0, 1, 2], [0, 2], [3, 2],
  [2, 3], [0, 2], [0, 1, 2, 3], [1, 2], [1, 3], [0, 1], [3, 0], null];
// One pass over the cells. The levels are sorted, so a cell only has to look
// at the ones between its lowest and highest corner (found by bisection); a
// pass per level over every cell was 20+ passes over the slice for a handful
// of crossings each.
function drawContours(ctx, vals, W, H, sx, sy, levels, style) {
  ctx.strokeStyle = style; ctx.lineWidth = 1;
  ctx.beginPath();
  const nL = levels.length;
  if (!nL) { ctx.stroke(); return; }
  let a, b, c, d, i, j, lv;
  const ex = e => (e === 0 ? (i + (lv - a) / (b - a)) * sx : e === 1 ? (i + 1) * sx : e === 2 ? (i + (lv - d) / (c - d)) * sx : i * sx);
  const ey = e => (e === 0 ? j * sy : e === 1 ? (j + (lv - b) / (c - b)) * sy : e === 2 ? (j + 1) * sy : (j + (lv - a) / (d - a)) * sy);
  const firstAbove = v => { let lo = 0, hi = nL; while (lo < hi) { const m = (lo + hi) >> 1; if (levels[m] <= v) lo = m + 1; else hi = m; } return lo; };
  for (j = 0; j < H - 1; j++) for (i = 0; i < W - 1; i++) {
    a = vals[j * W + i]; b = vals[j * W + i + 1]; c = vals[(j + 1) * W + i + 1]; d = vals[(j + 1) * W + i];
    if (Number.isNaN(a + b + c + d)) continue;
    const mn = Math.min(a, b, c, d), mx = Math.max(a, b, c, d);
    // a level crosses the cell iff mn < lv <= mx (the corner test is >=)
    for (let l = firstAbove(mn); l < nL && levels[l] <= mx; l++) {
      lv = levels[l];
      const code = (a >= lv) | ((b >= lv) << 1) | ((c >= lv) << 2) | ((d >= lv) << 3);
      const segs = MS_SEGS[code];
      if (!segs) continue;
      for (let k = 0; k < segs.length; k += 2) { ctx.moveTo(ex(segs[k]), ey(segs[k])); ctx.lineTo(ex(segs[k + 1]), ey(segs[k + 1])); }
    }
  }
  ctx.stroke();
}
function contourLevels() {
  const d = V(), out = [];
  for (let v = Math.ceil(d.lo / d.contour) * d.contour; v <= d.hi; v += d.contour) out.push(v);
  return out;
}
// paint one field image: small (W x H) samples -> big canvas with holes and contours
// the support ramp, quantised to 256 steps and kept: the continuous varColor
// allocates per call, and the slice asks for it per pixel
const confTab = [];
const confColor = c => {
  const q = Math.round(clamp01((c - VARS.conf.lo) / (VARS.conf.hi - VARS.conf.lo)) * 255);
  return confTab[q] ??= varColor(VARS.conf.lo + (VARS.conf.hi - VARS.conf.lo) * q / 255, 'conf');
};
// per-target scratch: the two small canvases, their ImageData and the sample
// arrays are made once per (canvas, W, H) and refilled, not rebuilt per tick
const PF = new Map();
function paintField(big, W, H, sample) {
  let S = PF.get(big);
  if (!S || S.W !== W || S.H !== H) {
    const small = document.createElement('canvas'), mask = document.createElement('canvas');
    small.width = mask.width = W; small.height = mask.height = H;
    const sctx = small.getContext('2d'), mctx = mask.getContext('2d');
    S = {
      W, H, small, sctx, img: sctx.createImageData(W, H), mask, mctx, mimg: mctx.createImageData(W, H),
      vals: new Float32Array(W * H), sup: new Float32Array(W * H), insideAny: new Uint8Array(W * H)
    };
    PF.set(big, S);
  }
  const { vals, sup, insideAny, img, mimg, sctx } = S;
  vals.fill(NaN); sup.fill(0); insideAny.fill(0);
  const data = img.data, mdata = mimg.data, isConf = curVar === 'conf';
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const s = sample(i, j);
    const o = (j * W + i) * 4;
    // fade: the whole layer (fill, contours, marks) is multiplied by
    // supportAlpha through the mask below; unsampled pixels are cleared
    if (!s || s.conf < CUT) { data[o + 3] = 0; mdata[o + 3] = 0; continue; }
    const p = j * W + i;
    vals[p] = s.conf > 0.01 ? s.v : NaN; sup[p] = s.conf; insideAny[p] = 1;
    const rgb = isConf ? confColor(s.conf) : vsupColor(s.v, s.conf);
    data[o] = rgb[0] * 255; data[o + 1] = rgb[1] * 255; data[o + 2] = rgb[2] * 255;
    data[o + 3] = 255;
    mdata[o + 3] = supportAlpha(s.conf) * 255;
  }
  sctx.putImageData(img, 0, 0);
  const ctx = big.getContext('2d');
  ctx.clearRect(0, 0, big.width, big.height);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(S.small, 0, 0, big.width, big.height);
  const sx = big.width / (W - 1), sy = big.height / (H - 1);
  if (contoursOn)
    drawContours(ctx, vals, W, H, sx, sy, contourLevels(), 'rgba(0,0,0,0.55)');
  if (stippleOn) {
    // screen-door: hole size grows as support falls from SOLID to CUT
    ctx.globalCompositeOperation = 'destination-out';
    ctx.fillStyle = '#000';
    const cell = 7;
    for (let y = 0; y < big.height; y += cell) for (let x = 0; x < big.width; x += cell) {
      const i = Math.round(x / big.width * (W - 1)), j = Math.round(y / big.height * (H - 1));
      const s = sup[j * W + i];
      if (Number.isNaN(vals[j * W + i]) || s >= SOLID) continue;
      const hole = cell * Math.pow(clamp01((SOLID - s) / (SOLID - CUT)), 0.7) * 0.9;
      ctx.fillRect(x + (cell - hole) / 2, y + (cell - hole) / 2, hole, hole);
    }
    ctx.globalCompositeOperation = 'source-over';
  }
  S.mctx.putImageData(mimg, 0, 0);
  ctx.globalCompositeOperation = 'destination-in';
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(S.mask, 0, 0, big.width, big.height);
  ctx.globalCompositeOperation = 'source-over';
  return ctx;
}

const slice = (() => {
  const canvas = document.createElement('canvas');
  canvas.width = 1024; canvas.height = 1024;
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.NoColorSpace;
  tex.anisotropy = renderer.capabilities.getMaxAnisotropy();
  const mesh = new THREE.Mesh(new THREE.CircleGeometry(RC, 96),
    new THREE.MeshBasicMaterial({ map: tex, side: THREE.DoubleSide, transparent: true, alphaTest: 0.02, depthWrite: true }));
  mesh.rotation.x = -Math.PI / 2;
  world.add(mesh);
  const label = addLabel('', -RC * 0.72, 0, RC * 0.72, 'axis');
  return { canvas, tex, mesh, label, pres: 300 };
})();
function updateSlice(pres) {
  slice.pres = pres;
  slice.dirty = false;
  const N = 200;
  const ctx = paintField(slice.canvas, N, N, (i, j) => {
    // CircleGeometry UV: +v = +y before rotation = north (-z)
    const x = (i / (N - 1) * 2 - 1) * RC, z = -((1 - j / (N - 1)) * 2 - 1) * RC;
    if (Math.hypot(x, z) > RC) return null;
    return fieldSample(lonOfX(x), latOfZ(z), pres);
  });
  // cast positions that reach this depth: black-outlined dots = the measurements
  const S = slice.canvas.width;
  for (const c of casts) {
    if (c.pres[c.pres.length - 1] < pres) continue;
    const u = (xOf(c.lon) / RC + 1) / 2, v = 1 - (-zOf(c.lat) / RC + 1) / 2;
    ctx.beginPath(); ctx.arc(u * S, v * S, 4, 0, Math.PI * 2);
    ctx.fillStyle = '#fff'; ctx.fill(); ctx.lineWidth = 1.5; ctx.strokeStyle = '#000'; ctx.stroke();
  }
  slice.tex.needsUpdate = true;
  slice.mesh.position.y = yOf(pres);
  slice.label.pos.y = yOf(pres);
  slice.label.el.textContent = `slice ${pres.toFixed(0)} dbar (interpolated)`;
}

const sect = (() => {
  const canvas = document.createElement('canvas');
  canvas.width = 1024; canvas.height = 512;
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.NoColorSpace;
  tex.anisotropy = renderer.capabilities.getMaxAnisotropy();
  const mesh = new THREE.Mesh(new THREE.PlaneGeometry(2 * RC, DEEP),
    new THREE.MeshBasicMaterial({ map: tex, side: THREE.DoubleSide, transparent: true, alphaTest: 0.02, depthWrite: true }));
  world.add(mesh);
  const label = addLabel('', 0, 0.02, 0, 'axis');
  return { canvas, tex, mesh, label, axis: 'lon', val: lonMid };
})();
function updateSection(axis, val) {
  sect.axis = axis; sect.val = val;
  sect.dirty = false;
  const W = 240, H = 100;
  const along = axis === 'lon';                    // N-S plane at fixed longitude
  const c0 = along ? xOf(val) : zOf(val);          // fixed coordinate
  // plane u axis: -RC..RC along z (N-S; u=0 north) or along x (E-W; u=0 west)
  const ctx = paintField(sect.canvas, W, H, (i, j) => {
    const u = (i / (W - 1) * 2 - 1) * RC;
    const x = along ? c0 : u, z = along ? -u : c0;
    if (Math.hypot(x, z) > RC) return null;
    return fieldSample(lonOfX(x), latOfZ(z), j / (H - 1) * G.presMax);
  });
  // casts within 40 km of the plane: dashed verticals = where the section is
  // actually constrained by measurements
  const CW = sect.canvas.width, CH = sect.canvas.height;
  ctx.setLineDash([6, 5]); ctx.lineWidth = 2;
  for (const c of casts) {
    const cx = xOf(c.lon), cz = zOf(c.lat);
    const dist = along ? Math.abs(cx - c0) : Math.abs(cz - c0);
    if (dist > 40) continue;
    const u = along ? -cz : cx;
    const px = (u / RC + 1) / 2 * CW;
    const p1 = c.pres[c.pres.length - 1] / G.presMax * CH;
    ctx.strokeStyle = `rgba(20,24,30,${0.95 - dist / 40 * 0.6})`;
    ctx.beginPath(); ctx.moveTo(px, 0); ctx.lineTo(px, p1); ctx.stroke();
  }
  ctx.setLineDash([]);
  sect.tex.needsUpdate = true;
  // orient the plane. PlaneGeometry: +u = +x, +v = +y, faces +z.
  sect.mesh.position.set(along ? c0 : 0, -DEEP / 2, along ? 0 : c0);
  sect.mesh.rotation.y = along ? Math.PI / 2 : 0;   // rotated: local +x -> world -z (north at u=0)
  sect.label.pos.set(along ? c0 : 0, 0.03, along ? 0 : c0);
  sect.label.el.textContent = along
    ? `section ${Math.abs(val).toFixed(1)}°W, N–S (interpolated)`
    : `section ${val.toFixed(1)}°N, E–W (interpolated)`;
}

// ---------- 2D linked views: T(p), S(p), T-S ----------

const DPR2 = Math.min(devicePixelRatio || 1, 2);
const panels = {};
for (const id of ['p-temp', 'p-psal', 'p-ts']) panels[id] = document.getElementById(id);
// the plot box is set in CSS (--plot-w / --plot-h, per breakpoint) and the
// bitmap follows the box; frameOf/axes read the bitmap, so nothing below needs
// to know the size. A hidden panel (figure mode) measures 0 and keeps its bitmap.
function sizeInsets() {
  for (const cv of Object.values(panels)) {
    const r = cv.getBoundingClientRect();
    if (!r.width || !r.height) continue;
    const w = Math.round(r.width * DPR2), h = Math.round(r.height * DPR2);
    if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
  }
}
sizeInsets();
const PM = { l: 34, r: 8, t: 8, b: 20 };           // css px margins
const PIN_COLORS = ['#d97706', '#db2777', '#059669'];
const pinLabels = [];

const tsSig = decodeU16(OCEAN.tsSigma.b64, OCEAN.tsSigma.lo, OCEAN.tsSigma.hi);
const tsN = OCEAN.tsSigma.n;
const sigAt = (u, v) => {
  const fx = u * (tsN - 1), fy = v * (tsN - 1);
  const x0 = Math.floor(fx), y0 = Math.floor(fy);
  const x1 = Math.min(x0 + 1, tsN - 1), y1 = Math.min(y0 + 1, tsN - 1);
  const gx2 = fx - x0, gy2 = fy - y0;
  return tsSig[y0 * tsN + x0] * (1 - gx2) * (1 - gy2) + tsSig[y0 * tsN + x1] * gx2 * (1 - gy2)
    + tsSig[y1 * tsN + x0] * (1 - gx2) * gy2 + tsSig[y1 * tsN + x1] * gx2 * gy2;
};

function frameOf(cv) {
  const l = PM.l * DPR2, r = PM.r * DPR2, t = PM.t * DPR2, b = PM.b * DPR2;
  return { l, t, w: cv.width - l - r, h: cv.height - t - b };
}
function axes(ctx, cv, xlo, xhi, ylo, yhi, xlab, ylab, xdec = 0, ydec = 0, yflip = false) {
  const f = frameOf(cv);
  ctx.fillStyle = '#f2f4f6'; ctx.fillRect(0, 0, cv.width, cv.height);
  ctx.fillStyle = '#ffffff'; ctx.fillRect(f.l, f.t, f.w, f.h);
  ctx.strokeStyle = '#d5dae0'; ctx.lineWidth = 1;
  ctx.fillStyle = '#2a3038'; ctx.font = `${9 * DPR2}px Cascadia Code, monospace`;
  const X = v => f.l + f.w * (v - xlo) / (xhi - xlo);
  const Y = v => yflip ? f.t + f.h * (v - ylo) / (yhi - ylo) : f.t + f.h * (1 - (v - ylo) / (yhi - ylo));
  for (const v of niceTicks(xlo, xhi, 4).ticks) {
    ctx.beginPath(); ctx.moveTo(X(v), f.t); ctx.lineTo(X(v), f.t + f.h); ctx.stroke();
    ctx.textAlign = 'center'; ctx.fillText(v.toFixed(xdec), X(v), cv.height - 6 * DPR2);
  }
  for (const v of niceTicks(ylo, yhi, 4).ticks) {
    ctx.beginPath(); ctx.moveTo(f.l, Y(v)); ctx.lineTo(f.l + f.w, Y(v)); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(v.toFixed(ydec), f.l - 3 * DPR2, Y(v) + 3 * DPR2);
  }
  ctx.textAlign = 'right'; ctx.fillStyle = '#5b6470';
  ctx.fillText(xlab, f.l + f.w, f.t + 9 * DPR2);
  ctx.textAlign = 'left'; ctx.fillText(ylab, f.l + 3 * DPR2, f.t + 9 * DPR2);
  return { X, Y };
}
function polyline(ctx, pts, style, width) {
  ctx.strokeStyle = style; ctx.lineWidth = width; ctx.beginPath();
  pts.forEach(([x, y], i) => i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y));
  ctx.stroke();
}
function drawProfilePanel(cv, key, xlab, dec) {
  const ctx = cv.getContext('2d');
  const d = VARS[key];
  const { X, Y } = axes(ctx, cv, d.lo, d.hi, 0, R.pres[1], xlab, 'dbar', dec, 0, true);
  const line = (c, style, w) => polyline(ctx, Array.from(c.pres, (p, i) => [X(c[key][i]), Y(p)]), style, w);
  for (const c of casts) if (castVisible(c)) line(c, 'rgba(60,70,80,0.22)', 1);
  pins.forEach((ci, k) => line(casts[ci], PIN_COLORS[k], 2 * DPR2));
  if (hoveredCast != null) line(casts[hoveredCast], '#000000', 1.6 * DPR2);
}
const contourCache = new Map();
function tsContour(level, X, Y) {
  if (contourCache.has(level)) return contourCache.get(level);
  const path = [];
  for (let k = 0; k <= 60; k++) {
    const u = k / 60;
    if ((sigAt(u, 0) - level) * (sigAt(u, 1) - level) > 0) continue;
    let lo = 0, hi = 1;
    for (let it = 0; it < 18; it++) {
      const mid = (lo + hi) / 2;
      if ((sigAt(u, lo) - level) * (sigAt(u, mid) - level) <= 0) hi = mid; else lo = mid;
    }
    const s = R.psal[0] + u * (R.psal[1] - R.psal[0]);
    const t = R.temp[0] + (lo + hi) / 2 * (R.temp[1] - R.temp[0]);
    path.push([X(s), Y(t), s, t]);
  }
  contourCache.set(level, path);
  return path;
}
function drawTSPanel(cv) {
  const ctx = cv.getContext('2d');
  const { X, Y } = axes(ctx, cv, R.psal[0], R.psal[1], R.temp[0], R.temp[1], 'S', '°C', 1, 0);
  ctx.font = `${8.5 * DPR2}px Cascadia Code, monospace`;
  for (let lv = 24; lv <= 27.75; lv += 0.5) {
    const p = tsContour(+lv.toFixed(2), X, Y);
    polyline(ctx, p, '#b7bec7', 1);
    if (p.length && Math.abs(lv % 1) < 0.01) {
      const q = p[Math.floor(p.length * 0.15)];
      ctx.fillStyle = '#7a828c'; ctx.textAlign = 'left'; ctx.fillText(lv.toFixed(0), q[0] + 2, q[1] - 2);
    }
  }
  const scatter = (c, style, size) => {
    ctx.fillStyle = style;
    for (let i = 0; i < c.pres.length; i += 2) ctx.fillRect(X(c.psal[i]), Y(c.temp[i]), size, size);
  };
  for (const c of casts) if (castVisible(c)) scatter(c, 'rgba(60,70,80,0.35)', 1.2 * DPR2);
  pins.forEach((ci, k) => polyline(ctx, Array.from(casts[ci].pres, (p, i) => [X(casts[ci].psal[i]), Y(casts[ci].temp[i])]), PIN_COLORS[k], 2 * DPR2));
  if (hoveredCast != null) {
    const c = casts[hoveredCast];
    polyline(ctx, Array.from(c.pres, (p, i) => [X(c.psal[i]), Y(c.temp[i])]), '#000000', 1.6 * DPR2);
  }
}
function drawInsets() {
  drawProfilePanel(panels['p-temp'], 'temp', '°C', 0);
  drawProfilePanel(panels['p-psal'], 'psal', 'S', 1);
  drawTSPanel(panels['p-ts']);
}

// ---------- legend: VSUP tree ----------

// what the legend has to admit when the uncertainty weight is turned down
const uncNote = () => (UNC < 1
  ? ` Uncertainty weight ${UNC.toFixed(2)}: the tree is collapsed towards its bottom row, so these steps are not what is drawn.`
  : '');
function drawLegend() {
  const d = V();
  document.getElementById('lg-title').textContent = `${d.name} · ${d.unit}`;
  const cv = document.getElementById('lg-vsup');
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const W = Math.round(cv.clientWidth) || 282, H = 118;   // the panel's width; 282 only if hidden
  cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + 'px';
  const ctx = cv.getContext('2d');
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, W, H);
  ctx.font = '9px Cascadia Code, monospace';
  const LW = 74;                                     // left gutter for support labels
  const bw = W - LW - 2, rowH = 15, top = 4;
  if (curVar === 'conf') {                           // support itself: plain bar
    for (let x = 0; x < bw; x++) {
      ctx.fillStyle = css(lutColor(d.lut, x / (bw - 1)));
      ctx.fillRect(LW + x, top, 1, 14);
    }
    ctx.fillStyle = '#2a3038'; ctx.textAlign = 'center';
    for (const v of niceTicks(d.lo, d.hi, 5).ticks)
      ctx.fillText(v.toFixed(1), LW + (v - d.lo) / (d.hi - d.lo) * bw, top + 26);
    document.getElementById('lg-note').textContent =
      'Support = normalised kernel weight (0 far from any cast, 1 on a cast). Geometry proxy only.' + uncNote();
    return;
  }
  document.getElementById('lg-note').textContent = (SHOWCASE
    ? 'A measured profile gets the full colour range. Anything filled in between profiles loses colours as the evidence for it thins — fewer, greyer steps, then nothing.'
    : 'Bottom row = full ramp (measured tubes use it continuously). Interpolated layers climb the tree as support falls — fewer, greyer steps, then nothing.') + uncNote();
  // rows: layer L-1 (1 bin, neutral) at the top ... layer 0 (16 bins) at the bottom
  for (let layer = VSUP_LAYERS - 1; layer >= 0; layer--) {
    const y = top + (VSUP_LAYERS - 1 - layer) * rowH;
    const bins = 1 << (VSUP_LAYERS - 1 - layer);
    const [a, b2] = vsupRange(layer);
    for (let b = 0; b < bins; b++) {
      const v = d.lo + (b + 0.5) / bins * (d.hi - d.lo);
      ctx.fillStyle = css(vsupColor(v, layer === 0 ? 1 : (a + b2) / 2));
      ctx.fillRect(LW + b / bins * bw, y, bw / bins + 0.5, rowH - 1);
    }
    ctx.fillStyle = '#2a3038'; ctx.textAlign = 'right';
    ctx.fillText(layer === 0 ? `≥ ${a.toFixed(2)}` : `${a.toFixed(2)}–${b2.toFixed(2)}`, LW - 5, y + 11);
  }
  ctx.save(); ctx.translate(8, top + VSUP_LAYERS * rowH / 2); ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center'; ctx.fillStyle = '#5b6470'; ctx.fillText('support', 0, 0); ctx.restore();
  const yT = top + VSUP_LAYERS * rowH + 12;
  ctx.fillStyle = '#2a3038'; ctx.textAlign = 'center';
  for (const v of niceTicks(d.lo, d.hi, 5).ticks) {
    const x = LW + (v - d.lo) / (d.hi - d.lo) * bw;
    ctx.fillRect(x, yT - 10, 1, 4);
    ctx.fillText(v.toFixed(d.dec > 1 ? 1 : 0), x, yT + 2);
  }
  ctx.fillText(d.unit, LW + bw / 2, yT + 14);
}

// ---------- time: the field as it was known at an instant ----------
//
// The static grid weights every cast equally, so it merges 30 days into one
// picture of a month that never existed. With the timeline on, each cast's
// kernel weight gains a causal temporal factor
//
//     a_c(p) = exp(-(t - t_c)^2 / 2 L_t(p)^2)
//
// so a profile counts most near its own time and fades either side of it.
// The factor is two-sided on purpose: this record was collected a year ago
// and is complete, so the field at t is a smoother over everything, not a
// filter that has only seen the past. A causal version made every cast pop
// into the field at full weight on the frame it was taken. A real-time
// deployment would have to be causal, and would be worse. L_t(p)
// is fitted from these casts (scripts/fit_time_scale.py): ~6.6 d in the top
// 150 dbar, where repeat visits disagree more the longer the gap; deeper, a
// 30-day record resolves no decay at all, so L_t is pinned at the record
// length and flagged unresolved — a limit of the test, not deep-ocean
// certainty. The field is rebuilt here, in the same separable form used by
// build_demo_data.py, so time weighting cannot drift from the static maths.

// static field, kept pristine so "all June" is a restore, not a recompute
const F0 = { temp: F.temp.slice(), psal: F.psal.slice(), sigma0: F.sigma0.slice(), conf: F.conf.slice() };

let fieldMs = 0;                                     // last rebuild cost, shown in expert mode
function rebuildTimeField() {
  const t0 = performance.now();
  if (!TIME.on) {
    F.temp.set(F0.temp); F.psal.set(F0.psal); F.sigma0.set(F0.sigma0); F.conf.set(F0.conf);
    fieldMs = performance.now() - t0;
    return;
  }
  // these cutoffs are for speed only. Turning the uncertainty weight down is
  // a claim that everything the kernel reaches is worth drawing, so they have
  // to move with it or the drawn area shrinks as the playhead runs.
  const awMin = 1e-5 * Math.max(UNC, 0.01);
  const liveMin = 1e-4 * Math.max(UNC, 0.01);
  accumulateField((c, a) => {
    const dt = TIME.t - c.t;
    // L_t grows with depth, so the deepest level is this cast's best case:
    // if even that has decayed away the cast contributes nothing anywhere
    let live = false;
    for (let k = 0; k < NZ; k++) {
      const w = timeWeight(dt, LT[k] * TIME.scale);
      a[k] = w;
      if (w > liveMin) live = true;
    }
    return live;
  }, awMin);
  fieldMs = performance.now() - t0;
}

// ---------- state / controls ----------

const $ = id => document.getElementById(id);
const fmt = (v, d = 1) => v.toFixed(d);

let exag = 1, exagBefore = null;
function setExag(e) {
  exag = e;
  world.scale.y = e;
  for (const m of markGroup.children) m.position.y = MARK_LIFT / e;
  const s = e < 10 ? fmt(e, 1) : fmt(e, 0);
  $('v-exag').textContent = '×' + s;
  $('exag-val').textContent = s + (e < 1.05 ? ' — true scale' : ' vertical exaggeration');
}
const exagFromSlider = v => Math.exp((v / 1000) * Math.log(500));
$('sl-exag').addEventListener('input', ev => {
  exagBefore = null; $('btn-true').textContent = 'show true scale';
  setExag(exagFromSlider(+ev.target.value));
});
// The OS reduced-motion setting is ignored on purpose: the moves here are the
// content (the fall is how the chart becomes a volume), and with the setting on
// a viewer got the short cue sheet and no transitions, which read as broken.
const reducedMotion = false;
let tween = null;
function toggleTrueScale() {
  const goingTrue = exagBefore === null;
  const target = goingTrue ? 1 : exagBefore;
  exagBefore = goingTrue ? exag : null;
  $('btn-true').textContent = goingTrue ? 'restore ×' + fmt(exagBefore, 0) : 'show true scale';
  if (reducedMotion) { setExag(target); tween = null; return; }
  tween = { from: exag, to: target, t0: performance.now(), ms: 900 };
}
$('btn-true').addEventListener('click', toggleTrueScale);

// ---------- animated view changes ----------
// The exaggeration tween above moves one number; showing the sea floor also has
// to move the camera, because the column becomes 2.7x deeper and the framing
// that suited a 2000-dbar box puts the floor off screen. Both run at once and
// read the live exag, so the target stays on the same water while the scale
// changes under it. The azimuth is left alone: whatever the viewer has rotated
// to is theirs.
let viewTween = null;
function camNow() {
  const t = controls.target;
  return {
    pres: -t.y / exag * 1000,
    hor: Math.hypot(camera.position.x - t.x, camera.position.z - t.z),
    hgt: camera.position.y - t.y,
    az: Math.atan2(camera.position.x - t.x, camera.position.z - t.z),
  };
}
function applyView(v) {
  if (v.fov != null && persp.fov !== v.fov) { persp.fov = v.fov; persp.updateProjectionMatrix(); }
  const y = yOf(v.pres) * exag;
  controls.target.set(0, y, 0);
  camera.position.set(Math.sin(v.az) * v.hor, y + v.hgt, Math.cos(v.az) * v.hor);
  controls.update();
}
// v may omit az, in which case the viewer's own heading is kept: a layer toggle
// has no business turning the scene round. The intro is the one caller that
// asks for a specific one.
function tweenTo(v, e, ms = 1100, force = false) {
  camTween = null;                                 // two things steering the camera is one too many
  if (camera === ortho) goPersp();
  const from = camNow();
  if (v.az == null) v = { ...v, az: from.az };
  // force: the sea floor toggle animates even under prefers-reduced-motion. The
  // move IS the message there — the scale drops and the column shrinks into the
  // top of a much deeper box — and a cut leaves no way to see that happen.
  if (reducedMotion && !force) { setExag(e); applyView(v); tween = viewTween = null; return; }
  tween = { from: exag, to: e, t0: performance.now(), ms };
  viewTween = { from, to: v, t0: performance.now(), ms };
}

// ---------- sea floor toggle ----------
const BATHY_VIEW = { pres: 2900, hor: 2.7, hgt: 0.85, exag: 100 };
let bathyBack = null;
let bathyFade = null;
function setBathy(on) {
  if (on && !bathyGroup.children.length) bathyBuild();
  if (on) bathyGroup.visible = true;                 // hidden again when the fade reaches 0
  bathyFade = { from: bathyMat.opacity, to: on ? 1 : 0, t0: performance.now(), ms: 1100 };
  if (on) {
    bathyBack = { v: camNow(), exag };
    tweenTo({ pres: BATHY_VIEW.pres, hor: BATHY_VIEW.hor * RC, hgt: BATHY_VIEW.hgt * RC },
      BATHY_VIEW.exag, 1100, true);
  } else if (bathyBack) {
    tweenTo(bathyBack.v, bathyBack.exag, 1100, true);
    bathyBack = null;
  }
}
bathyGroup.visible = false;

const castVisible = c => castExists(c);
function onDay(skipInsets = false) {
  let n = 0;
  castMeshes.forEach((m, i) => {
    const vis = castVisible(casts[i]);
    m.visible = vis && $('tg-casts').checked && !introHiding;
    markGroup.children[i].visible = vis && $('tg-casts').checked;
    if (vis) n++;
  });
  syncIsoRings();
  if (!skipInsets) drawInsets();
  return n;
}
$('tg-casts').addEventListener('change', onDay);

function onDepth() {
  const p = (+$('sl-depth').value / 1000) * G.presMax;
  $('v-depth').textContent = fmt(p, 0) + ' dbar';
  updateSlice(p);
}
$('sl-depth').addEventListener('input', onDepth);

let sectAxis = 'lon';
function onSect() {
  const t = +$('sl-sect').value / 1000;
  if (sectAxis === 'lon') {
    const lon = G.lon0 + t * (G.lon1 - G.lon0);
    $('v-sect').textContent = fmt(-lon, 2) + '°W';
    updateSection('lon', lon);
  } else {
    const lat = G.lat0 + t * (G.lat1 - G.lat0);
    $('v-sect').textContent = fmt(lat, 2) + '°N';
    updateSection('lat', lat);
  }
}
$('sl-sect').addEventListener('input', onSect);
for (const b of $('sect-axis').querySelectorAll('button'))
  b.addEventListener('click', () => {
    sectAxis = b.dataset.ax;
    for (const o of $('sect-axis').querySelectorAll('button')) o.classList.toggle('on', o === b);
    onSect();
  });

$('tg-slice').addEventListener('change', () => {
  slice.mesh.visible = $('tg-slice').checked; slice.label.visible = slice.mesh.visible;
  if (slice.mesh.visible && slice.dirty) updateSlice(slice.pres);
});
$('tg-sect').addEventListener('change', () => {
  sect.mesh.visible = $('tg-sect').checked; sect.label.visible = sect.mesh.visible;
  if (sect.mesh.visible && sect.dirty) onSect();
});
// iso-surface: the slider range follows whichever variable is being coloured
const isoRange = () => [VARS[isoVar()].lo, VARS[isoVar()].hi];
function syncIsoUI() {
  const key = isoVar(), [lo, hi] = isoRange();
  ISO.level[key] = Math.max(lo, Math.min(hi, ISO.level[key]));
  $('sl-iso').value = Math.round((ISO.level[key] - lo) / (hi - lo) * 1000);
  $('iso-lab').textContent = ISO_NAME[key][0] + ' (interpolated)';
  $('v-iso').textContent = ISO.level[key].toFixed(isoDec());
}
function rebuildIso() { syncIsoUI(); buildIsoSheets(); buildIsoRings(); }
function setIso(level) { ISO.level[isoVar()] = level; rebuildIso(); }
$('tg-iso').addEventListener('change', () => { ISO.on = $('tg-iso').checked; rebuildIso(); });
$('sl-iso').addEventListener('input', () => {
  const [lo, hi] = isoRange();
  setIso(lo + (+$('sl-iso').value / 1000) * (hi - lo));
});
$('tg-thermo').addEventListener('change', () => { THERMO.on = $('tg-thermo').checked; buildThermo(); });
$('sl-thermo').addEventListener('input', () => {
  THERMO.frac = +$('sl-thermo').value / 100;         // 25..80 -> 0.25..0.80 of peak N^2
  $('v-thermo').textContent = THERMO.frac.toFixed(2);
  buildThermo();
});
$('tg-contours').addEventListener('change', () => { contoursOn = $('tg-contours').checked; refreshField(); });
$('tg-stipple').addEventListener('change', () => { stippleOn = $('tg-stipple').checked; refreshField(); });
$('sl-cut').addEventListener('input', () => {
  CUT = +$('sl-cut').value / 1000; $('v-cut').textContent = CUT.toFixed(2);
  refreshField(); buildIsoSheets(); buildThermo();
});
$('v-cut').textContent = CUT.toFixed(2);
// uncertainty weight: one knob over both halves of the encoding. It changes
// what is drawn, never what is known, so anything below 1 raises a badge.
function syncUncUI() {
  $('sl-unc').value = Math.round(UNC * 1000);
  $('v-unc').textContent = UNC.toFixed(2);
  document.body.classList.toggle('unc-off', UNC < 1);
  $('unc-badge').textContent = UNC < 1
    ? (SHOWCASE ? 'uncertainty hidden \u2014 the same data, drawn as if it were measured everywhere'
      : UNC === 0 ? 'uncertainty weight 0 \u2014 every point any cast reaches, at every instant, drawn as if measured'
        : `uncertainty weight ${UNC.toFixed(2)} \u2014 colour and reach overstate the data`)
    : '';
  $('prov-unc').textContent = UNC < 1
    ? `Currently ${UNC.toFixed(2)}: support is drawn as 1 \u2212 ${UNC.toFixed(2)}(1 \u2212 s). The kernel itself is untouched, so the field still varies in space and time \u2014 only its admission of how little it rests on is suppressed.`
    : '';
}
function setUnc(u) {
  UNC = clamp01(u);
  syncUncUI();
  computeCurves();
  drawLegend();
  recolorCasts();
  if (TIME.on) applyTime();
  else { refreshField(); buildIsoSheets(); buildThermo(); drawTimeline(); }
}
$('sl-unc').addEventListener('input', () => setUnc(+$('sl-unc').value / 1000));
$('sl-rad').addEventListener('input', () => {
  TUBE_R = +$('sl-rad').value / 10; $('v-rad').textContent = TUBE_R.toFixed(1) + ' km'; rebuildCasts(); buildIsoRings();
});
$('v-rad').textContent = TUBE_R.toFixed(1) + ' km';
$('btn-clear').addEventListener('click', () => { pins.length = 0; syncPins(); });

function setVar(key) {
  curVar = key;
  for (const b of $('var-seg').querySelectorAll('button')) b.classList.toggle('on', b.dataset.var === key);
  recolorCasts(); refreshField(); drawLegend();
  rebuildIso();                                      // the surface follows the variable
}
for (const b of $('var-seg').querySelectorAll('button')) b.addEventListener('click', () => setVar(b.dataset.var));
// a hidden plane is not repainted: it is marked and painted when switched on
function refreshField() {
  if (slice.mesh.visible) updateSlice(slice.pres); else slice.dirty = true;
  if (sect.mesh.visible) onSect(); else sect.dirty = true;
}

function setMode(m) {
  document.body.classList.toggle('explain', m === 'explain');
  document.body.classList.toggle('expert', m === 'expert');
  $('mode-explain').classList.toggle('on', m === 'explain');
  $('mode-expert').classList.toggle('on', m === 'expert');
}
$('mode-explain').addEventListener('click', () => setMode('explain'));
$('mode-expert').addEventListener('click', () => setMode('expert'));

// ---------- timeline ----------
//
// Three lanes over the same time axis: when each cast was taken, how much of
// the box is well constrained at two depths as casts arrive and age, and the
// playhead. The two curves are the argument: at 1000 dbar a month of profiles
// accumulates into knowledge, at 50 dbar it does not — surface coverage is a
// treadmill, because L_t there is shorter than the gap between visits.

const CURVE_P = [50, 1000];                          // dbar, one shallow one deep
const CURVE_COL = ['#8a939e', '#1c2128'];
const NT = 160;
const tl = { cv: $('tl'), t: new Float32Array(NT), y: [], stat: [], drag: false };
for (let i = 0; i < NT; i++) tl.t[i] = T.tMin + (T.tMax - T.tMin) * i / (NT - 1);

// fraction of the horizontal plane at pressure p whose support reaches SOLID,
// as a function of playhead time. Uses the same footprints and the same
// normalisation as the field itself.
function computeCurves() {
  tl.y = []; tl.stat = [];
  for (const p of CURVE_P) {
    const k = Math.round(clamp01(p / G.presMax) * (NZ - 1));
    const den = new Float32Array(NX * NY);
    const out = new Float32Array(NT);
    for (let i = 0; i < NT; i++) {
      den.fill(0);
      const tNow = tl.t[i];
      for (let ci = 0; ci < casts.length; ci++) {
        const aw = casts[ci].sv[k] * timeWeight(tNow - casts[ci].t, LT[k] * TIME.scale);
        if (aw < 1e-5) continue;
        const fp = foot[ci];
        for (let y = fp.y0, wi = 0; y <= fp.y1; y++)
          for (let x = fp.x0; x <= fp.x1; x++, wi++) den[y * NX + x] += aw * fp.w[wi];
      }
      let n = 0;
      for (let j = 0; j < den.length; j++) if (den[j] / G.ref >= SOLID) n++;
      out[i] = n / den.length;
    }
    tl.y.push(out);
    // the static all-June field at the same level, for the dashed reference
    let n = 0;
    for (let j = 0; j < NX * NY; j++) if (F0.conf[k * NY * NX + j] >= SOLID) n++;
    tl.stat.push(n / (NX * NY));
  }
}

const TL_H = 96, TL_PADL = 44, TL_PADR = 10;
const tlX = (t, W) => TL_PADL + (t - T.tMin) / (T.tMax - T.tMin) * (W - TL_PADL - TL_PADR);
const tlT = (px, W) => T.tMin + (px - TL_PADL) / (W - TL_PADL - TL_PADR) * (T.tMax - T.tMin);

function drawTimeline() {
  const cv = tl.cv, dpr = Math.min(devicePixelRatio || 1, 2);
  const W = Math.max(200, cv.clientWidth || cv.parentElement.clientWidth - 20);
  if (cv.width !== W * dpr || cv.height !== TL_H * dpr) { cv.width = W * dpr; cv.height = TL_H * dpr; }
  cv.style.height = TL_H + 'px';
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, TL_H);
  ctx.font = '9px Cascadia Code, monospace';
  const yTick0 = 4, yTick1 = 17;                      // cast lane
  const yC0 = 24, yC1 = 74;                           // curve lane
  const yAx = 86;                                     // date axis

  // date axis: a tick every 5 days
  ctx.strokeStyle = '#c4cad1'; ctx.fillStyle = '#5b6470'; ctx.textAlign = 'center';
  ctx.beginPath();
  for (let d = Math.ceil(T.tMin); d <= T.tMax; d++) {
    const x = tlX(d, W), big = d % 5 === 0;
    ctx.moveTo(x + 0.5, yAx - (big ? 6 : 3)); ctx.lineTo(x + 0.5, yAx);
    if (big) ctx.fillText(dayLabel(d).split(' ').slice(0, 2).join(' '), x, yAx + 10);
  }
  ctx.stroke();

  // lane 1: one tick per cast, shaded by how much it counts right now. The
  // weight is two-sided, so casts on both sides of the playhead darken.
  for (const c of casts) {
    const x = tlX(c.t, W);
    const w = TIME.on ? ageWeight(c, 0) : 1;
    ctx.fillStyle = `rgb(${[0.11, 0.13, 0.16].map(v => Math.round(255 + (v * 255 - 255) * Math.max(w, 0.08))).join(',')})`;
    ctx.fillRect(x - 0.75, yTick0, 1.5, yTick1 - yTick0);
  }
  // a pinned cast carries its colour here too. Without it there is no way to
  // find a pinned profile's date among 83 identical ticks, and the date is the
  // whole reason the 2D plots and this lane are worth looking at together.
  // The weight still shades it: being pinned does not make a cast current.
  const markTick = (c, col, half) => {
    const x = tlX(c.t, W);
    ctx.globalAlpha = Math.max(TIME.on ? ageWeight(c, 0) : 1, 0.3);
    ctx.fillStyle = col;
    ctx.fillRect(x - half, yTick0 - 3, 2 * half, yTick1 - yTick0 + 6);
    ctx.globalAlpha = 1;
  };
  pins.forEach((ci, k) => markTick(casts[ci], PIN_COLORS[k], 1.5));
  if (hoveredCast != null && !pins.includes(hoveredCast)) markTick(casts[hoveredCast], '#1c2128', 1.25);

  // lane 2: supported-area curves + dashed static reference
  const maxY = Math.max(0.02, ...tl.stat, ...tl.y.map(a => Math.max(...a))) * 1.15;
  const yv = v => yC1 - clamp01(v / maxY) * (yC1 - yC0);
  ctx.strokeStyle = '#dfe3e8'; ctx.beginPath();
  ctx.moveTo(TL_PADL, yC1 + 0.5); ctx.lineTo(W - TL_PADR, yC1 + 0.5); ctx.stroke();
  CURVE_P.forEach((p, ci) => {
    ctx.setLineDash([2, 3]); ctx.strokeStyle = CURVE_COL[ci]; ctx.globalAlpha = 0.5;
    ctx.beginPath(); ctx.moveTo(TL_PADL, yv(tl.stat[ci])); ctx.lineTo(W - TL_PADR, yv(tl.stat[ci])); ctx.stroke();
    ctx.setLineDash([]); ctx.globalAlpha = 1;
    ctx.strokeStyle = CURVE_COL[ci]; ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (let i = 0; i < NT; i++) {
      const x = tlX(tl.t[i], W), y = yv(tl.y[ci][i]);
      i ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
    }
    ctx.stroke(); ctx.lineWidth = 1;
    ctx.fillStyle = CURVE_COL[ci]; ctx.textAlign = 'left';
    ctx.fillText(`${p} dbar`, W - TL_PADR - 46, yv(tl.y[ci][NT - 1]) - 4);
  });
  ctx.fillStyle = '#5b6470'; ctx.textAlign = 'right';
  ctx.fillText((maxY * 100).toFixed(0) + '%', TL_PADL - 5, yC0 + 4);
  ctx.fillText('0', TL_PADL - 5, yC1 + 3);
  ctx.save();
  ctx.translate(9, (yC0 + yC1) / 2); ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center'; ctx.fillText('area support ≥ 0.5', 0, 0);
  ctx.restore();

  // playhead
  if (TIME.on) {
    const x = tlX(TIME.t, W);
    // shade +-L_t(surface) around the playhead: the window the surface field
    // is actually drawn from, on both sides, not "everything so far"
    const xa = tlX(Math.max(T.tMin, TIME.t - ltAt(0)), W), xb = tlX(Math.min(T.tMax, TIME.t + ltAt(0)), W);
    ctx.fillStyle = 'rgba(180,83,9,0.10)';
    ctx.fillRect(xa, yTick0, xb - xa, yC1 - yTick0);
    ctx.strokeStyle = '#b45309'; ctx.beginPath();
    ctx.moveTo(x + 0.5, yTick0 - 2); ctx.lineTo(x + 0.5, yC1 + 3); ctx.stroke();
    ctx.fillStyle = '#b45309';
    ctx.beginPath(); ctx.moveTo(x, yTick0 - 2); ctx.lineTo(x - 4, yTick0 - 8); ctx.lineTo(x + 4, yTick0 - 8); ctx.fill();
  } else {
    ctx.fillStyle = 'rgba(28,33,40,0.05)';
    ctx.fillRect(TL_PADL, yTick0, W - TL_PADL - TL_PADR, yC1 - yTick0);
  }
}

function timeStatus() {
  if (!TIME.on) {
    $('tl-now').textContent = 'all June (no time weighting)';
    $('tl-note').innerHTML = 'Every cast counts equally, so this is a month of data drawn as one instant &mdash; a snapshot that never existed. Press play.';
    return;
  }
  let live = 0, before = 0;
  for (const c of casts) { if (c.t <= TIME.t) before++; if (ageWeight(c, 0) >= 0.5) live++; }
  $('tl-now').textContent = dayLabel(TIME.t);
  const pct = (a) => (a * 100).toFixed(1);
  const iS = Math.min(NT - 1, Math.round((TIME.t - T.tMin) / (T.tMax - T.tMin) * (NT - 1)));
  if (SHOWCASE) {
    $('tl-note').innerHTML =
      `At this instant, <b>${live}</b> of the ${casts.length} profiles still carry half their weight or more near the surface. ` +
      `The rest were taken too long before or after it to say much about the water now.`;
    return;
  }
  $('tl-note').innerHTML =
    `<b>${live}</b> of ${casts.length} profiles at half weight or better at the surface ` +
    `(${before} of them taken before this instant; the record is complete, so the estimate looks both ways). ` +
    `Area with support &ge; 0.5: <b>${pct(tl.y[0][iS])} %</b> at 50 dbar (all-June ${pct(tl.stat[0])} %), ` +
    `<b>${pct(tl.y[1][iS])} %</b> at 1000 dbar (all-June ${pct(tl.stat[1])} %). ` +
    `L<sub>t</sub> = ${ltAt(0).toFixed(1)} d at the surface, ` +
    `${ltAt(G.presMax).toFixed(0)} d below ${T.rampTo} dbar (unresolved in a 30-day record).` +
    (UNC < 1 ? ` Uncertainty weight ${UNC.toFixed(2)}: the field still moves with the playhead, but the render draws its support as 1 &minus; ${UNC.toFixed(2)}(1 &minus; s).` : '');
}

let thermoAtT = NaN;                                // playhead the layer was last built at
function applyTime(quick = false) {
  const t0 = performance.now();
  rebuildTimeField();
  onDay(quick);                                      // the 2D insets redraw 83
  recolorCasts();                                    // profiles; too slow per frame
  refreshField();
  buildIsoSheets();                                  // one surface, so it can keep up
  // The layer's temporal covariance runs on L_t(800) = the record length, so
  // a frame of playhead (0.04 d at 3 d/s) moves its weights by under 0.1 %.
  // While the playhead runs it is rebuilt every 0.25 d and on every full
  // tick; a rebuild is a third of a tick's cost, and this one draws nothing
  // a frame later could not.
  // While the playhead runs, the full ticks (quick = false, ~6 Hz) are the
  // frames the insets redraw on, and the layer stays off them unless it has
  // fallen 0.5 d behind: the two heaviest rebuilds on one frame is the stutter.
  const lag = Math.abs(TIME.t - thermoAtT);
  const due = TIME.playing
    ? (quick ? !(lag < 0.25) : !(lag < 0.5))
    : true;
  if (due) { thermoAtT = TIME.t; buildThermo(); }
  drawTimeline();
  timeStatus();
  $('tl-ms').textContent = `field ${fieldMs.toFixed(0)} / all ${(performance.now() - t0).toFixed(0)} ms`;
}
function setTime(t, quick = false) {
  TIME.t = Math.max(T.tMin, Math.min(T.tMax, t));
  applyTime(quick);
}
function setTimeOn(on) {
  TIME.on = on;
  $('btn-alltime').textContent = on ? 'show all June' : 'back to timeline';
  $('btn-play').disabled = !on;
  document.body.classList.toggle('timed', on);
  applyTime();
}
$('btn-alltime').addEventListener('click', () => { TIME.playing = false; syncPlayBtn(); setTimeOn(!TIME.on); });
function syncPlayBtn() { $('btn-play').textContent = TIME.playing ? '❚❚ pause' : '▶ play'; }
$('btn-play').addEventListener('click', () => {
  if (!TIME.on) setTimeOn(true);
  TIME.playing = !TIME.playing;
  if (TIME.playing && TIME.t >= T.tMax - 1e-6) setTime(T.tMin, true);
  else if (!TIME.playing) applyTime();               // a paused frame is a full tick: nothing left stale
  syncPlayBtn();
});
$('sl-lt').addEventListener('input', () => {
  TIME.scale = Math.exp((+$('sl-lt').value / 1000 * 2 - 1) * Math.log(4));   // x1/4 .. x4
  $('v-lt').textContent = '×' + TIME.scale.toFixed(2);
  syncUncUI();
  computeCurves();
  if (TIME.on) applyTime(); else drawTimeline();
});
{
  const pos = ev => {
    const r = tl.cv.getBoundingClientRect();
    return tlT((ev.clientX - r.left) * (tl.cv.clientWidth / r.width), tl.cv.clientWidth);
  };
  tl.cv.addEventListener('pointerdown', ev => {
    tl.drag = true; tl.cv.setPointerCapture(ev.pointerId);
    TIME.playing = false; syncPlayBtn();
    if (!TIME.on) setTimeOn(true);
    setTime(pos(ev), true);
  });
  tl.cv.addEventListener('pointermove', ev => { if (tl.drag) setTime(pos(ev), true); });
  tl.cv.addEventListener('pointerup', ev => { tl.drag = false; applyTime(); });
}

// ---------- view gizmo: E / N / Up balls, click to snap; centre = isometric ----------

const gizmo = $('gizmo');
const GZ = 96, GR = 34, GB = 9;
gizmo.width = GZ * DPR2; gizmo.height = GZ * DPR2;
gizmo.style.width = gizmo.style.height = GZ + 'px';
const AXES = [
  { n: 'E', dir: new THREE.Vector3(1, 0, 0), col: '#c94b4b' },
  { n: 'N', dir: new THREE.Vector3(0, 0, -1), col: '#3f8f4f' },
  { n: 'U', dir: new THREE.Vector3(0, 1, 0), col: '#3f6fbf' },
];
const _q = new THREE.Quaternion(), _d = new THREE.Vector3();
let gizmoHits = [];
function drawGizmo() {
  if (SHOWCASE) return;                            // display:none there; the compass has the cell
  const ctx = gizmo.getContext('2d');
  ctx.setTransform(DPR2, 0, 0, DPR2, 0, 0);
  ctx.clearRect(0, 0, GZ, GZ);
  ctx.fillStyle = 'rgba(255,255,255,0.75)';
  ctx.beginPath(); ctx.arc(GZ / 2, GZ / 2, GZ / 2 - 1, 0, Math.PI * 2); ctx.fill();
  _q.copy(camera.quaternion).invert();
  const items = [];
  for (const a of AXES) for (const sgn of [1, -1]) {
    _d.copy(a.dir).multiplyScalar(sgn).applyQuaternion(_q);
    items.push({ a, sgn, x: GZ / 2 + _d.x * GR, y: GZ / 2 - _d.y * GR, z: _d.z });
  }
  items.sort((p, q) => p.z - q.z);                 // back to front
  gizmoHits = items;
  ctx.font = `bold ${10}px Segoe UI, sans-serif`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  for (const it of items) {
    const dim = 0.55 + 0.45 * (it.z + 1) / 2;      // fade the far side
    ctx.globalAlpha = dim;
    if (it.sgn > 0) {
      ctx.strokeStyle = it.a.col; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(GZ / 2, GZ / 2); ctx.lineTo(it.x, it.y); ctx.stroke();
      ctx.fillStyle = it.a.col;
      ctx.beginPath(); ctx.arc(it.x, it.y, GB, 0, Math.PI * 2); ctx.fill();
      ctx.fillStyle = '#fff'; ctx.fillText(it.a.n, it.x, it.y + 0.5);
    } else {
      ctx.fillStyle = 'rgba(255,255,255,0.9)'; ctx.strokeStyle = it.a.col; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(it.x, it.y, GB - 1, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    }
  }
  ctx.globalAlpha = 1;
  ctx.fillStyle = '#5b6470';
  ctx.beginPath(); ctx.arc(GZ / 2, GZ / 2, 2.5, 0, Math.PI * 2); ctx.fill();
}
let camTween = null;
function snapView(dir) {
  if (camera === ortho) goPersp();
  // keep the orbit distance, move the camera onto the axis through the target
  const off = camera.position.clone().sub(controls.target);
  const r = off.length();
  const from = off.clone().normalize();
  const to = dir.clone().normalize();
  if (Math.abs(to.y) > 0.999) to.x = to.z = 0.02 * Math.sign(to.y || 1);   // keep OrbitControls' polar clamp happy
  to.normalize();
  if (reducedMotion) { camera.position.copy(controls.target).addScaledVector(to, r); controls.update(); goOrtho(to); return; }
  camTween = { from, to, r, r1: r, t0: performance.now(), ms: 450, then: () => goOrtho(to) };
}
// Snapped views are orthographic (as in Blender); any rotation away returns to
// perspective. Frustum height is matched at the orbit distance so the swap is
// visually continuous; OrbitControls drives ortho.zoom for the wheel.
let snapDir = null;
let compassBack = null;                            // the orbit the map view was entered from
function goOrtho(dir) {
  const r = camera.position.distanceTo(controls.target);
  const halfH = r * Math.tan(THREE.MathUtils.degToRad(persp.fov / 2));
  ortho.top = halfH; ortho.bottom = -halfH;
  ortho.left = -halfH * persp.aspect; ortho.right = halfH * persp.aspect;
  ortho.zoom = 1;
  ortho.position.copy(persp.position); ortho.quaternion.copy(persp.quaternion);
  ortho.updateProjectionMatrix();
  snapDir = dir.clone();                           // before update(): it fires 'change'
  camera = ortho; controls.object = ortho; controls.update();
}
function goPersp() {
  const dir = ortho.position.clone().sub(controls.target).normalize();
  const halfH = ortho.top / ortho.zoom;            // current visible half-height
  const r = halfH / Math.tan(THREE.MathUtils.degToRad(persp.fov / 2));
  persp.position.copy(controls.target).addScaledVector(dir, r);
  snapDir = null;
  compassBack = null;                              // any way out of the map view is a way out
  camera = persp; controls.object = persp; controls.update();
}
controls.addEventListener('change', () => {
  if (camera !== ortho || camTween || !snapDir) return;
  const dir = ortho.position.clone().sub(controls.target).normalize();
  if (dir.angleTo(snapDir) > 0.02) goPersp();     // rotated off the snapped axis
});
gizmo.addEventListener('pointerdown', ev => {
  const rect = gizmo.getBoundingClientRect();
  const x = ev.clientX - rect.left, y = ev.clientY - rect.top;
  let best = null;
  for (const it of gizmoHits) {                    // front-most wins: iterate back->front, keep last hit
    if (Math.hypot(it.x - x, it.y - y) <= GB + 1) best = it;
  }
  if (best) snapView(best.a.dir.clone().multiplyScalar(best.sgn));
  else if (Math.hypot(x - GZ / 2, y - GZ / 2) <= 12) snapView(new THREE.Vector3(1, 1, 1));   // isometric
});
gizmo.addEventListener('pointermove', ev => {
  const rect = gizmo.getBoundingClientRect();
  const x = ev.clientX - rect.left, y = ev.clientY - rect.top;
  const hit = gizmoHits.some(it => Math.hypot(it.x - x, it.y - y) <= GB + 1) || Math.hypot(x - GZ / 2, y - GZ / 2) <= 12;
  gizmo.style.cursor = hit ? 'pointer' : 'default';
});

// ---------- compass (showcase) ----------
// The showcase has no axis gizmo: it is a demonstration, not an instrument,
// and a set of labelled axis balls is a modelling control. It still turns,
// though, and a viewer who has turned it has no way of telling which way
// north is or of getting back. So: a compass rose lying flat on the sea
// surface, drawn through the live camera rather than as a 2D dial, which is
// what makes it read the tilt as well as the heading -- a circle seen from
// overhead, a line seen from the side, and the ellipse between. Clicking it
// flies to the map view the opening holds on; clicking it again returns the
// orbit it left.
const compass = $('compass');
const CS = 84, CR = CS * 0.30;                     // canvas units; style.css sizes the box
compass.width = compass.height = CS * DPR2;
const C_INK = '28,33,40', C_ACC = '180,83,9';      // --ink, --accent
let compassOver = false;
const _cw = new THREE.Vector3(), _cz = new THREE.Vector3(), _cq = new THREE.Quaternion();
// a point of the ground plane at bearing `a` (0 = north, + towards east) and
// radius `r` in ring units, projected through the camera into canvas px; the
// third component is its depth, positive towards the viewer
function cpt(a, r) {
  _cw.set(Math.sin(a) * r, 0, -Math.cos(a) * r).applyQuaternion(_cq);
  return [CS / 2 + _cw.x * CR, CS / 2 - _cw.y * CR, _cw.z];
}
function drawCompass() {
  if (!SHOWCASE || !compass.classList.contains('in')) return;   // nothing to draw before the step that reveals it
  const ctx = compass.getContext('2d');
  ctx.setTransform(DPR2, 0, 0, DPR2, 0, 0);
  ctx.clearRect(0, 0, CS, CS);
  const c = CS / 2;
  const lift = compassOver ? 1.25 : 1;             // no plate to brighten on hover; the ink takes it
  _cq.copy(camera.quaternion).invert();

  // No plate. The rose is drawn straight onto the scene, which is the whole
  // point of a compass that lies in it -- but the ground is light and a tube
  // is not, so every stroke carries a tight white halo instead. That is a
  // legibility floor, not a background: it follows the ink rather than
  // boxing it, so nothing is hidden behind the widget.
  ctx.shadowColor = `rgba(255,255,255,${compassOver ? 0.95 : 0.8})`;
  ctx.shadowBlur = 3;

  // the ring, split at the horizon of its own plane. The depth of a
  // ground-plane point at bearing a is sin(a)*Ez + cos(a)*Nz, zero at
  // atan2(-Nz, Ez) and positive over the half turn above it, and drawing that
  // half stronger is the only depth cue a flat ellipse has. How much stronger
  // is `tilt` = hypot(Ez, Nz) = cos(elevation): overhead the two halves are
  // the same ring and any split would be an artefact of a rounding error.
  _cz.set(1, 0, 0).applyQuaternion(_cq); const ez = _cz.z;
  _cz.set(0, 0, -1).applyQuaternion(_cq); const nz = _cz.z;
  const horizon = Math.atan2(-nz, ez);
  const tilt = Math.min(1, Math.hypot(ez, nz));
  const arc = (from, to, alpha, w) => {
    ctx.beginPath();
    for (let i = 0; i <= 48; i++) {
      const p = cpt(from + (to - from) * i / 48, 1);
      if (i) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]);
    }
    ctx.strokeStyle = `rgba(${C_INK},${Math.min(1, alpha * lift).toFixed(3)})`; ctx.lineWidth = w; ctx.stroke();
  };
  arc(horizon - Math.PI, horizon, 0.44 - 0.28 * tilt, 1.3 - 0.3 * tilt);
  arc(horizon, horizon + Math.PI, 0.44, 1.3);

  // E, S, W as ticks; north is the needle
  for (const b of [Math.PI / 2, Math.PI, -Math.PI / 2]) {
    const p = cpt(b, 1);
    ctx.fillStyle = `rgba(${C_INK},${Math.min(1, (0.42 - 0.28 * Math.max(0, -p[2])) * lift).toFixed(3)})`;
    ctx.beginPath(); ctx.arc(p[0], p[1], 1.7, 0, Math.PI * 2); ctx.fill();
  }

  // the needle is two triangles in the ground plane, so it lies down with the
  // rose instead of standing up off it. Tail first: on an edge-on view the
  // two overlap and the head is the half that has to survive.
  const b1 = cpt(Math.PI / 2, 0.2), b2 = cpt(-Math.PI / 2, 0.2);
  const tri = (p, fill) => {
    ctx.beginPath(); ctx.moveTo(p[0], p[1]); ctx.lineTo(b1[0], b1[1]); ctx.lineTo(b2[0], b2[1]);
    ctx.closePath(); ctx.fillStyle = fill; ctx.fill();
  };
  tri(cpt(Math.PI, 0.6), `rgba(${C_INK},0.28)`);
  tri(cpt(0, 0.94), `rgba(${C_ACC},0.92)`);
  ctx.fillStyle = `rgba(${C_INK},0.55)`;
  ctx.beginPath(); ctx.arc(c, c, 1.7, 0, Math.PI * 2); ctx.fill();

  // the letter stays upright: it is a label, not part of the rose
  const lp = cpt(0, 1.26);
  ctx.font = '600 9.5px Segoe UI, sans-serif';
  ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillStyle = `rgba(${C_INK},0.72)`;
  ctx.fillText('N', lp[0], lp[1]);

  // in the map view the rose is a circle and the needle is straight up, which
  // says nothing about there being a way back: the rim does
  if (compassBack) {
    ctx.strokeStyle = `rgba(${C_ACC},0.55)`; ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.arc(c, c, c - 1.6, 0, Math.PI * 2); ctx.stroke();
  }
}
// The map view is the one the opening holds on: straight down, north up, the
// domain circle filling the frame, orthographic on arrival because a map has
// no vanishing point. Built from the same two numbers introStart() uses, so
// the two cannot drift apart. The exaggeration and the orbit centre are left
// alone -- a compass turns the camera, it does not rescale the scene.
function flyTo(to, r1, then) {
  const off = camera.position.clone().sub(controls.target);
  const r = off.length();
  if (reducedMotion) {
    camera.position.copy(controls.target).addScaledVector(to, r1);
    controls.update(); then?.(); return;
  }
  camTween = { from: off.normalize(), to: to.clone(), r, r1, t0: performance.now(), ms: 900, then };
}
function compassToggle() {
  const back = compassBack;                        // read before goPersp(), which clears it
  compassBack = null;
  if (camera === ortho) goPersp();                 // the flight is a perspective orbit either way
  if (back) { compass.title = 'map view'; flyTo(back.dir, back.r, null); return; }
  const off = camera.position.clone().sub(controls.target);
  compassBack = { dir: off.clone().normalize(), r: off.length() };
  compass.title = 'back to the last view';
  // 0.005 rad off vertical, as in introStart(): straight up is OrbitControls'
  // polar limit and the azimuth stops meaning anything there
  const to = new THREE.Vector3(0, 1, 0.005).normalize();
  flyTo(to, H0 / Math.tan(THREE.MathUtils.degToRad(FOV1 / 2)), () => goOrtho(to));
}
// Where it hangs is a measurement. The compass is meant to sit on the top edge
// of whatever holds the bottom right -- the timeline, or the colour ramp once a
// short viewport moves it there -- and that edge is set by the height of the
// plots and by the breakpoint, so there is no constant to write. It shares the
// guide's cell (card pinned to the top, this to the bottom) and a negative
// bottom margin drops it the rest of the way; the correction is relative, so it
// is a no-op once it is right and converges in one pass from anywhere.
// Positions are read with offsetTop, not getBoundingClientRect: the panels are
// mid-reveal transform for the first 600 ms of their life and offsets are not.
const C_HANG = 8;                                  // clear of the panel it hangs on
function placeCompass() {
  if (!SHOWCASE) return;
  if (MOBILE) { compass.style.marginBottom = ''; return; }
  let top = innerHeight;
  for (const id of ['timeline', 'legend']) {
    const el = $(id);
    if (!el || getComputedStyle(el).display === 'none') continue;
    if (el.offsetLeft + el.offsetWidth <= compass.offsetLeft) continue;   // a different column
    top = Math.min(top, el.offsetTop);
  }
  const cur = parseFloat(compass.style.marginBottom) || 0;
  const delta = (top - C_HANG) - (compass.offsetTop + compass.offsetHeight);
  compass.style.marginBottom = (cur - delta) + 'px';
}
compass.addEventListener('click', () => { if (controls.enabled) compassToggle(); });
compass.addEventListener('pointerenter', () => { compassOver = true; });
compass.addEventListener('pointerleave', () => { compassOver = false; });

// ---------- hover / click ----------

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();
const readout = $('readout');
// the readout keeps to the viewport: right of the pointer while there is room,
// else left of it, and the same below/above. Its content (so its size) is set
// after the position, so frame() re-runs this before the paint.
const readoutAt = { x: 0, y: 0 };
function placeReadout() {
  if (readout.style.display !== 'block') return;
  const w = readout.offsetWidth, h = readout.offsetHeight;
  let x = readoutAt.x + 14, y = readoutAt.y + 10;
  if (x + w > innerWidth - 6) x = Math.max(6, readoutAt.x - 14 - w);
  if (y + h > innerHeight - 6) y = Math.max(6, readoutAt.y - 10 - h);
  readout.style.left = x + 'px'; readout.style.top = y + 'px';
}
const kv = (k, v) => `<span class="k">${k}</span> ${v}`;
function setHover(ci) {
  if (hoveredCast === ci) return;
  const was = hoveredCast;
  hoveredCast = ci;
  if (was != null) applyHilite(was);
  if (ci != null) applyHilite(ci);
  drawInsets();
  drawTimeline();
}

function pickAt(ev) {
  pointer.x = (ev.clientX / innerWidth) * 2 - 1;
  pointer.y = -(ev.clientY / innerHeight) * 2 + 1;
  raycaster.setFromCamera(pointer, camera);
  const tubes = raycaster.intersectObjects(castMeshes.filter(m => m.visible));
  if (tubes.length) return { kind: 'cast', hit: tubes[0] };
  const planes = [slice.mesh, sect.mesh].filter(m => m.visible);
  const ph = raycaster.intersectObjects(planes);
  for (const h of ph) {
    // reject hits on transparent (unsupported) parts: sample the field
    const p = h.point, lon = lonOfX(p.x), lat = latOfZ(p.z), pres = -p.y / exag * 1000;
    const s = fieldAt(lon, lat, pres);
    if (s && s.conf >= CUT) return { kind: 'field', hit: h, lon, lat, pres, s };
  }
  return null;
}
function syncPins() {
  castMeshes.forEach((m, ci) => applyHilite(ci));
  drawTimeline();
  for (const L of pinLabels) { L.el.remove(); labels.splice(labels.indexOf(L), 1); }
  pinLabels.length = 0;
  pins.forEach((ci, k) => {
    const c = casts[ci];
    const L = addLabel(`▮ ${c.id}/${c.cyc}`, xOf(c.lon), 0.05, zOf(c.lat), 'pin');
    L.el.style.color = PIN_COLORS[k];
    pinLabels.push(L);
  });
  drawInsets();
}
// Filling the readout is its own function because a mouse raises it by
// hovering and a finger has to raise it by tapping — same content, two events.
function showReadout(ev) {
  const pk = pickAt(ev);
  if (!pk) {
    readout.style.display = 'none';
    setHover(null);
    return false;
  }
  readout.style.display = 'block';
  readoutAt.x = ev.clientX; readoutAt.y = ev.clientY;
  placeReadout();
  if (pk.kind === 'cast') {
    const ci = pk.hit.object.userData.ci, c = casts[ci];
    const pres = -pk.hit.point.y / exag * 1000;
    let bi = 0;
    while (bi < c.pres.length - 1 && c.pres[bi + 1] < pres) bi++;
    const age = TIME.t - c.t;
    // showcase: what it is, where on it the cursor is, and how far the
    // playhead is from it. The rest of the record is in the 2D plots; a
    // tooltip that follows the cursor is the wrong place for it.
    if (SHOWCASE) {
      const ag = Math.abs(age).toFixed(1);
      readout.innerHTML =
        `<b>MEASURED</b> float ${c.id} · ${c.stamp}\n` +
        kv('at', `${fmt(c.pres[bi], 0)} m`) + kv('  T', `${fmt(c.temp[bi], 2)} °C`) +
        (TIME.on ? `\n${age >= 0 ? ag + ' d before' : ag + ' d after'} the playhead` : '');
      setHover(ci);
      return;
    }
    readout.innerHTML =
      `<b>MEASURED</b> float ${c.id} cycle ${c.cyc} · ${c.stamp}\n` +
      (TIME.on ? kv('\u0394t', `${age >= 0 ? '+' : '\u2212'}${Math.abs(age).toFixed(1)} d ${age >= 0 ? 'before' : 'after'} playhead`) +
        kv('  counts', `${ageWeight(c, pres).toFixed(2)} (L_t ${ltAt(pres).toFixed(1)} d here)`) + '\n' : '') +
      kv('mode', modeName(c.mode)) + kv('  σ', `±${castSigma(c)} °C`) + '\n' +
      kv('pres', `${fmt(c.pres[bi], 0)} dbar`) + kv('  lat/lon', `${c.lat.toFixed(2)}, ${c.lon.toFixed(2)}`) + '\n' +
      kv('T', `${fmt(c.temp[bi], 2)} °C`) + kv('  S', fmt(c.psal[bi], 3)) + kv('  σ₀', `${fmt(c.sigma0[bi], 3)}`);
    setHover(ci);
  } else {
    const s = pk.s;
    const sup = s.conf >= SOLID ? 'full palette' : `VSUP layer ${vsupLayer(s.conf)} of ${VSUP_LAYERS - 1}`;
    readout.innerHTML =
      `<b>INTERPOLATED</b> (kernel field, not a measurement)\n` +
      kv('at', `${pk.lat.toFixed(2)}°N ${(-pk.lon).toFixed(2)}°W ${fmt(pk.pres, 0)} dbar`) + '\n' +
      kv('T', `${fmt(s.temp, 2)} °C`) + kv('  S', fmt(s.psal, 3)) + kv('  σ₀', fmt(s.sigma0, 3)) + '\n' +
      kv('support', `${s.conf.toFixed(2)} (${sup})`);
    setHover(null);
  }
  return true;
}
const hideReadout = () => { readout.style.display = 'none'; setHover(null); };
// A touch drag fires pointermove throughout, so hovering off it would raycast
// 83 tubes on every frame of every orbit — the one gesture a phone can least
// afford. Touch raises the readout on the tap instead (below).
renderer.domElement.addEventListener('pointermove', ev => {
  if (ev.pointerType === 'touch') return;
  showReadout(ev);
});
// the controls step clears itself when both gestures have actually been used
renderer.domElement.addEventListener('pointermove', ev => { if (ev.buttons) navDrag = true; });
renderer.domElement.addEventListener('wheel', () => { navZoom = true; }, { passive: true });
// pinch is the touch dolly and fires no wheel event, so the zoom half of that
// step is read off the camera distance instead of off the gesture
let navDist = null;
controls.addEventListener('change', () => {
  const d = controls.getDistance();
  if (navDist == null) { navDist = d; return; }
  if (Math.abs(d - navDist) / navDist > 0.08) { navZoom = true; navDist = d; }
});
let downPos = null;
renderer.domElement.addEventListener('pointerdown', ev => {
  downPos = [ev.clientX, ev.clientY];
  if (ev.pointerType === 'touch') hideReadout();
});
// A drag on the scene while a flight is in the air is the viewer taking the
// camera back. The flight stops where it is rather than fighting the orbit for
// the rest of its 900 ms, and since nothing has been snapped yet the camera is
// still perspective, so the drag simply continues from wherever the tween had
// got to -- no cut, and no lerp pulling against the pointer. If it had already
// landed, the map view is orthographic, and this leaves it for perspective on
// the gesture rather than after the 0.02 rad the 'change' handler waits for,
// so there is no moment where a drag is turning an orthographic camera.
// Threshold, not the first move: a click that pins a cast is not a rotation.
renderer.domElement.addEventListener('pointermove', ev => {
  if (!downPos || !ev.buttons || (!camTween && camera !== ortho)) return;
  if (Math.hypot(ev.clientX - downPos[0], ev.clientY - downPos[1]) < 3) return;
  camTween = null;
  compassBack = null;
  if (camera === ortho) goPersp();
});
renderer.domElement.addEventListener('pointerup', ev => {
  if (!downPos) return;
  const moved = Math.hypot(ev.clientX - downPos[0], ev.clientY - downPos[1]);
  downPos = null;
  // a finger is wider and shakier than a cursor: 5 px of slop rejects taps a
  // mouse would have made cleanly
  if (moved > (ev.pointerType === 'touch' ? 12 : 5)) return;
  // one tap does both jobs on touch, because there is no second gesture to
  // give the readout: it pins the cast and reads out where the tap landed
  if (ev.pointerType === 'touch') showReadout(ev);
  const pk = pickAt(ev);
  if (pk && pk.kind === 'cast') {
    const ci = pk.hit.object.userData.ci;
    const at = pins.indexOf(ci);
    if (at >= 0) pins.splice(at, 1); else { if (pins.length >= 3) pins.shift(); pins.push(ci); }
    syncPins();
  }
});

// ---------- guided tour: captions computed from the data at load ----------

function tourFacts() {
  // near-surface and 500 dbar temperature across casts
  const near = [], deep = [];
  for (const c of casts) {
    if (c.pres[0] <= 15) near.push(c.temp[0]);
    let bi = 0; while (bi < c.pres.length - 1 && c.pres[bi + 1] < 500) bi++;
    if (c.pres[c.pres.length - 1] >= 500) deep.push(c.temp[bi]);
  }
  const mm = a => [Math.min(...a), Math.max(...a)];
  // strongest horizontal contrast at 300 dbar between casts < 150 km apart
  const t300 = casts.map(c => {
    if (c.pres[c.pres.length - 1] < 300) return null;
    let bi = 0; while (bi < c.pres.length - 1 && c.pres[bi + 1] < 300) bi++;
    return c.temp[bi];
  });
  let front = { dT: 0 };
  for (let i = 0; i < casts.length; i++) for (let j = i + 1; j < casts.length; j++) {
    if (t300[i] == null || t300[j] == null) continue;
    const km = Math.hypot(xOf(casts[i].lon) - xOf(casts[j].lon), zOf(casts[i].lat) - zOf(casts[j].lat));
    if (km > 150 || km < 5) continue;
    const dT = Math.abs(t300[i] - t300[j]);
    if (dT > front.dT) front = { dT, km, i, j };
  }
  // sigma0 = 27.0 depth per cast (first crossing from the top)
  const d27 = [];
  for (const c of casts) for (let i = 1; i < c.pres.length; i++)
    if ((c.sigma0[i - 1] - 27) * (c.sigma0[i] - 27) <= 0) {
      d27.push(c.pres[i - 1] + (27 - c.sigma0[i - 1]) / (c.sigma0[i] - c.sigma0[i - 1]) * (c.pres[i] - c.pres[i - 1]));
      break;
    }
  let nSolid = 0, nAny = 0;
  for (let i = 0; i < F.conf.length; i++) { if (F.conf[i] >= SOLID) nSolid++; if (F.conf[i] >= 0.08) nAny++; }
  const modes = { D: 0, A: 0, R: 0 };
  for (const c of casts) modes[c.mode] = (modes[c.mode] || 0) + 1;
  // temporal structure function, the same quantity scripts/fit_time_scale.py
  // fits: how far apart two casts of the same water end up as the gap grows.
  // Reported for a shallow band (short L_t) and a deep one (none resolvable).
  const bandMean = (c, p0, p1) => {
    let s = 0, n = 0;
    for (let i = 0; i < c.pres.length; i++)
      if (c.pres[i] >= p0 && c.pres[i] < p1) { s += c.temp[i]; n++; }
    return n ? s / n : null;
  };
  const lagRms = (p0, p1) => {
    const v = casts.map(c => bandMean(c, p0, p1));
    const acc = [[0, 0], [0, 0]];                   // [near lag, far lag] -> [sum sq, n]
    for (let i = 0; i < casts.length; i++) for (let j = i + 1; j < casts.length; j++) {
      if (v[i] == null || v[j] == null) continue;
      const km = Math.hypot(xOf(casts[i].lon) - xOf(casts[j].lon), zOf(casts[i].lat) - zOf(casts[j].lat));
      if (km > 60) continue;
      const lag = Math.abs(casts[i].t - casts[j].t);
      const b = lag < 2 ? 0 : (lag >= 10 ? 1 : -1);
      if (b < 0) continue;
      acc[b][0] += (v[i] - v[j]) ** 2; acc[b][1]++;
    }
    return acc.map(([s, n]) => (n ? Math.sqrt(s / n) : NaN));
  };
  const [shNear, shFar] = lagRms(0, 150);
  // the deep claim is the permutation test from the fit, not a second eyeball
  // comparison: a two-bin rms over 300-1000 dbar spans too steep a gradient to
  // mean anything on its own
  const deepBand = T.bands.find(b => b.p0 === 300) || T.bands[T.bands.length - 1];
  return {
    near: mm(near), deep: mm(deep), front, d27: mm(d27), nSolid: nSolid / F.conf.length,
    nAny: nAny / F.conf.length, modes, shNear, shFar, deepBand
  };
}
const FACTS = tourFacts();
const f1 = v => v.toFixed(1), f0 = v => v.toFixed(0);
const setLayers = o => {
  for (const [id, on] of Object.entries(o)) {
    const el = $(id);
    if (el.checked !== on) { el.checked = on; el.dispatchEvent(new Event('change')); }
  }
};
const TOUR = [
  {
    name: 'Profiles', apply() {
      setVar('temp'); setLayers({ 'tg-casts': true, 'tg-slice': false, 'tg-sect': false, 'tg-iso': false, 'tg-thermo': false });
      setView('overview');
    }, cap: () => `Each tube is one Argo profile: a float rising from ~2000 dbar to the surface, measuring every metre or two. ` +
      `Colour along the tube is <b>temperature</b>. Surface water here spans <b>${f1(FACTS.near[0])}–${f1(FACTS.near[1])} °C</b>; ` +
      `at 500 dbar it is <b>${f1(FACTS.deep[0])}–${f1(FACTS.deep[1])} °C</b>. Tubes are ${f0(TUBE_R * 2)} km wide on the map ` +
      `but the instrument samples a column centimetres across — horizontally the data is almost nothing.` +
      `<span class="n">${FACTS.modes.D} profiles are delayed-mode calibrated (filled marker), ${FACTS.modes.R} are real-time and uncalibrated (ring marker), ${FACTS.modes.A} adjusted.</span>`
  },
  {
    name: 'Depth slice', apply() {
      setVar('temp'); setLayers({ 'tg-casts': true, 'tg-slice': true, 'tg-sect': false, 'tg-iso': false, 'tg-thermo': false });
      $('sl-depth').value = 150; onDepth(); setView('overview');
    }, cap: () => `A horizontal cut at <b>${f0(slice.pres)} dbar</b> through the interpolated field, with contour lines every ${V().contour} ${V().unit}. ` +
      `White dots are the profiles that pass this depth — the only places the colour is a measurement. Between them a kernel ` +
      `(25 km horizontal, 30 dbar vertical) fills in; as its support drops the colours collapse to <b>fewer, greyer steps</b> (the legend tree) and the layer fades out; far from any profile nothing is left. Contour lines follow the same fade — over grey they are the kernel's guess, not data. Drag the depth slider.`
  },
  {
    name: 'The front', apply() {
      setVar('temp'); setLayers({ 'tg-casts': false, 'tg-slice': false, 'tg-sect': true, 'tg-iso': false, 'tg-thermo': false });
      const fr = FACTS.front;
      sectAxis = 'lon';
      for (const o of $('sect-axis').querySelectorAll('button')) o.classList.toggle('on', o.dataset.ax === 'lon');
      const lon = fr.i != null ? (casts[fr.i].lon + casts[fr.j].lon) / 2 : lonMid;
      $('sl-sect').value = Math.round((lon - G.lon0) / (G.lon1 - G.lon0) * 1000); onSect();
      setView('section');
    }, cap: () => {
      const fr = FACTS.front; return `A north–south vertical section through the interpolated field, contours every 2 °C. The Gulf Stream front crosses this box, warm Sargasso water to the south, cold slope water to the north: ` +
        `at 300 dbar the sharpest measured contrast is <b>${f1(fr.dT)} °C over ${f0(fr.km)} km</b> (floats ${casts[fr.i].id} and ${casts[fr.j].id}). ` +
        `Tubes are hidden; dashed dark lines mark the profiles within 40 km of the plane, and the section is only trustworthy near them. Switch to E–W, move it with the slider, or turn the profiles back on.`;
    }
  },
  {
    name: 'Density layers', apply() {
      setVar('sigma0'); setLayers({ 'tg-casts': true, 'tg-slice': false, 'tg-sect': false, 'tg-iso': true, 'tg-thermo': false });
      setIso(27.0); setView('density');
    }, cap: () => `Water is stratified by density; flow follows surfaces of constant density (isopycnals) rather than depth. ` +
      `The collar on each tube is where that profile crosses <b>σ₀ = 27.0 kg m⁻³</b> — measured, between <b>${f0(FACTS.d27[0])} and ${f0(FACTS.d27[1])} dbar</b>: shallow on the cold side of the front, deep under the warm Sargasso Sea. ` +
      `The gridded sheet between collars is the same surface in the interpolated field; it greys and fades as support falls, by the same rule as the slices. ` +
      `<span class="n">Drag the value slider to sweep the surface through the water column — one surface at a time, so you can watch it move rather than read five fixed sheets. It follows whichever variable is coloured: an isotherm for temperature, an isohaline for salinity.</span>`
  },
  {
    name: 'Main pycnocline', apply() {
      setVar('temp'); setLayers({ 'tg-casts': true, 'tg-slice': false, 'tg-sect': true, 'tg-iso': false, 'tg-thermo': true });
      setView('density');
    }, cap: () => {
      const s = thermoStats(); return `Some structures are layers, not surfaces, and are drawn as volumes. This one is the <b>permanent pycnocline</b> — the main thermocline of the subtropical gyre, not the shallow seasonal one. These profiles have both: a sharp seasonal layer in the top ~100 dbar holding 22 % of the temperature drop above 1000 dbar, and this one, three times weaker per dbar but five times thicker, holding 59 % of it. ` +
        `It is found the way Feucher et al. (2016, 2019) find it, by the <b>order of the stratification extrema</b> rather than by a threshold: the seasonal N² peak, then the first N² minimum below it (the subtropical mode water), then this — the N² maximum below that. Core at <b>${f0(s.medCore)} dbar</b>, peak N² ${(s.medPk * 1e5).toFixed(1)}×10⁻⁵ s⁻², against 722 m and 2.3×10⁻⁵ in their 16-year climatology of the western North Atlantic. <b>${s.nCast} of ${casts.length}</b> casts have one; in the other ${casts.length - s.nCast}, both north of 39 °N, N² falls off monotonically — no mode water, so no permanent pycnocline, and they are dropped rather than forced. ` +
        `The edge is the one free choice, and it is the slider: top and base are where N² falls to <b>${(THERMO.frac * 100).toFixed(0)} % of its peak</b>, here a median ${f0(s.medTop)}–${f0(s.medBot)} dbar, <b>${f0(s.med)} dbar thick</b> (at 75 % it is 415 dbar, which is where this walk meets their half-Gaussian fit's 439 m). The two collars on each tube are those depths, measured. ` +
        `<b>Both sheets are surfaces this array can actually map.</b> Hold one cast out, map the other ${s.nCast - 1}, predict it: the base misses by ${f0(s.loo)} dbar against ${f0(s.sd)} for guessing the basin mean — <b>${(s.r2 * 100).toFixed(0)} % of the variance</b> — and the top ${f0(s.tLoo)} against ${f0(s.tSd)}, <b>${(s.tR2 * 100).toFixed(0)} %</b>. So the volume is drawn bowed: down past 1200 dbar in the southwest, up toward 800 in the northeast. That is the contrast with the seasonal layer, which the same test resolves at 12 % — it is measured at 83 points and its shape cannot be mapped from them, and this one can. ` +
        `<span class="n">Each surface gets its own cross-validated covariance (base L = ${f0(s.L)} km, nugget ${s.nug.toFixed(2)}; top L = ${f0(s.tL)} km, nugget ${s.tNug.toFixed(2)}), chosen by leave-one-out over the casts and not by fitting a variogram — a variogram least-squares fit is dominated here by the far separation bins, which outnumber the near ones 30:1, and on the seasonal layer it picked a covariance that predicted a held-out cast worse than the mean did while looking smooth and confident over 93 % of the box. Drag <b>uncertainty weight</b> to 0 and both nuggets go to ${LAY.nugMin}: the same equations become an interpolant through every measured depth, covering the whole box, and nothing about the extra relief cross-validates. Shading is each analysis's own error variance, <b>calibrated</b>: uncalibrated it claimed ${f0(s.predSd)} dbar of error at a held-out cast while actually missing ${f0(s.loo)}, so the error variance is scaled by ${s.cal.toFixed(1)}× on the base and ${s.tCal.toFixed(1)}× on the top to report the error it is measured to make. That is why the volume fades to grey and then to nothing away from the casts instead of sitting solid over the whole domain — the raw analysis error says it knows this field almost everywhere, and cross-validation says it does not. ${s.nTrunc ? `<b>${s.nTrunc}</b> casts stop before N² has fallen back to the edge, so their base is set by where the profile ends, not by the ocean.` : ''}</span>`;
    }
  },
  {
    name: 'What was measured', apply() {
      setVar('conf'); setLayers({ 'tg-casts': true, 'tg-slice': true, 'tg-sect': true, 'tg-iso': false, 'tg-thermo': false });
      $('sl-depth').value = 250; onDepth(); setView('overview');
    }, cap: () => `Colour is now <b>support</b>: how much measured data reaches each point of the interpolated grid (normalised kernel weight — a geometry proxy, not a formal error). ` +
      `Only <b>${(FACTS.nSolid * 100).toFixed(0)} %</b> of the ${NX}×${NY}×${NZ} grid volume has support ≥ 0.5 and only <b>${(FACTS.nAny * 100).toFixed(0)} %</b> has more than 0.08; the rest is the neutral top row of the legend. ` +
      `A smooth full-volume render of this month would be mostly invention. ` +
      `<span class="n">The uncertainty-weight slider is that render: at 1 every layer is drawn as above, at 0 the greying and the fading switch off and the same field is painted as if it were measured everywhere the kernel reaches, at every instant. The kernels are untouched, so it still varies across the box and still moves with the playhead — it just stops admitting what it rests on. Nothing about the data changes between the two.</span>`
  },
  {
    name: 'Time', time: true, apply() {
      setVar('temp'); setLayers({ 'tg-casts': true, 'tg-slice': true, 'tg-sect': false, 'tg-iso': false, 'tg-thermo': false });
      $('sl-depth').value = 25; onDepth();          // 50 dbar, where L_t is short
      setTimeOn(true); setTime(T.tMax);
      setView('overview');
    }, cap: () => `Everything above collapses a month into one picture. It shouldn't, and these casts say so themselves: ` +
      `in the top 150 dbar two of them within 60 km disagree by <b>${f1(FACTS.shNear)} °C</b> when taken within 2 days and <b>${f1(FACTS.shFar)} °C</b> when taken 10+ days apart. ` +
      `Deeper, the same test finds no trend with lag at all (${FACTS.deepBand.p0}–${FACTS.deepBand.p1} dbar: rank correlation ${FACTS.deepBand.rho}, p = ${FACTS.deepBand.p}). So each cast now carries a temporal weight too, ` +
      `<b>L<sub>t</sub> = ${f1(LT[0])} d</b> near the surface and ${f0(LT[NZ - 1])} d in the deep (that deep figure is the record length: no decay was resolvable, which is not the same as none). ` +
      `Drag the playhead or press play. Tubes grey from the top down as the playhead moves away from them, in either direction: the record is complete, so the estimate at an instant is a smoother over the casts either side of it rather than a filter that has only seen the past. ` +
      `<span class="n">Watch the two curves: at 1000 dbar L<sub>t</sub> is longer than the record, so almost every profile counts at every instant and the curve sits near the all-June line; at 50 dbar only the casts within about a week count, so the coverage at any real instant stays a fraction of what the pooled month claims — and no amount of waiting fixes it, because each cast expires before the next arrives.</span>`
  },
  {
    name: 'True scale', trueScale: true, apply() {
      setVar('temp'); setLayers({ 'tg-casts': true, 'tg-slice': false, 'tg-sect': true, 'tg-iso': false, 'tg-thermo': false });
      if (exagBefore === null) toggleTrueScale(); setView('low');
    }, cap: () => `The same volume with no vertical exaggeration: ~${f0(2 * RC)} km across and 2 km deep. Everything above is a sheet of ` +
      `paper. That is why oceanographers stretch the vertical (the badge at the top always says by how much) — and why a few dozen ` +
      `profiles can describe a basin's vertical structure but not its horizontal detail.`
  },
];

// ---------- showcase: opening move, then a five-step introduction ----------
// The panel used to carry three long blocks of prose that nobody reads before
// touching anything. This asks for one action at a time instead, and each step
// clears itself when the viewer does the thing rather than when they click
// "next" -- so the text is a prompt, not a wall.
const AZ0 = Math.atan2(0.95, 1);                   // the scene's default heading

// Nothing is on screen before the step that asks for it. Panels carry .rv and
// gain .in when their moment arrives, and this sequence is the only thing that
// adds them, so a panel cannot appear early by accident.
const REVEAL = {
  descent: ['title', 'exag-badge'],
  '-1': ['guide'],                                 // the disclaimer, before any panel
  0: ['guide', 'tour', 'row-casts', 'legend'],   // the profiles and the ramp
  1: ['compass'],                                  // drag to turn -> a heading
  2: ['insets'],                                   // click a tube -> plots
  3: ['timeline'],                                 // drag the playhead
  4: ['row-thermo'],                               // switch the layer on
  end: ['row-bathy'],                              // the closing card
};
function reveal(key) {
  for (const id of REVEAL[key] || []) $(id).classList.add('in');
  placeCompass();
}
function revealAll() {
  for (const k of Object.keys(REVEAL)) reveal(k);
  $('introcard').classList.add('out');
}
let introHiding = false;                           // tubes held back for the map shot
let ftueStep = -1, ftueT0 = 0;

let navDrag = false, navZoom = false;              // set by the canvas listeners below

const FTUE = [
  {
    text: () => `Each cast is shown by a <b>tube</b>: its temperature profile from two kilometres down to the ` +
      `surface; <b>warm at the top, cold below</b>. ${casts.length} Argo casts over one month in the North Atlantic.`
  },
  {
    text: () => (MOBILE
      ? `<b>Drag with one finger</b> to turn the volume. <b>Pinch</b> to move closer or further out. `
      : `<b>Drag</b> to turn the volume. <b>Scroll</b> to move closer or further out. `) +
      `Depth is stretched <kbd>×200</kbd>.`,
    done: () => navDrag && navZoom, linger: 2500
  },
  {
    text: () => (MOBILE
      ? `<b>Tap a tube</b> to plot it against all ${casts.length} and read what it measured where you touched it. ` +
      `Up to three at once; tap a plotted one again to drop it.`
      : `<b>Click a tube</b> to plot it against all ${casts.length}. Up to three at once. ` +
      `Hover one to read what it measured at that depth.`),
    done: () => pins.length > 0, linger: 2500
  },
  {
    text: () => `<b>Drag the playhead</b> <i class="ph-glyph"></i> along the timeline, left and right. ` +
      `A cast fades as the playhead moves away from it.`,
    done: () => Math.abs(TIME.t - ftueT0) > 0.4, linger: 2500
  },
  {
    text: () => `<b>Switch on Main thermocline.</b> Most of the temperature drop happens inside that one layer.`,
    done: () => $('tg-thermo').checked
  },
];
// Step 00. The register is deliberately flatter than the steps that follow:
// this is the one card that is not selling the instrument. The split it draws
// is data against rendering: the casts are real and cited, and every way they
// are drawn is a proposal for a tool, not a reading of this ocean. Erring
// wide on purpose -- it costs nothing here and a scientist who finds one
// over-claim stops believing the rest.
const DISCLAIMER =
  `<b>This is a concept, not a working tool.</b> The data is real and public, but every ` +
  `representation of it here is conceptual. <b>No OceanX data is used.</b>`;

// what the ring should be pointing at on each step; a cast is projected from
// the scene, everything else is a control that can be found in the DOM
const FOCUS = { 2: { ci: () => focusCast }, 3: { playhead: true }, 4: { sel: '#row-thermo' } };
// the playhead's extent on the timeline canvas, in CSS px from its top: the
// triangle's apex to the foot of the line (drawTimeline's yTick0 - 8 .. yC1 + 3)
const PH_Y0 = -4, PH_Y1 = 77;
let focusCast = 0;
for (let i = 1; i < casts.length; i++) {
  const d = c => Math.hypot(xOf(c.lon), zOf(c.lat));
  if (d(casts[i]) < d(casts[focusCast])) focusCast = i;
}
let focusOn = null;
function setFocus(i) {
  focusOn = FOCUS[i] || null;
  $('focus').classList.toggle('on', !!focusOn);
  $('draghint').classList.toggle('on', !!(focusOn && focusOn.playhead));
}
function updateFocus() {
  if (!focusOn) return;
  const st = $('focus').style;
  if (focusOn.playhead) {
    // the ring sits on the marker itself, not the panel: the marker is the
    // thing to drag, and it moves, so this follows it every frame
    const cv = tl.cv, r = cv.getBoundingClientRect();
    const x = r.left + tlX(TIME.t, cv.clientWidth);
    st.left = (x - 12) + 'px'; st.top = (r.top + PH_Y0 - 6) + 'px';
    st.width = '24px'; st.height = (PH_Y1 - PH_Y0 + 12) + 'px';
    const hs = $('draghint').style;
    hs.left = x + 'px'; hs.top = (r.top + PH_Y0 - 34) + 'px';
    return;
  }
  if (focusOn.sel) {
    const el = document.querySelector(focusOn.sel);
    const r = el && el.getBoundingClientRect();
    if (!r || !r.width) { $('focus').classList.remove('on'); return; }
    st.left = (r.left - 7) + 'px'; st.top = (r.top - 6) + 'px';
    st.width = (r.width + 14) + 'px'; st.height = (r.height + 12) + 'px';
    return;
  }
  const c = casts[focusOn.ci()];
  _v.set(xOf(c.lon), MARK_LIFT, zOf(c.lat)).project(camera);
  if (_v.z > 1) { $('focus').classList.remove('on'); return; }
  $('focus').classList.add('on');
  const S = 48, x = (_v.x + 1) / 2 * innerWidth, y = (1 - _v.y) / 2 * innerHeight;
  st.left = (x - S / 2) + 'px'; st.top = (y - S / 2) + 'px';
  st.width = st.height = S + 'px';
}

// The card is off in the right column while the eye is in the middle of the
// scene, so it has to be able to say "this changed". Two rings shed from its
// edge, three times, on a new step -- and again if the step is waiting on an
// action that has not come. Not while the focus ring is up: that ring is
// already pointing at the control the step is about, and two accent pulses at
// once point at nothing.
const NUDGE_MS = 7800, NUDGE_IDLE = 11000;
let nudgeTimer = 0, ftueShownAt = 0;
function guideNudge() {
  ftueShownAt = performance.now();                 // set first: the idle clock
  const g = $('guide');                            // restarts whether or not
  if (!SHOWCASE || MOBILE || !g) return;            // there is a card to ring
  g.classList.remove('nudge');
  void g.offsetWidth;                              // restart, not queue
  g.classList.add('nudge');
  clearTimeout(nudgeTimer);
  nudgeTimer = setTimeout(() => g.classList.remove('nudge'), NUDGE_MS);
}

// manual: reached by the back button. The step's action has usually been done
// already, so auto-advance is off until the viewer clicks next.
let ftueManual = false, ftueDoneAt = null;
function ftueShow(i, manual = false) {
  ftueStep = i;
  ftueManual = manual;
  ftueDoneAt = null;
  reveal(i);
  setFocus(i);
  $('g-num').textContent = String(Math.min(i + 1, FTUE.length)).padStart(2, '0');
  $('ftue-back').style.visibility = i > 0 ? 'visible' : 'hidden';
  if (i < 0) {
    // step 00: what this is, before anything is claimed by it. It gates the
    // walkthrough rather than sitting beside it, because a viewer who reads
    // one line of a demo reads the first one.
    $('g-num').textContent = '00';
    $('ftue-dots').innerHTML = FTUE.map(() => '<i></i>').join('');
    $('ftue-text').innerHTML = DISCLAIMER;
    $('ftue-skip').style.display = '';
    $('ftue-next').textContent = 'start';
    guideNudge();
    return;
  }
  if (i >= FTUE.length) {
    // the closing card: what is left to do, with no step attached to it
    reveal('end');
    $('g-num').textContent = 'END';
    $('ftue-dots').innerHTML = FTUE.map(() => '<i class="past"></i>').join('');
    $('ftue-text').innerHTML = 'That is the walkthrough. <b>Press play</b> to run the month, ' +
      'switch on <b>Sea floor</b> to see how deep the water under these profiles goes, and keep turning it.';
    $('ftue-back').style.visibility = 'hidden';
    $('ftue-skip').style.display = 'none';
    $('ftue-next').textContent = 'done';
    guideNudge();
    return;
  }
  if (i === 3) ftueT0 = TIME.t;
  $('ftue-text').innerHTML = FTUE[i].text();
  $('ftue-dots').innerHTML = FTUE.map((_, k) =>
    `<i class="${k === i ? 'on' : k < i ? 'past' : ''}"></i>`).join('');
  $('ftue-next').textContent = i === FTUE.length - 1 ? 'done' : 'next';
  guideNudge();
}
// a step clears when its action is done; steps with `linger` wait that long
// after it so the viewer sees the result before the text moves on
function ftueCheck() {
  if (ftueStep < 0 || ftueStep >= FTUE.length) return;
  if (!focusOn && performance.now() - ftueShownAt > NUDGE_IDLE) guideNudge();
  if (ftueManual) return;
  const st = FTUE[ftueStep];
  if (!st.done || !st.done()) { ftueDoneAt = null; return; }
  if (!st.linger) { ftueShow(ftueStep + 1); return; }
  if (ftueDoneAt == null) ftueDoneAt = performance.now();
  else if (performance.now() - ftueDoneAt >= st.linger) ftueShow(ftueStep + 1);
}

// ---------- choreography: one clock, pure functions of it ----------
// The opening used to be a chain of timers, each firing against whatever state
// the last one had left. It is now the model Remotion uses (useCurrentFrame +
// interpolate + spring), ported to plain functions because this scene is
// Three.js under a page that cannot load React: every channel -- a dot's
// scale, the card's opacity, the camera -- is computed from the one clock and
// never accumulated. That is what makes a skip a seek to the end, a still a
// seek to a time (?intro=SECONDS), and a timer unable to fire late.
const bezier = (x1, y1, x2, y2) => {               // CSS cubic-bezier(), Newton on x
  const A = (a, b) => 1 - 3 * b + 3 * a, B = (a, b) => 3 * b - 6 * a, C = a => 3 * a;
  const at = (t, a, b) => ((A(a, b) * t + B(a, b)) * t + C(a)) * t;
  const slope = (t, a, b) => 3 * A(a, b) * t * t + 2 * B(a, b) * t + C(a);
  return x => {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    let t = x;
    for (let i = 0; i < 8; i++) {
      const s = slope(t, x1, x2);
      if (Math.abs(s) < 1e-6) break;
      t -= (at(t, x1, x2) - x) / s;
    }
    return at(t, y1, y2);
  };
};
const EASE = {
  push: bezier(0.16, 1, 0.3, 1),                   // an entrance: quick in, long settle
  crane: bezier(0.6, 0, 0.25, 1),                  // a camera move: off the mark slowly, lands soft
  out: t => 1 - Math.pow(1 - t, 3),
};
// value of x over keyframes xs -> ys, clamped at both ends
function interp(x, xs, ys, ease = t => t) {
  if (x <= xs[0]) return ys[0];
  const n = xs.length - 1;
  if (x >= xs[n]) return ys[n];
  let i = 0;
  while (x > xs[i + 1]) i++;
  return ys[i] + (ys[i + 1] - ys[i]) * ease((x - xs[i]) / (xs[i + 1] - xs[i]));
}
// damped spring, closed form, `sec` seconds after release. Same parameter names
// and defaults as Remotion's spring(): damping 10 / stiffness 100 / mass 1
// overshoots by about a tenth; a heavy damping is its bounce-free "push".
function spring(sec, { damping = 10, stiffness = 100, mass = 1, from = 0, to = 1 } = {}) {
  if (sec <= 0) return from;
  const z = damping / (2 * Math.sqrt(stiffness * mass)), w0 = Math.sqrt(stiffness / mass);
  const x0 = from - to;
  let x;
  if (z < 1) {
    const w1 = w0 * Math.sqrt(1 - z * z);
    x = Math.exp(-z * w0 * sec) * (x0 * Math.cos(w1 * sec) + z * w0 * x0 / w1 * Math.sin(w1 * sec));
  } else if (z === 1) {
    x = Math.exp(-w0 * sec) * x0 * (1 + w0 * sec);
  } else {
    const w2 = w0 * Math.sqrt(z * z - 1);
    x = Math.exp(-z * w0 * sec) * (x0 * Math.cosh(w2 * sec) + z * w0 * x0 / w2 * Math.sinh(w2 * sec));
  }
  return to + x;
}

// ---------- the opening ----------
// A station chart first: the domain circle, empty, seen straight down through
// an orthographic camera -- a map, not a scene. The casts land on it in the
// order they were taken, each one pinged, the first one named, while the
// header counts them in. Then the camera falls: a narrow perspective camera
// matched to the orthographic frame tilts and opens to the working field of
// view while the profiles extrude down from their dots. Nothing on the HUD
// exists before the moment that introduces it.
// Cue sheet, in seconds. Two gates: the bare chart holds until the viewer
// casts the CTDs, the full chart holds until they descend. The clock stops at
// each and the one button carries both actions, so the map is read for as
// long as it takes and the casts land because the viewer asked for them.
const CUES = {
  eyebrow: 0.5, title: 0.8, about: 1.15, rule: 1.5, log: 1.7, skip: 1.7, btn: 1.9,
  gate1: 2.6,                                      // the bare chart; holds for "Cast CTDs"
  first: 2.9,                                      // the first cast lands, and is named
  swarm: 4.1, swarmDur: 3.0,                       // the other 82, chronological
  go: 7.5,                                         // the button returns as "Descend"
  gate: 8.1, fall: 8.1, fallDur: 3.6,              // holds; on click the card leaves, the camera falls
  cardOut: 8.7, descent: 9.5,                     // title and badge, mid-fall
  step: 10.4, end: 11.8,                           // step 01; then the clock stops
};
const POPCFG = { damping: 16, stiffness: 180 };    // ~10 % overshoot, settled in 0.6 s
const PING = { sec: 0.75, grow: 3.4, alpha: 0.42 };
const FOV0 = 5;                                    // near-orthographic start of the fall
const FOV1 = persp.fov;
const H0 = RC * 1.22;                              // half-height of the chart frame
const INTRO_SHIFT = { x: 0.10, y: 0.0 };           // chart sits right of the card, centred vertically
// the working view the fall lands on, and that ?ftue=N cuts to
const WORK = { pres: 820, hor: RC * 2.23, hgt: RC * 0.9, az: AZ0 };
const EL1 = Math.atan2(WORK.hgt, WORK.hor);
const H1 = Math.hypot(WORK.hor, WORK.hgt) * Math.tan(THREE.MathUtils.degToRad(FOV1 / 2));

// landing rings, one per cast: a ping where a float surfaces
const pingGroup = new THREE.Group();
world.add(pingGroup);
{
  const g = new THREE.RingGeometry(7.5, 8.6, 32);
  for (const c of casts) {
    const m = new THREE.Mesh(g, new THREE.MeshBasicMaterial({
      color: 0x1c2128, transparent: true, opacity: 0, side: THREE.DoubleSide, depthWrite: false
    }));
    m.rotation.x = -Math.PI / 2;
    m.position.set(xOf(c.lon), MARK_LIFT / 200, zOf(c.lat));
    m.visible = false;
    pingGroup.add(m);
  }
}
const popOrder = casts.map((c, i) => i).sort((a, b) => casts[a].t - casts[b].t);
const landAt = k => k === 0 ? CUES.first : CUES.swarm + (k - 1) / (popOrder.length - 2) * CUES.swarmDur;
const day = t => dayLabel(t).split(' ').slice(0, 2).join(' ');

let introT0 = null, introFreeze = null, introDone = false, introStage = 0, introTc = 0;
const GATES = [CUES.gate1, CUES.gate];
function introClick() {
  // only at a gate: a click while the casts are still landing is not a descent
  if (introStage < 2 && introTc >= GATES[introStage] - 0.01) introStage++;
}
const fired = new Set();
function cue(name, T, fn) {
  if (T >= CUES[name] && !fired.has(name)) { fired.add(name); fn(); }
}
function introLog() {
  const row = (k, v) => `<dt>${k}</dt><dd>${v}</dd>`;
  $('ic-log').innerHTML =
    row('station set', `<span id="ic-count">0 / ${casts.length} casts</span>`) +
    row('period', `<span id="ic-period">${day(T.tMin)} 2023</span>`) +
    row('bounds', `${G.lat0.toFixed(2)}–${G.lat1.toFixed(2)}°N · ` +
      `${(-G.lon0).toFixed(2)}–${(-G.lon1).toFixed(2)}°W`) +
    row('depth', `0–${G.presMax.toFixed(0)} dbar`) +
    row('source', 'Argo / Ifremer ERDDAP');
  const c = casts[popOrder[0]];
  $('callout').innerHTML = `<span class="co-lead"></span><span class="co-text">` +
    `<span class="co-k">cast 1 of ${casts.length}</span>` +
    `<span class="co-v">float ${c.id} · ${day(c.t)} · ${c.mode === 'D' ? 'delayed-mode' : 'real-time'}</span></span>`;
}
function introStart(now) {
  introLog();
  introT0 = now;
  introHiding = true;                              // no tube to hover under the chart
  onDay();
  setExag(200);
  setDepthAxis(false);
  for (const m of markGroup.children) m.scale.setScalar(0);
  castGroup.scale.y = 1e-4;
  tubeMat.opacity = tubeMatHi.opacity = tubeMatFade.opacity = tubeMatFadeHi.opacity = 0;
  for (const id of ['ic-eyebrow', 'ic-title', 'ic-about', 'ic-rule', 'ic-log', 'ic-skip', 'ic-go'])
    $(id).style.opacity = 0;
  $('ic-rule').style.transform = 'scaleX(0)';
  // straight down, north up, through the orthographic camera: the frustum is
  // sized from the perspective one at this distance, so set that view first
  scShift.x = INTRO_SHIFT.x; scShift.y = INTRO_SHIFT.y; applyViewShift();
  const hgt = H0 / Math.tan(THREE.MathUtils.degToRad(FOV1 / 2));
  applyView({ pres: 0, hor: hgt * 0.005, hgt, az: 0 });
  goOrtho(persp.position.clone().sub(controls.target).normalize());
  controls.enabled = false;                        // on rails until step 01
}
function introTick(now) {
  if (introDone || introT0 == null) return;
  let tc = introFreeze ?? (now - introT0) / 1000;
  // the gates: a still (?intro=) may look past them, a viewer waits at them
  if (introFreeze == null && introStage < 2 && tc > GATES[introStage]) {
    introT0 = now - GATES[introStage] * 1000; tc = GATES[introStage];
  }
  introTc = tc;
  const rm = reducedMotion;
  // --- the card: each line enters on its cue ---
  const enter = (id, at, dy = 10, dur = 0.9) => {
    const u = interp(tc, [at, at + dur], [0, 1], EASE.push);
    const st = $(id).style;
    st.opacity = u;
    st.transform = rm ? 'none' : `translateY(${((1 - u) * dy).toFixed(2)}px)`;
  };
  enter('ic-eyebrow', CUES.eyebrow, 6);
  enter('ic-title', CUES.title, 14, 1.1);
  enter('ic-about', CUES.about, 8);
  // the button: in with the log as "Cast CTDs", out as they land, back as "Descend"
  {
    const u = tc < CUES.go
      ? interp(tc, [CUES.btn, CUES.btn + 0.9, CUES.gate1 + 0.05, CUES.gate1 + 0.35], [0, 1, 1, 0], EASE.push)
      : interp(tc, [CUES.go, CUES.go + 0.9], [0, 1], EASE.push);
    const st = $('ic-go').style;
    st.opacity = u; st.transform = `translateY(${((1 - u) * 6).toFixed(2)}px)`;
    st.pointerEvents = u > 0.9 ? 'auto' : 'none';
  }
  cue('go', tc, () => { $('ic-go').innerHTML = 'Descend &darr;'; });
  enter('ic-log', CUES.log);
  enter('ic-skip', CUES.skip, 4);
  $('ic-rule').style.opacity = interp(tc, [CUES.rule, CUES.rule + 0.3], [0, 1]);
  $('ic-rule').style.transform = `scaleX(${interp(tc, [CUES.rule, CUES.rule + 0.9], [0, 1], EASE.push).toFixed(4)})`;
  // --- the casts land, chronologically; the header counts them in ---
  let shown = 0, last = -1;
  for (let k = 0; k < popOrder.length; k++) {
    const dt = tc - landAt(k);
    const m = markGroup.children[popOrder[k]], p = pingGroup.children[popOrder[k]];
    m.scale.setScalar(dt <= 0 ? 0 : rm ? 1 : spring(dt, POPCFG));
    const q = dt / PING.sec;
    p.visible = !rm && q > 0 && q < 1;
    if (p.visible) { p.scale.setScalar(1 + (PING.grow - 1) * EASE.out(q)); p.material.opacity = PING.alpha * (1 - q); }
    if (dt > 0.04) { shown = k + 1; last = k; }
  }
  if (shown !== introTick.shown) {
    introTick.shown = shown;
    $('ic-count').textContent = `${shown} / ${casts.length} casts`;
    $('ic-period').textContent = shown < 2 ? `${day(T.tMin)} 2023`
      : `${day(T.tMin)} → ${day(casts[popOrder[last]].t)} 2023`;
  }
  // --- the first cast is named: a leader from the dot ---
  {
    const co = $('callout');
    const a = interp(tc, [CUES.first + 0.25, CUES.first + 0.7, CUES.swarm + 1.5, CUES.swarm + 2.0], [0, 1, 1, 0]);
    co.style.opacity = a;
    if (a > 0) {
      const c = casts[popOrder[0]];
      _v.set(xOf(c.lon), MARK_LIFT, zOf(c.lat)).project(camera);
      co.style.left = ((_v.x + 1) / 2 * innerWidth).toFixed(1) + 'px';
      co.style.top = ((1 - _v.y) / 2 * innerHeight).toFixed(1) + 'px';
    }
  }
  // --- the card leaves; the camera falls; the profiles come down with it ---
  const out = interp(tc, [CUES.fall, CUES.fall + 0.55], [0, 1], EASE.out);
  $('introcard').style.opacity = 1 - out;
  $('introcard').style.transform = rm ? 'none' : `translateX(${(-16 * out).toFixed(2)}px)`;
  cue('cardOut', tc, () => $('introcard').classList.add('out'));
  cue('fall', tc, () => {
    introHiding = false;
    onDay();
    // the perspective camera opens narrow and far, on the orthographic frame:
    // the first frames of the fall are the same picture, tilting
    persp.fov = FOV0; persp.updateProjectionMatrix();
    goPersp();
    setDepthAxis(true);
  });
  if (tc >= CUES.fall) {
    const f0 = CUES.fall, D = CUES.fallDur;
    const u = interp(tc, [f0, f0 + D], [0, 1], EASE.crane);
    const uf = interp(tc, [f0 + 0.25 * D, f0 + D], [0, 1], EASE.crane);   // perspective blooms after the tilt begins
    const fov = FOV0 + (FOV1 - FOV0) * uf;
    const H = H0 + (H1 - H0) * u;
    const el = Math.min(Math.PI / 2 - 0.005, Math.PI / 2 + (EL1 - Math.PI / 2) * u);
    const r = H / Math.tan(THREE.MathUtils.degToRad(fov / 2));
    applyView({ pres: WORK.pres * u, hor: r * Math.cos(el), hgt: r * Math.sin(el), az: WORK.az * u, fov });
    scShift.x = INTRO_SHIFT.x + (scBase.x - INTRO_SHIFT.x) * u;
    scShift.y = INTRO_SHIFT.y + (scBase.y - INTRO_SHIFT.y) * u;
    applyViewShift();
    castGroup.scale.y = Math.max(1e-4, interp(tc, [f0 + 0.15 * D, f0 + 0.8 * D], [0, 1], EASE.out));
    // the depth axis is edge-on until the tilt opens it: its labels would pile
    // up in one spot, so they come in once there is a depth to label
    const dl = interp(tc, [f0 + 0.45 * D, f0 + 0.8 * D], [0, 1]);
    for (const L of depthLabels) L.el.style.opacity = dl;
    tubeMat.opacity = tubeMatHi.opacity = tubeMatFade.opacity = tubeMatFadeHi.opacity =
      interp(tc, [f0 + 0.1 * D, f0 + 0.45 * D], [0, 1]);
  }
  cue('descent', tc, () => reveal('descent'));
  cue('step', tc, () => { controls.enabled = true; ftueShow(-1); });
  cue('end', tc, () => {
    introDone = true;
    for (const m of markGroup.children) m.scale.setScalar(1);
    for (const p of pingGroup.children) p.visible = false;
    castGroup.scale.y = 1;
    applyView({ ...WORK, fov: FOV1 });
  });
}
// The skip is a seek: the next frame evaluates the whole sheet at its end and
// lands on the working view with everything in place. A viewer who skips the
// opening is not asking to be walked through it either.
function introSkip() {
  if (introT0 == null || introDone) return;
  introFreeze = null;
  introStage = 2;
  introT0 = performance.now() - CUES.end * 1000 - 1;
  revealAll();
}

// one way in: the rail's checkboxes are the state, everything the showcase
// shows is derived from them, so ?layers= and a click land in the same place
function scSync() {
  if (!SHOWCASE) return;
  $('sc-casts').checked = $('tg-casts').checked;
  $('sc-thermo').checked = $('tg-thermo').checked;
  setGraticule(!$('tg-thermo').checked);
}
if (SHOWCASE) {
  for (const [id, tg] of [['sc-casts', 'tg-casts'], ['sc-thermo', 'tg-thermo']])
    $(id).addEventListener('change', () => { setLayers({ [tg]: $(id).checked }); scSync(); });
  // the sea floor is not one of the rail's layers: it has no field behind it
  // and no support to encode, so it owns its own state
  $('sc-bathy').addEventListener('change', () => {
    setBathy($('sc-bathy').checked);
  });
  $('ftue-next').addEventListener('click', () => {
    if (ftueStep >= FTUE.length) { $('guide').classList.remove('in'); return; }
    ftueShow(ftueStep + 1);
  });
  $('ftue-back').addEventListener('click', () => ftueShow(Math.max(0, ftueStep - 1), true));
  // skip means no more guide: every panel comes out and the card goes, rather
  // than landing on the coda as if the walk had been taken
  $('ftue-skip').addEventListener('click', () => {
    revealAll();
    $('guide').classList.remove('in');
  });
  $('ic-skip').addEventListener('click', introSkip);
  $('ic-go').addEventListener('click', introClick);
} else {
  const wrap = $('tour-steps');
  TOUR.forEach((st, i) => {
    const b = document.createElement('button');
    b.textContent = `${i + 1}. ${st.name}`;
    b.addEventListener('click', () => runTour(i));
    wrap.appendChild(b);
  });
}
const STEPS = TOUR;
function runTour(i) {
  const st = STEPS[i];
  if (!st.trueScale && exagBefore !== null) toggleTrueScale();
  // every other step's caption describes the all-June field, so don't leave
  // time weighting switched on underneath it
  if (!st.time && TIME.on) { TIME.playing = false; syncPlayBtn(); setTimeOn(false); }
  st.apply();
  $('tour-cap').innerHTML = st.cap();
  [...$('tour-steps').children].forEach((b, k) => b.classList.toggle('on', k === i));
}

// ---------- init ----------

setExag(exagFromSlider(+$('sl-exag').value));
// showcase opens deeper-stretched; keep the slider in step so "show true scale"
// restores to this rather than to the console default
if (SHOWCASE) { $('sl-exag').value = Math.round(1000 * Math.log(200) / Math.log(500)); setExag(200); }
$('sl-thermo').value = Math.round(THERMO.frac * 100);
$('v-thermo').textContent = THERMO.frac.toFixed(2);
ISO.on = $('tg-iso').checked; THERMO.on = $('tg-thermo').checked;
onDay(); onDepth(); onSect(); rebuildIso(); buildThermo(); drawLegend(); drawInsets();
computeCurves(); drawTimeline(); timeStatus(); syncPlayBtn();
$('v-lt').textContent = '×' + TIME.scale.toFixed(2);
syncUncUI();
// no hover on a touch screen, so the hint must not ask for one
if (MOBILE) document.querySelector('#insets .hint').textContent =
  'tap a tube to trace and pin it (up to 3). Grey = all 83 profiles.';

// URL params for scripted screenshots: ?var=psal&step=3&mode=expert&exag=1&view=map
{
  const q = new URLSearchParams(location.search);
  if (SHOWCASE) {
    setVar('temp');
    // profiles only: the thermocline is something the introduction asks the
    // viewer to switch on, not something already on when they arrive
    setLayers({ 'tg-casts': true, 'tg-slice': false, 'tg-sect': false, 'tg-iso': false, 'tg-thermo': false });
    // the showcase has no pooled-month state to return to: the playhead is the
    // view. Opens mid-record, so casts fade in both directions from the start
    // rather than only behind it.
    setTimeOn(true); setTime((T.tMin + T.tMax) / 2);
    if (q.get('bathy') === '1') { $('sc-bathy').checked = true; setBathy(true); }
    // ?ftue=skip lands in the working view with no box; ?ftue=note on the
    // disclaimer that gates the walk; ?ftue=N lands there on step N. All exist
    // so a still can be captured of any state without waiting out the opening.
    const fq = q.get('ftue');
    // ?intro=SECONDS freezes the opening at that instant, for stills of any
    // frame of it; ?ftue=card is the frame with every cast in and the header full
    if (fq === 'card' || q.has('intro')) {
      introFreeze = fq === 'card' ? CUES.gate : +q.get('intro');
      introStart(performance.now());
    } else if (fq != null) {
      introDone = true;
      applyView(WORK);
      revealAll();
      ftueShow(fq === 'skip' ? FTUE.length : fq === 'note' ? -1 : +fq);
    } else if (MOBILE) {
      // the opening is a camera move with constants tuned to a wide viewport,
      // and it costs a phone several seconds of the worst frames it will draw
      // before anything responds. Land in the working view and start the
      // walkthrough, which is the part that carries the argument anyway.
      introDone = true;
      applyView(WORK);
      revealAll();
      ftueShow(-1);
    } else introStart(performance.now());
  } else if (q.has('step')) runTour(+q.get('step')); else runTour(0);
  if (q.has('var')) setVar(q.get('var'));
  if (q.has('mode')) setMode(q.get('mode'));
  if (q.has('exag')) setExag(Math.max(1, +q.get('exag')));
  if (q.has('view')) setView(q.get('view'));
  if (q.has('map') && controls.enabled) compassToggle();   // stills of the map view; not while the opening is on rails
  if (q.has('snap')) { const m = { iso: [1, 1, 1], top: [0, 1, 0], e: [1, 0, 0], n: [0, 0, -1] }[q.get('snap')]; if (m) snapView(new THREE.Vector3(...m)); }
  if (q.has('unc')) setUnc(+q.get('unc'));
  if (q.has('ltscale')) {
    TIME.scale = Math.max(0.25, Math.min(4, +q.get('ltscale')));
    $('sl-lt').value = Math.round((Math.log(TIME.scale) / Math.log(4) + 1) / 2 * 1000);
    $('v-lt').textContent = '×' + TIME.scale.toFixed(2);
    syncUncUI();
    computeCurves();
  }
  // ?t=<days since Jun 1> turns the timeline on at that instant; ?t=all is static
  if (q.has('t')) {
    if (q.get('t') === 'all') setTimeOn(false);
    else { TIME.on = true; document.body.classList.add('timed'); $('btn-alltime').textContent = 'show all June'; setTime(+q.get('t')); }
  }
  // ?iso=27.0 sets the surface value on the current variable; ?thermo=1&grad=0.012
  if (q.has('iso')) { setLayers({ 'tg-iso': true }); setIso(+q.get('iso')); }
  if (q.has('frac')) { $('sl-thermo').value = Math.round(+q.get('frac') * 100); $('sl-thermo').dispatchEvent(new Event('input')); }
  // ?kind=seasonal draws the steep upper thermocline instead of the permanent
  // pycnocline; ?render=sheet|cloud changes what is drawn of it (see buildSheet
  // and buildCloud). Both must be set before the layer is switched on.
  if (q.has('kind')) THERMO.kind = q.get('kind');
  if (q.has('render')) THERMO.render = q.get('render');
  if (q.get('mask') === 'off') THERMO.mask = false;
  if (q.get('walls') === 'off') THERMO.walls = false;
  if (q.has('psize')) cloudMat.uniforms.size.value = +q.get('psize');
  if (q.has('colorby')) THERMO.colorBy = q.get('colorby');
  // ?apeak scales the layer's peak opacity. It rescales the whole alpha ramp,
  // so how transparency encodes uncertainty is unchanged — it is a print
  // exposure control, for a layer whose skill is low enough that the on-screen
  // setting renders it almost invisible.
  if (q.has('apeak')) LAY.aPeak = Math.max(0.05, Math.min(1, +q.get('apeak')));
  if (q.has('cloudres')) CLOUD.res = Math.max(0.5, Math.min(3, +q.get('cloudres')));
  // ?grid=off removes the reference frame — graticule, rim, depth axis, scale
  // bar — and every projected label, leaving only the data objects
  if (q.get('grid') === 'off') { frameGroup.visible = false; document.body.classList.add('nogrid'); }
  // ?bg=ffffff overrides the ground colour, for stills that must sit on a
  // white page. The default paper-grey exists so the ramp's bright end and the
  // VSUP neutral stay visible; on white the faintest support fades out sooner,
  // which a print figure can accept.
  if (q.has('bg')) {
    const c = '#' + q.get('bg').replace(/[^0-9a-fA-F]/g, '');
    renderer.setClearColor(c, 1);
    document.body.style.background = c;
  }
  // ?dens scales the cloud's keep-probability; ?synth=1 swaps the layer maps
  // for the invented concept surface (see synthLayer — figures must say so)
  if (q.has('dens')) CLOUD.dens = Math.max(0.1, Math.min(6, +q.get('dens')));
  if (q.get('synth') === '1') SYNTH = true;
  if (q.get('thermo') === '1') setLayers({ 'tg-thermo': true });
  if (q.has('pins')) { pins.push(...q.get('pins').split(',').map(Number).slice(0, 3)); syncPins(); }
  // ?play=1 starts the playhead at t0 and runs it, for capturing clips
  if (q.get('play') === '1') {
    if (!TIME.on) setTimeOn(true);
    setTime(q.has('t') ? TIME.t : T.tMin, true);
    TIME.playing = true; syncPlayBtn();
  }
  // ?layers=casts,thermo turns exactly these on and everything else off; it runs
  // last, so it overrides the layers ?iso= and ?thermo= switch on.
  if (q.has('layers')) {
    const want = new Set(q.get('layers').split(','));
    setLayers(Object.fromEntries(['casts', 'slice', 'sect', 'iso', 'thermo'].map(k => ['tg-' + k, want.has(k)])));
  }
  // the tour step and the layer toggles have both already run by here, and
  // neither rebuilds when only the layer's kind or its rendering changed
  if (['kind', 'render', 'mask', 'colorby', 'cloudres', 'apeak', 'dens', 'synth', 'walls'].some(k => q.has(k))) {
    try { buildThermo(); } catch (e) { console.log('buildThermo failed: ' + e.stack); }
  }
  // ?cam=az,el,r[,dbar] — azimuth ° clockwise from north, elevation ° above the
  // horizon, radius in domain radii, and optionally the depth to look at, which
  // is what centres a still: the named views aim at a fixed height, so raising
  // the exaggeration slides the object out of frame.
  if (q.has('cam')) {
    const [az, el, r, look] = q.get('cam').split(',').map(Number);
    const a = az * Math.PI / 180, e = el * Math.PI / 180, R = (r || 2) * RC;
    if (camera === ortho) goPersp();
    if (Number.isFinite(look)) controls.target.y = yOf(look) * exag;
    camera.position.set(
      controls.target.x + R * Math.cos(e) * Math.sin(a),
      controls.target.y + R * Math.sin(e),
      controls.target.z + R * Math.cos(e) * Math.cos(a));
    controls.update();
  }
  // ?ui=min hides the console but keeps the legend and badges; ?ui=none hides
  // those too. Stills for the proposal figures; see .figure in style.css.
  scSync();
  if (q.has('ui')) {
    document.body.classList.add('figure');
    if (q.get('ui') === 'none') document.body.classList.add('nolegend');
  }
}

// ---------- small screen: the notice and the sheet ----------
// The panels that have nowhere to sit on a phone's bottom edge are moved into
// #sheet and slide up over the scene on demand. The nodes are moved, not
// copied: every listener, every canvas and every id lookup already wired up
// stays valid, and the move reverses if the viewport grows back past the
// breakpoint. The showcase keeps its walkthrough card and its toggles in the
// grid — those are its content, not its controls — so only the plots and the
// colour key go in.
const SHEET_PANELS = SHOWCASE
  ? ['insets', 'legend']
  : ['tour', 'rail', 'insets', 'legend', 'prov-box'];
const SHEET_LABEL = SHOWCASE ? 'Plots & key' : 'Controls';
const sheetHome = new Map();                       // id -> where it came from
let sheetMoved = false;
function setSheet(open) {
  document.body.classList.toggle('sheet-open', open);
  $('sheet-btn').setAttribute('aria-expanded', String(open));
  $('sheet-btn').textContent = open ? 'Close' : SHEET_LABEL;
  if (open) { sizeInsets(); drawInsets(); drawLegend(); }   // sized from a box that was 0 wide while closed
}
function syncMobileDom() {
  if (MOBILE === sheetMoved) return;
  const sheet = $('sheet');
  if (MOBILE) {
    for (const id of SHEET_PANELS) {
      const el = $(id);
      sheetHome.set(id, [el.parentNode, el.nextSibling]);
      sheet.appendChild(el);
    }
  } else {
    for (const id of SHEET_PANELS) {
      const [parent, next] = sheetHome.get(id);
      parent.insertBefore($(id), next);
    }
    setSheet(false);
  }
  sheetMoved = MOBILE;
  setSheet(document.body.classList.contains('sheet-open') && MOBILE);
}
$('sheet-btn').addEventListener('click', () =>
  setSheet(!document.body.classList.contains('sheet-open')));

// The notice is dismissed for the tab, not for good: a reload after a rotation
// should not nag again, but a link opened tomorrow should still say it. Its
// height feeds --note-h, which is the HUD's extra top padding, so it pushes the
// header down rather than covering it.
const NOTE_KEY = 'wc-small-screen-note';
const NOTE_FLOAT_MQ = matchMedia('(max-height: 420px)');
function syncNote() {
  let dismissed = false;
  try { dismissed = sessionStorage.getItem(NOTE_KEY) === '1'; } catch { /* private mode */ }
  const on = MOBILE && !dismissed && !document.body.classList.contains('figure');
  document.body.classList.toggle('mn-off', !on);
  // on a viewport this short there is no height to give it: it stops reserving
  // a row and overlays the header instead, which it can afford to do because
  // it is the one thing on the page whose whole job is to be dismissed
  const reserve = on && !NOTE_FLOAT_MQ.matches;
  document.body.classList.toggle('mn-float', on && !reserve);
  document.documentElement.style.setProperty(
    '--note-h', (reserve ? $('mobilenote').offsetHeight : 0) + 'px');
}
$('mn-close').addEventListener('click', () => {
  try { sessionStorage.setItem(NOTE_KEY, '1'); } catch { /* private mode */ }
  syncNote();
});

// A view tuned on a wide screen frames the domain by its *vertical* field, so a
// portrait phone gets the same height and a third of the width. Pull the camera
// back by the square root of the shortfall — the whole shortfall would put the
// column too far off to read — once, at load. After that the framing is the
// viewer's and a rotation has no business undoing their zoom.
function frameForViewport() {
  const k = Math.min(1.9, Math.sqrt((16 / 9) / Math.max(0.2, innerWidth / innerHeight)));
  if (k <= 1.001) return;
  camera.position.sub(controls.target).multiplyScalar(k).add(controls.target);
  controls.update();
}

// ---------- layout ----------
// One pass for everything that depends on the viewport: the cameras and the
// renderer; the method note's fold; every 2D canvas, whose bitmap is its CSS
// box, so a resize re-allocates and redraws it rather than letting CSS stretch
// a stale one (the timeline's playhead used to be drawn at one width and
// dragged at another); and the showcase view shift. Coalesced to one run per
// frame, and run once here after the URL parameters have set the mode.
// Where the showcase's colour ramp goes is a measurement, like placeCompass().
// The left column is title / free row / toggles / ramp / plots, and every one
// of those is sized in rem while the root font-size follows the viewport
// *width* -- the stack wants ~660 px of height at 1400 across and ~870 at
// 2700. A single height breakpoint (859) was therefore right only at the wide
// end and pushed the ramp into the bottom right of every ordinary laptop
// window, which has room for it. Sum the intrinsic heights instead and move it
// only when the free row would be squeezed under LEGEND_SLACK.
//
// It is not monotonic in height, and that is the plots' 699 px breakpoint, not
// this: over a ~60 px band the column is short of room, and just under 699 the
// plots shrink and it has room again.
//
// The sum does not depend on the answer: nothing in that column stretches
// (#title is align-self:start, the toggles, ramp and plots are align-self:end),
// so each panel measures the same in either layout and the class cannot
// oscillate. Below 1200 px across the ramp is hidden outright and the class
// must come off, or its template would reserve a column for a panel that is
// not drawn.
const LEGEND_SLACK = 48;
function placeLegend() {
  if (!SHOWCASE) return;
  const legend = $('legend');
  if (MOBILE || getComputedStyle(legend).display === 'none') {
    document.body.classList.remove('legend-bottom');
    return;
  }
  const hs = getComputedStyle($('hud'));
  const gap = parseFloat(hs.rowGap) || 0;
  const need = parseFloat(hs.paddingTop) + parseFloat(hs.paddingBottom) + 4 * gap
    + $('title').offsetHeight + $('tour').offsetHeight
    + legend.offsetHeight + $('insets').offsetHeight;
  document.body.classList.toggle('legend-bottom', innerHeight - need < LEGEND_SLACK);
}

const PROV_OPEN_MQ = matchMedia('(min-width: 1440px) and (min-height: 860px)');
let provOpenWas = null, layoutQueued = false;
function relayout() {
  layoutQueued = false;
  // "small" is re-read here rather than latched at load, so a desktop window
  // dragged across the breakpoint, or a phone rotated, lands on the right
  // layout and the right input model
  if (MOBILE_FORCE == null) MOBILE = MOBILE_MQ.matches;
  document.body.classList.toggle('mobile', MOBILE);
  syncMobileDom();
  syncNote();
  persp.aspect = innerWidth / innerHeight;
  persp.updateProjectionMatrix();
  if (camera === ortho) { const h = ortho.top; ortho.left = -h * persp.aspect; ortho.right = h * persp.aspect; ortho.updateProjectionMatrix(); }
  renderer.setSize(innerWidth, innerHeight);
  // the note is open on a big screen and folded on a small one; a manual
  // toggle stands until the viewport crosses the line again
  if (PROV_OPEN_MQ.matches !== provOpenWas) { provOpenWas = PROV_OPEN_MQ.matches; $('prov-box').open = provOpenWas; }
  sizeInsets(); placeLegend();
  drawInsets(); drawLegend(); drawTimeline();
  placeCompass();
  computeViewShift();
  if (introT0 == null || introDone) { scShift.x = scBase.x; scShift.y = scBase.y; }
  applyViewShift();
  updateFocus();
}
addEventListener('resize', () => {
  if (layoutQueued) return;
  layoutQueued = true;
  requestAnimationFrame(relayout);
});
relayout();
if (MOBILE_AT_LOAD) frameForViewport();
// ?sheet=1 opens the mobile sheet, so a capture can land on it like any other state
if (MOBILE && new URLSearchParams(location.search).get('sheet') === '1') setSheet(true);

let lastFrame = 0;
function frame(now) {
  requestAnimationFrame(frame);
  const dts = lastFrame ? Math.min(0.1, (now - lastFrame) / 1000) : 0;
  lastFrame = now;
  if (TIME.playing) {
    // the 83 profile traces in the 2D insets are the expensive part of a
    // rebuild, so they are held back to ~6 Hz while the playhead runs
    const t = TIME.t + dts * TIME.speed;
    const quick = (now - (frame.insAt || 0)) < 160;
    if (!quick) frame.insAt = now;
    if (t >= T.tMax) { TIME.playing = false; setTime(T.tMax); syncPlayBtn(); }   // stopped first: the last tick is a full one
    else setTime(t, quick);
  }
  if (tween) {
    const k = Math.min(1, (now - tween.t0) / tween.ms);
    const ease = 1 - Math.pow(1 - k, 3);
    setExag(Math.exp(Math.log(tween.from) + (Math.log(tween.to) - Math.log(tween.from)) * ease));
    if (k >= 1) tween = null;
  }
  introTick(now);
  if (bathyFade) {
    const k = Math.min(1, (now - bathyFade.t0) / bathyFade.ms);
    bathyMat.opacity = bathyFade.from + (bathyFade.to - bathyFade.from) * (1 - Math.pow(1 - k, 3));
    if (k >= 1) { bathyGroup.visible = bathyFade.to > 0; bathyFade = null; }
  }
  if (viewTween) {
    const k = Math.min(1, (now - viewTween.t0) / viewTween.ms);
    const e = 1 - Math.pow(1 - k, 3);
    const f = viewTween.from, t = viewTween.to;
    let daz = t.az - f.az;
    while (daz > Math.PI) daz -= 2 * Math.PI;
    while (daz < -Math.PI) daz += 2 * Math.PI;
    applyView({
      pres: f.pres + (t.pres - f.pres) * e,
      hor: f.hor + (t.hor - f.hor) * e,
      hgt: f.hgt + (t.hgt - f.hgt) * e,
      az: f.az + daz * e,
    });
    if (k >= 1) viewTween = null;
  }
  if (camTween) {
    const k = Math.min(1, (now - camTween.t0) / camTween.ms);
    const e = 1 - Math.pow(1 - k, 3);
    _d.copy(camTween.from).lerp(camTween.to, e).normalize();   // nlerp is fine for < 180°
    camera.position.copy(controls.target)
      .addScaledVector(_d, camTween.r + (camTween.r1 - camTween.r) * e);
    if (k >= 1) { const t = camTween; camTween = null; controls.update(); t.then?.(); }
  }
  ftueCheck();
  updateFocus();
  controls.update();
  renderer.render(scene, camera);
  updateLabels();
  drawGizmo();
  drawCompass();
  placeReadout();
}
requestAnimationFrame(frame);

// ?bench=1: per-stage timings of a timeline tick and of a frame, written to
// document.title so a headless capture can read them off --dump-dom. Throwaway
// measurement, not a feature.
if (new URLSearchParams(location.search).has('bench')) {
  const out = [];
  const say = s => { out.push(s); console.log('BENCH ' + s); document.title = 'BENCH|' + out.join('|'); };
  const med = a => { const s = [...a].sort((x, y) => x - y); return s[s.length >> 1]; };
  setTimeout(() => {
    setLayers({ 'tg-casts': true, 'tg-slice': true, 'tg-sect': true, 'tg-iso': true, 'tg-thermo': true });
    if (!TIME.on) setTimeOn(true);
    const tOf = i => T.tMin + (T.tMax - T.tMin) * ((0.13 + i * 0.37) % 1);
    const stages = [
      ['field', () => rebuildTimeField()], ['recolor', () => recolorCasts()],
      ['slice', () => updateSlice(slice.pres)], ['sect', () => onSect()],
      ['iso', () => buildIsoSheets()], ['thermo', () => buildThermo()],
      ['timeline', () => drawTimeline()], ['status', () => timeStatus()], ['insets', () => drawInsets()],
    ];
    const sub = {};
    for (const nm of ['solveLayer', 'layerField', 'heightSurface', 'fishnet', 'median3x3', 'expandFine', 'paintField', 'drawContours', 'drawProfilePanel', 'drawTSPanel', 'oiSetup', 'layerCasts', 'buildSheet']) {
      const orig = eval(nm);
      const w = function (...a) { const t0 = performance.now(); const r = orig.apply(this, a); sub[nm] = (sub[nm] || 0) + performance.now() - t0; return r; };
      eval(nm + ' = w');
    }
    {
      const proto = THREE.BufferGeometry.prototype, orig = proto.computeVertexNormals;
      proto.computeVertexNormals = function () { const t0 = performance.now(); orig.call(this); sub.normals = (sub.normals || 0) + performance.now() - t0; };
    }
    const acc = stages.map(() => []);
    for (let i = 0; i < 10; i++) {
      TIME.t = tOf(i);
      stages.forEach(([, fn], k) => { const t0 = performance.now(); fn(); acc[k].push(performance.now() - t0); });
    }
    say('tick ' + stages.map(([n], k) => `${n}=${med(acc[k]).toFixed(1)}`).join(' '));
    say('sub ' + Object.entries(sub).map(([k, v]) => `${k}=${(v / 10).toFixed(1)}`).join(' '));
    let nObj = 0, nTri = 0; scene.traverse(o => { if (o.isMesh || o.isLine || o.isLineSegments || o.isPoints) nObj++; });
    say(`objects=${nObj} groups cast=${castGroup.children.length} mark=${markGroup.children.length} isoRing=${isoRingGroup.children.length} thermo=${thermoGroup.children.length} iso=${isoGroup.children.length}`);
    // frame: wrap the per-frame calls and orbit for 90 frames
    const tm = { render: [], labels: [], gizmo: [], frame: [] };
    const wrap = (fn, key) => function (...a) { const t0 = performance.now(); const r = fn.apply(this, a); tm[key].push(performance.now() - t0); return r; };
    renderer.render = wrap(renderer.render.bind(renderer), 'render');
    updateLabels = wrap(updateLabels, 'labels');
    drawGizmo = wrap(drawGizmo, 'gizmo');
    controls.autoRotate = true; controls.autoRotateSpeed = 30;
    TIME.playing = false;
    let n = 0, last = 0;
    const tickF = now => {
      if (last) tm.frame.push(now - last); last = now;
      if (++n < 90) { requestAnimationFrame(tickF); return; }
      const inf = renderer.info.render;
      say(`nav render=${med(tm.render).toFixed(1)} labels=${med(tm.labels).toFixed(2)} gizmo=${med(tm.gizmo).toFixed(2)} frame=${med(tm.frame).toFixed(1)} calls=${inf.calls} tris=${inf.triangles} lines=${inf.lines} px=${renderer.domElement.width}x${renderer.domElement.height} dpr=${renderer.getPixelRatio()} msaa=${renderer.getContext().getParameter(renderer.getContext().SAMPLES)}`);
      // play: same, with the timeline running
      for (const k in tm) tm[k].length = 0;
      TIME.playing = true; setTime(T.tMin, true); n = 0; last = 0;
      const tickP = now => {
        if (last) tm.frame.push(now - last); last = now;
        if (++n < 90) { requestAnimationFrame(tickP); return; }
        say(`play render=${med(tm.render).toFixed(1)} frame=${med(tm.frame).toFixed(1)} calls=${renderer.info.render.calls}`);
        say('done');
      };
      requestAnimationFrame(tickP);
    };
    requestAnimationFrame(tickF);
  }, 1500);
}

// ?layoutcheck=1: HUD panels that overlap each other or leave the viewport,
// written to document.title for a headless capture to read off --dump-dom.
// Throwaway measurement, like ?bench=1.
if (new URLSearchParams(location.search).has('layoutcheck')) {
  setTimeout(() => {
    const ids = ['title', 'gizmo', 'compass', 'tour', 'guide', 'rail', 'insets', 'timeline', 'legend', 'prov-box', 'exag-badge'];
    const R = [];
    for (const id of ids) {
      const el = $(id);
      if (!el || getComputedStyle(el).display === 'none') continue;
      // the mobile sheet is a scroll container that overlays the bottom band on
      // purpose: what is inside it is neither out of place nor out of bounds
      if (el.closest('#sheet')) continue;
      const r = el.getBoundingClientRect();
      if (r.width && r.height) R.push({ id, r });
    }
    const bad = [];
    for (let i = 0; i < R.length; i++) {
      const a = R[i].r;
      if (a.left < -1 || a.top < -1 || a.right > innerWidth + 1 || a.bottom > innerHeight + 1) bad.push(R[i].id + '>viewport');
      for (let j = i + 1; j < R.length; j++) {
        const b = R[j].r;
        const ox = Math.min(a.right, b.right) - Math.max(a.left, b.left);
        const oy = Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top);
        if (ox > 4 && oy > 4) bad.push(`${R[i].id}~${R[j].id}`);
      }
    }
    document.title = `LAYOUT|${innerWidth}x${innerHeight}|` + (bad.length ? bad.join(',') : 'ok');
  }, 1500);
}





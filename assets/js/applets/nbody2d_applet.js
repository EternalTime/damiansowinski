(function () {
'use strict';

const _cs  = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();
const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

// ── Inject CSS for power spectrum panel ───────────────────────────────────────
(function () {
  if (document.getElementById('nb-styles')) return;
  const s = document.createElement('style');
  s.id = 'nb-styles';
  s.textContent = `
    #nb-ctrl-panel { display:flex; flex-direction:column; overflow:hidden; }
    .nb-pspec-section { flex:1; min-height:0; display:flex; flex-direction:column; padding:8px 12px 10px; }
    #nb-pspec { flex:1; min-height:0; width:100%; display:block; }
    #nb-xi    { flex:1; min-height:0; width:100%; display:block; }
  `;
  document.head.appendChild(s);
})();

// ── Simulation parameters ─────────────────────────────────────────────────────
const N       = 10000;
let G         = 1.0 / N;
let softening = 0.004;
let theta     = 0.3;
let dt        = 0.002;
const STEPS_PER_FRAME = 1;
let expanding = false;
let lambda    = 0.0;
let spinFrac  = 0.2;
let zoom      = 1.0;

// ── Power spectrum & correlation function ─────────────────────────────────────
const GRID     = 256;   // NGP grid size
const N_KBINS  = 256;   // radial k bins
const N_RBINS  = 256;   // radial r bins for xi(r)
let pspecCanvas = null, pspecCtx = null;
let xiCanvas    = null, xiCtx    = null;
let frameCount  = 0;

// Radix-2 FFT (in-place, re/im interleaved in two separate arrays)
function fft(re, im, n) {
  // bit-reversal permutation
  let j = 0;
  for (let i = 1; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      let t = re[i]; re[i] = re[j]; re[j] = t;
          t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = -2 * Math.PI / len;
    const wRe = Math.cos(ang), wIm = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let curRe = 1, curIm = 0;
      for (let k = 0; k < len >> 1; k++) {
        const uRe = re[i+k], uIm = im[i+k];
        const vRe = re[i+k+len/2]*curRe - im[i+k+len/2]*curIm;
        const vIm = re[i+k+len/2]*curIm + im[i+k+len/2]*curRe;
        re[i+k]        = uRe + vRe; im[i+k]        = uIm + vIm;
        re[i+k+len/2]  = uRe - vRe; im[i+k+len/2]  = uIm - vIm;
        const newRe = curRe*wRe - curIm*wIm;
        curIm = curRe*wIm + curIm*wRe; curRe = newRe;
      }
    }
  }
}

function computePowerSpectrum() {
  const M = GRID;
  const re = new Float32Array(M * M);
  const im = new Float32Array(M * M);

  // NGP deposit — use actual bounding box
  let x0 = px[0], y0 = py[0], x1 = px[0], y1 = py[0];
  for (let i = 1; i < N; i++) {
    if (px[i] < x0) x0 = px[i]; if (px[i] > x1) x1 = px[i];
    if (py[i] < y0) y0 = py[i]; if (py[i] > y1) y1 = py[i];
  }
  const rx = x1 - x0 || 1, ry = y1 - y0 || 1;
  for (let i = 0; i < N; i++) {
    const gi = Math.min(Math.floor((px[i] - x0) / rx * M), M-1);
    const gj = Math.min(Math.floor((py[i] - y0) / ry * M), M-1);
    re[gj * M + gi]++;
  }

  // Subtract mean (sets k=0 mode to zero)
  const mean = N / (M * M);
  for (let k = 0; k < M * M; k++) re[k] -= mean;

  // 2D FFT via row then column 1D FFTs
  const row = new Float32Array(M), rim = new Float32Array(M);
  for (let j = 0; j < M; j++) {
    for (let i = 0; i < M; i++) { row[i] = re[j*M+i]; rim[i] = 0; }
    fft(row, rim, M);
    for (let i = 0; i < M; i++) { re[j*M+i] = row[i]; im[j*M+i] = rim[i]; }
  }
  const col = new Float32Array(M), cim = new Float32Array(M);
  for (let i = 0; i < M; i++) {
    for (let j = 0; j < M; j++) { col[j] = re[j*M+i]; cim[j] = im[j*M+i]; }
    fft(col, cim, M);
    for (let j = 0; j < M; j++) { re[j*M+i] = col[j]; im[j*M+i] = cim[j]; }
  }

  // Radial binning of P(k)
  const pk  = new Float32Array(N_KBINS);
  const pkN = new Int32Array(N_KBINS);
  const kMax = M / 2;
  for (let j = 0; j < M; j++) {
    const ky = j < M/2 ? j : j - M;
    for (let i = 0; i < M; i++) {
      const kx = i < M/2 ? i : i - M;
      const k  = Math.sqrt(kx*kx + ky*ky);
      const bin = Math.floor(k / kMax * N_KBINS);
      if (bin >= 0 && bin < N_KBINS) {
        pk[bin]  += re[j*M+i]*re[j*M+i] + im[j*M+i]*im[j*M+i];
        pkN[bin]++;
      }
    }
  }
  for (let b = 0; b < N_KBINS; b++) if (pkN[b] > 0) pk[b] /= pkN[b];

  // Inverse FFT of P(k) field to get xi(r): set each mode to its power, im=0, then IFFT
  const xire = new Float32Array(M * M);
  const xiim = new Float32Array(M * M);
  for (let k = 0; k < M * M; k++) {
    xire[k] = re[k]*re[k] + im[k]*im[k];
    xiim[k] = 0;
  }
  // IFFT = conjugate, FFT, conjugate, divide by M²
  for (let k = 0; k < M * M; k++) xiim[k] = -xiim[k]; // already 0
  const xrow = new Float32Array(M), xrim = new Float32Array(M);
  for (let j = 0; j < M; j++) {
    for (let i = 0; i < M; i++) { xrow[i] = xire[j*M+i]; xrim[i] = xiim[j*M+i]; }
    fft(xrow, xrim, M);
    for (let i = 0; i < M; i++) { xire[j*M+i] = xrow[i] / (M*M); xiim[j*M+i] = -xrim[i] / (M*M); }
  }
  const xcol = new Float32Array(M), xcim = new Float32Array(M);
  for (let i = 0; i < M; i++) {
    for (let j = 0; j < M; j++) { xcol[j] = xire[j*M+i]; xcim[j] = xiim[j*M+i]; }
    fft(xcol, xcim, M);
    for (let j = 0; j < M; j++) { xire[j*M+i] = xcol[j] / (M*M); }
  }

  // Radial binning of xi(r) — only out to M/4
  const xi   = new Float32Array(N_RBINS);
  const xiN  = new Int32Array(N_RBINS);
  const rMax = M / 4;
  for (let j = 0; j < M; j++) {
    const dy = j < M/2 ? j : j - M;
    for (let i = 0; i < M; i++) {
      const dx = i < M/2 ? i : i - M;
      const r  = Math.sqrt(dx*dx + dy*dy);
      const bin = Math.floor(r / rMax * N_RBINS);
      if (bin >= 0 && bin < N_RBINS) { xi[bin] += xire[j*M+i]; xiN[bin]++; }
    }
  }
  for (let b = 0; b < N_RBINS; b++) if (xiN[b] > 0) xi[b] /= xiN[b];

  drawTwoPoint(xi);
  drawPowerSpectrum(pk);
}

function drawTwoPoint(xi) {
  if (!xiCanvas) return;
  const W = xiCanvas.width, H = xiCanvas.height;
  if (W === 0 || H === 0) return;
  const PAD_L = 8, PAD_R = 8, PAD_T = 8, PAD_B = 8;
  const pw = W - PAD_L - PAD_R;
  const ph = H - PAD_T - PAD_B;

  xiCtx.clearRect(0, 0, W, H);

  // semilog-x: linear xi on y, log10(r) on x
  let maxXi = -Infinity, minXi = Infinity;
  for (let b = 1; b < N_RBINS; b++) {
    if (xi[b] > maxXi) maxXi = xi[b];
    if (xi[b] < minXi) minXi = xi[b];
  }
  const rangeXi = maxXi - minXi || 1;

  const logRMin = Math.log10(0.5);
  const logRMax = Math.log10(N_RBINS - 0.5);
  const rangeR  = logRMax - logRMin || 1;

  const xOf = b => PAD_L + (Math.log10(b + 0.5) - logRMin) / rangeR * pw;
  const yOf = v => PAD_T + ph - (v - minXi) / rangeXi * ph;

  const y0line = yOf(0);

  // filled curve
  xiCtx.beginPath();
  xiCtx.moveTo(xOf(1), yOf(xi[1]));
  for (let b = 2; b < N_RBINS; b++) xiCtx.lineTo(xOf(b), yOf(xi[b]));
  xiCtx.lineTo(xOf(N_RBINS - 1), y0line);
  xiCtx.lineTo(xOf(1),           y0line);
  xiCtx.closePath();
  xiCtx.fillStyle = _rgba('--pink-dark', 0.25);
  xiCtx.fill();

  // zero line
  xiCtx.beginPath();
  xiCtx.moveTo(PAD_L, y0line); xiCtx.lineTo(PAD_L + pw, y0line);
  xiCtx.strokeStyle = _rgba('--text-dim', 0.4);
  xiCtx.lineWidth = 1; xiCtx.setLineDash([4, 6]);
  xiCtx.stroke(); xiCtx.setLineDash([]);

  // line
  xiCtx.beginPath();
  xiCtx.moveTo(xOf(1), yOf(xi[1]));
  for (let b = 2; b < N_RBINS; b++) xiCtx.lineTo(xOf(b), yOf(xi[b]));
  xiCtx.strokeStyle = _rgba('--pink-light', 0.9);
  xiCtx.lineWidth   = 1.5;
  xiCtx.stroke();
}

function drawPowerSpectrum(pk) {
  if (!pspecCanvas) return;
  const W = pspecCanvas.width, H = pspecCanvas.height;
  if (W === 0 || H === 0) return;
  const PAD_L = 8, PAD_R = 8, PAD_T = 8, PAD_B = 8;
  const pw = W - PAD_L - PAD_R;
  const ph = H - PAD_T - PAD_B;

  pspecCtx.clearRect(0, 0, W, H);

  // bin centres in k, log10
  const kMax = GRID / 2;
  const dk   = kMax / N_KBINS;
  const logK  = new Float32Array(N_KBINS);
  const logPk = new Float32Array(N_KBINS);
  let maxLogK = -Infinity, minLogK = Infinity;
  let maxLogP = -Infinity;
  for (let b = 0; b < N_KBINS; b++) {
    const kc = (b + 0.5) * dk;
    logK[b]  = Math.log10(kc);
    logPk[b] = pk[b] > 0 ? Math.log10(pk[b]) : NaN;
    if (logK[b] > maxLogK) maxLogK = logK[b];
    if (logK[b] < minLogK) minLogK = logK[b];
    if (!isNaN(logPk[b]) && logPk[b] > maxLogP) maxLogP = logPk[b];
  }
  minLogK = Math.log10(0.5 * dk);
  const minLogP = maxLogP - 6;
  const rangeK = maxLogK - minLogK || 1;
  const rangeP = maxLogP - minLogP || 1;

  const xOf = lk => PAD_L + (lk - minLogK) / rangeK * pw;
  const yOf = lp => PAD_T + ph - (lp - minLogP) / rangeP * ph;

  // filled curve
  pspecCtx.beginPath();
  let started = false;
  for (let b = 0; b < N_KBINS; b++) {
    if (isNaN(logPk[b])) continue;
    const x = xOf(logK[b]), y = yOf(logPk[b]);
    if (!started) { pspecCtx.moveTo(x, y); started = true; }
    else pspecCtx.lineTo(x, y);
  }
  pspecCtx.lineTo(xOf(maxLogK), PAD_T + ph);
  pspecCtx.lineTo(xOf(minLogK), PAD_T + ph);
  pspecCtx.closePath();
  pspecCtx.fillStyle = _rgba('--teal-dark', 0.3);
  pspecCtx.fill();

  // line
  pspecCtx.beginPath();
  started = false;
  for (let b = 0; b < N_KBINS; b++) {
    if (isNaN(logPk[b])) continue;
    const x = xOf(logK[b]), y = yOf(logPk[b]);
    if (!started) { pspecCtx.moveTo(x, y); started = true; }
    else pspecCtx.lineTo(x, y);
  }
  pspecCtx.strokeStyle = _rgba('--teal-light', 0.9);
  pspecCtx.lineWidth   = 1.5;
  pspecCtx.stroke();
}

// ── Simulation arrays ─────────────────────────────────────────────────────────
let px, py, vx, vy, ax, ay;

// ── WebGL state ───────────────────────────────────────────────────────────────
let gl = null, simCanvas = null;
let prog = null;
let posBuf = null, spdBuf = null;
let posArr = null;   // Float32Array [x0,y0, x1,y1, ...]
let spdArr = null;   // Float32Array [speed0, speed1, ...]
let uPointSize, uOffset, uZoom;
let uColSlow, uColFast;
let followCoM = false;

let running = false, frameId = null;

// ── Barnes-Hut quadtree (flat node pool) ─────────────────────────────────────
const MAX_NODES = 1 << 20;
const nd_cx    = new Float64Array(MAX_NODES);
const nd_cy    = new Float64Array(MAX_NODES);
const nd_mass  = new Float64Array(MAX_NODES);
const nd_x0    = new Float64Array(MAX_NODES);
const nd_y0    = new Float64Array(MAX_NODES);
const nd_x1    = new Float64Array(MAX_NODES);
const nd_y1    = new Float64Array(MAX_NODES);
const nd_child = new Int32Array(MAX_NODES * 4);
let nodeCount  = 0;

function allocNode(x0, y0, x1, y1) {
  const i = nodeCount++;
  nd_cx[i] = 0; nd_cy[i] = 0; nd_mass[i] = 0;
  nd_x0[i] = x0; nd_y0[i] = y0; nd_x1[i] = x1; nd_y1[i] = y1;
  nd_child[i*4] = nd_child[i*4+1] = nd_child[i*4+2] = nd_child[i*4+3] = -1;
  return i;
}

function insertParticle(node, ix, iy, im) {
  if (nd_mass[node] === 0) {
    nd_cx[node] = ix; nd_cy[node] = iy; nd_mass[node] = im;
    return;
  }
  const mx = (nd_x0[node] + nd_x1[node]) * 0.5;
  const my = (nd_y0[node] + nd_y1[node]) * 0.5;
  const base = node * 4;
  const hasChildren = nd_child[base] !== -1 || nd_child[base+1] !== -1 ||
                      nd_child[base+2] !== -1 || nd_child[base+3] !== -1;
  if (!hasChildren) {
    const ex = nd_cx[node], ey = nd_cy[node], em = nd_mass[node];
    nd_child[base]   = allocNode(nd_x0[node], nd_y0[node], mx,          my);
    nd_child[base+1] = allocNode(mx,          nd_y0[node], nd_x1[node], my);
    nd_child[base+2] = allocNode(nd_x0[node], my,          mx,          nd_y1[node]);
    nd_child[base+3] = allocNode(mx,          my,          nd_x1[node], nd_y1[node]);
    const eq = (ex < mx ? 0 : 1) + (ey < my ? 0 : 2);
    insertParticle(nd_child[base + eq], ex, ey, em);
  }
  const totalMass = nd_mass[node] + im;
  nd_cx[node] = (nd_cx[node] * nd_mass[node] + ix * im) / totalMass;
  nd_cy[node] = (nd_cy[node] * nd_mass[node] + iy * im) / totalMass;
  nd_mass[node] = totalMass;
  const q = (ix < mx ? 0 : 1) + (iy < my ? 0 : 2);
  insertParticle(nd_child[base + q], ix, iy, im);
}

function buildTree() {
  nodeCount = 0;
  // Compute bounding box of all particles
  let x0 = px[0], y0 = py[0], x1 = px[0], y1 = py[0];
  for (let i = 1; i < N; i++) {
    if (px[i] < x0) x0 = px[i]; if (px[i] > x1) x1 = px[i];
    if (py[i] < y0) y0 = py[i]; if (py[i] > y1) y1 = py[i];
  }
  const pad = 0.001;
  const s = Math.max(x1 - x0, y1 - y0) + pad;
  allocNode(x0 - pad, y0 - pad, x0 - pad + s, y0 - pad + s);
  for (let i = 0; i < N; i++) insertParticle(0, px[i], py[i], 1.0);
}

// ── Tree walk (iterative) ─────────────────────────────────────────────────────
const _stack = new Int32Array(128);

function accelFromTree(i) {
  const xi = px[i], yi = py[i];
  const eps2 = softening * softening;
  let fx = 0, fy = 0, top = 0;
  _stack[top++] = 0;
  while (top > 0) {
    const node = _stack[--top];
    if (nd_mass[node] === 0) continue;
    const dx = nd_cx[node] - xi;
    const dy = nd_cy[node] - yi;
    const r2 = dx*dx + dy*dy;
    if (r2 < eps2) continue;
    const base = node * 4;
    const isLeaf = nd_child[base] === -1 && nd_child[base+1] === -1 &&
                   nd_child[base+2] === -1 && nd_child[base+3] === -1;
    if (isLeaf) {
      if (r2 < eps2) {
        // hard-core repulsion: constant magnitude 1/(N*ε²), pointing away
        const r = Math.sqrt(r2);
        if (r > 0) { const mag = 1.0 / (N * eps2 * r); fx -= mag * dx; fy -= mag * dy; }
      } else {
        const inv = G * nd_mass[node] / r2;
        fx += inv * dx; fy += inv * dy;
      }
    } else {
      const s = nd_x1[node] - nd_x0[node];
      if (s / Math.sqrt(r2) < theta) {
        const inv = G * nd_mass[node] / r2;
        fx += inv * dx; fy += inv * dy;
      } else {
        for (let c = 0; c < 4; c++) {
          const ch = nd_child[base + c];
          if (ch !== -1 && nd_mass[ch] > 0) _stack[top++] = ch;
        }
      }
    }
  }
  return [fx, fy];
}

// ── Init ──────────────────────────────────────────────────────────────────────
function init() {
  px = new Float64Array(N); py = new Float64Array(N);
  vx = new Float64Array(N); vy = new Float64Array(N);
  ax = new Float64Array(N); ay = new Float64Array(N);

  const v0 = Math.sqrt(0.623 * G * N) * 0.25;
  const R  = 1.0 / Math.sqrt(2.0);
  for (let i = 0; i < N; i++) {
    // uniform disk via rejection sampling
    let dx, dy;
    do { dx = Math.random() * 2 - 1; dy = Math.random() * 2 - 1; }
    while (dx*dx + dy*dy > 1.0);
    px[i] = 0.5 + dx * R; py[i] = 0.5 + dy * R;
    const angle = Math.random() * 2 * Math.PI;
    const rx = Math.cos(angle), ry = Math.sin(angle);
    // tangential direction (y, -x) relative to displacement from centre
    const tx = (py[i] - 0.5), ty = -(px[i] - 0.5);
    const tn2 = tx*tx + ty*ty || 1;
    const wx = (1 - spinFrac) * rx + spinFrac * tx/tn2;
    const wy = (1 - spinFrac) * ry + spinFrac * ty/tn2;
    const wn = Math.sqrt(wx*wx + wy*wy) || 1;
    vx[i] = v0 * wx/wn;
    vy[i] = v0 * wy/wn;
  }

  buildTree();
  for (let i = 0; i < N; i++) {
    const [fx, fy] = accelFromTree(i);
    ax[i] = fx; ay[i] = fy;
  }

  posArr = new Float32Array(N * 2);
  spdArr = new Float32Array(N);
  syncBuffers();
}

// ── Sync CPU → GPU buffers ────────────────────────────────────────────────────
const SPEED_MAX = 2.0;  // speed mapped to full pink

function syncBuffers() {
  // Compute CoM velocity
  let cvx = 0, cvy = 0;
  for (let i = 0; i < N; i++) { cvx += vx[i]; cvy += vy[i]; }
  cvx /= N; cvy /= N;

  const inv = 1.0 / SPEED_MAX;
  for (let i = 0; i < N; i++) {
    const dvx = vx[i] - cvx, dvy = vy[i] - cvy;
    const spd = Math.sqrt(dvx*dvx + dvy*dvy);
    posArr[i*2]   = px[i];
    posArr[i*2+1] = py[i];
    spdArr[i]     = Math.min(spd * inv, 1.0);
  }
}

// ── Leapfrog step ─────────────────────────────────────────────────────────────
function step() {
  const scale = 1.0 + lambda * dt;
  let cx = 0.5, cy = 0.5;
  if (expanding && followCoM) {
    for (let i = 0; i < N; i++) { cx += px[i]; cy += py[i]; }
    cx /= N; cy /= N;
  }
  for (let i = 0; i < N; i++) {
    vx[i] += 0.5 * dt * ax[i]; vy[i] += 0.5 * dt * ay[i];
    px[i] += dt * vx[i];       py[i] += dt * vy[i];
    if (expanding) {
      px[i] = cx + (px[i] - cx) * scale;
      py[i] = cy + (py[i] - cy) * scale;
    }
  }
  buildTree();
  for (let i = 0; i < N; i++) {
    const [fx, fy] = accelFromTree(i);
    ax[i] = fx; ay[i] = fy;
  }
  for (let i = 0; i < N; i++) {
    vx[i] += 0.5 * dt * ax[i]; vy[i] += 0.5 * dt * ay[i];
  }
  syncBuffers();
}

// ── WebGL shaders ─────────────────────────────────────────────────────────────
// aSpeed in [0,1]: 0 = slow (teal), 1 = fast (pink)
const VS = `
  attribute vec2  aPos;
  attribute float aSpeed;
  uniform float   uSize;
  uniform vec2    uOffset;
  uniform float   uZoom;
  varying float   vSpeed;
  void main() {
    vec2 p    = aPos - uOffset;
    // zoom about (0.5, 0.5) in particle space
    p = (p - 0.5) * uZoom + 0.5;
    gl_Position  = vec4(p * 2.0 - 1.0, 0.0, 1.0);
    gl_PointSize = uSize;
    vSpeed       = aSpeed;
  }
`;

const FS = `
  precision mediump float;
  uniform vec3 uColSlow;
  uniform vec3 uColFast;
  varying float vSpeed;
  void main() {
    vec2  d = gl_PointCoord - 0.5;
    float r = dot(d, d);
    if (r > 0.25) discard;
    float core  = clamp(1.0 - 4.0 * sqrt(r), 0.0, 1.0);
    vec3  col   = mix(uColSlow, uColFast, vSpeed);
    float alpha = core * 0.9 + 0.08;
    gl_FragColor = vec4(col, alpha);
  }
`;

function compileShader(type, src) {
  const sh = gl.createShader(type);
  gl.shaderSource(sh, src); gl.compileShader(sh);
  return sh;
}

function initGL() {
  simCanvas = document.getElementById('nb-canvas');
  gl = simCanvas.getContext('webgl', { antialias: false, alpha: false });
  if (!gl) return;

  const vs = compileShader(gl.VERTEX_SHADER,   VS);
  const fs = compileShader(gl.FRAGMENT_SHADER, FS);
  prog = gl.createProgram();
  gl.attachShader(prog, vs); gl.attachShader(prog, fs);
  gl.linkProgram(prog);
  gl.useProgram(prog);

  uPointSize = gl.getUniformLocation(prog, 'uSize');
  uColSlow   = gl.getUniformLocation(prog, 'uColSlow');
  uColFast   = gl.getUniformLocation(prog, 'uColFast');
  uOffset    = gl.getUniformLocation(prog, 'uOffset');
  uZoom      = gl.getUniformLocation(prog, 'uZoom');
  gl.uniform1f(uZoom, zoom);

  // Position buffer
  posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, posArr, gl.DYNAMIC_DRAW);
  const aPos = gl.getAttribLocation(prog, 'aPos');
  gl.enableVertexAttribArray(aPos);
  gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);

  // Speed buffer
  spdBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, spdBuf);
  gl.bufferData(gl.ARRAY_BUFFER, spdArr, gl.DYNAMIC_DRAW);
  const aSpd = gl.getAttribLocation(prog, 'aSpeed');
  gl.enableVertexAttribArray(aSpd);
  gl.vertexAttribPointer(aSpd, 1, gl.FLOAT, false, 0, 0);

  // Additive blending
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE);

  // Colours from palette
  const [sr, sg, sb] = _rgb('--teal-light');
  const [fr, fg, fb] = _rgb('--pink-light');
  gl.uniform3f(uColSlow, sr/255, sg/255, sb/255);
  gl.uniform3f(uColFast, fr/255, fg/255, fb/255);

  resizeGL();
}

function resizeGL() {
  if (!gl || !simCanvas) return;
  const w = simCanvas.clientWidth, h = simCanvas.clientHeight;
  simCanvas.width = w; simCanvas.height = h;
  gl.viewport(0, 0, w, h);
  const size = Math.max(7.2, Math.min(18.0, w / 500 * 6.0));
  gl.uniform1f(uPointSize, size);
}

function render() {
  if (!gl) return;
  const [br, bg, bb] = _rgb('--bg-void');
  gl.clearColor(br/255, bg/255, bb/255, 1);
  gl.clear(gl.COLOR_BUFFER_BIT);

  // CoM offset for camera follow — shift so CoM maps to centre of screen
  let ox = 0, oy = 0;
  if (followCoM) {
    let cx = 0, cy = 0;
    for (let i = 0; i < N; i++) { cx += px[i]; cy += py[i]; }
    ox = cx / N - 0.5; oy = cy / N - 0.5;
  }
  gl.uniform2f(uOffset, ox, oy);
  gl.uniform1f(uZoom, zoom);

  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferSubData(gl.ARRAY_BUFFER, 0, posArr);

  gl.bindBuffer(gl.ARRAY_BUFFER, spdBuf);
  gl.bufferSubData(gl.ARRAY_BUFFER, 0, spdArr);

  // Re-bind attributes after buffer switches
  const aPos = gl.getAttribLocation(prog, 'aPos');
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.vertexAttribPointer(aPos, 2, gl.FLOAT, false, 0, 0);

  const aSpd = gl.getAttribLocation(prog, 'aSpeed');
  gl.bindBuffer(gl.ARRAY_BUFFER, spdBuf);
  gl.vertexAttribPointer(aSpd, 1, gl.FLOAT, false, 0, 0);

  gl.drawArrays(gl.POINTS, 0, N);
}

// ── Animation loop ────────────────────────────────────────────────────────────
function loop() {
  if (running) {
    for (let s = 0; s < STEPS_PER_FRAME; s++) step();
    render();
    computePowerSpectrum();
  }
  frameId = requestAnimationFrame(loop);
}

// ── Shell wiring ──────────────────────────────────────────────────────────────
const shell = new AppletShell({
  id:     'nb',
  title:  'N-Body Gravity &mdash; Barnes&ndash;Hut',
  gap:    0,

  headerBtns: `<button class="applet-shell-header-btn" onclick="nbReset()">Reset</button><button class="applet-shell-header-btn active" id="nb-pause-btn" onclick="nbTogglePause()">Resume</button>`,

  ctrlHTML: `
    <div class="applet-shell-ctrl-section" style="flex:0 0 auto;">
      <div class="applet-shell-ctrl-title">Initial Spin</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">0%</span>
        <input type="range" id="nb-spin" min="0" max="1" step="0.01" value="0.2">
        <span class="applet-shell-side">100%</span>
        <span class="applet-shell-val" id="nb-spin-val">20%</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section" style="flex:0 0 auto;">
      <div class="applet-shell-ctrl-title">Expansion</div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" id="nb-expand-btn" onclick="nbToggleExpand()">Expand</button>
      </div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">0</span>
        <input type="range" id="nb-lambda" min="0" max="1" step="0.01" value="0">
        <span class="applet-shell-side">1</span>
        <span class="applet-shell-val" id="nb-lambda-val">0.00</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section" style="flex:0 0 auto;">
      <div class="applet-shell-ctrl-title">View</div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" id="nb-follow-btn" onclick="nbToggleFollow()">Follow CoM</button>
      </div>
    </div>
    <div class="applet-shell-ctrl-section" style="flex:0 0 auto;">
      <div class="applet-shell-ctrl-title">Tree &theta;</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Exact</span>
        <input type="range" id="nb-theta" min="0.2" max="1.2" step="0.05" value="0.3">
        <span class="applet-shell-side">Fast</span>
        <span class="applet-shell-val" id="nb-theta-val">0.30</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section nb-pspec-section">
      <div class="applet-shell-ctrl-title">Two-Point Correlation &xi;(r)</div>
      <canvas id="nb-xi"></canvas>
    </div>
    <div class="applet-shell-ctrl-section nb-pspec-section">
      <div class="applet-shell-ctrl-title">Power Spectrum P(k)</div>
      <canvas id="nb-pspec"></canvas>
    </div>
  `,

  onOpen: function ({ canvas: c, S, W, H }) {
    // Start paused
    running = false;
    const pb = document.getElementById('nb-pause-btn');
    if (pb) { pb.textContent = 'Resume'; pb.classList.add('active'); }

    zoom = 1.0;
    frameCount = 0;
    init();
    setTimeout(() => {
      initGL();
      xiCanvas = document.getElementById('nb-xi');
      xiCtx    = xiCanvas.getContext('2d');
      const xiSection = xiCanvas.closest('.nb-pspec-section');
      xiCanvas.width  = xiSection.clientWidth  - 24;
      xiCanvas.height = xiSection.clientHeight - 28;

      pspecCanvas = document.getElementById('nb-pspec');
      pspecCtx    = pspecCanvas.getContext('2d');
      const pspecSection = pspecCanvas.closest('.nb-pspec-section');
      pspecCanvas.width  = pspecSection.clientWidth  - 24;
      pspecCanvas.height = pspecSection.clientHeight - 28;
      simCanvas.addEventListener('wheel', function (e) {
        e.preventDefault();
        zoom *= e.deltaY > 0 ? 0.98 : 1.02;
        zoom = Math.max(0.1, Math.min(20.0, zoom));
        if (!running) render();
      }, { passive: false });
      render(); // draw initial state while paused
      if (!frameId) frameId = requestAnimationFrame(loop);
    }, 80);
  },

  onClose: function () {
    running = false;
    if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
  },

  onResize: function () {
    resizeGL();
    if (!running) render();
  },
});

// ── Global entry points ───────────────────────────────────────────────────────
window.nbOpen  = () => shell.open();
window.nbClose = () => shell.close();

window.nbReset = function () {
  running = false;
  const pb = document.getElementById('nb-pause-btn');
  if (pb) { pb.textContent = 'Resume'; pb.classList.add('active'); }
  init();
  if (gl) {
    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.bufferData(gl.ARRAY_BUFFER, posArr, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, spdBuf);
    gl.bufferData(gl.ARRAY_BUFFER, spdArr, gl.DYNAMIC_DRAW);
  }
  render();
};

window.nbToggleExpand = function () {
  expanding = !expanding;
  const btn = document.getElementById('nb-expand-btn');
  if (btn) btn.classList.toggle('active', expanding);
};

window.nbToggleFollow = function () {
  followCoM = !followCoM;
  const btn = document.getElementById('nb-follow-btn');
  if (btn) btn.classList.toggle('active', followCoM);
  if (!running) render();
};

window.nbTogglePause = function () {
  running = !running;
  const pb = document.getElementById('nb-pause-btn');
  if (pb) {
    pb.textContent = running ? 'Pause' : 'Resume';
    pb.classList.toggle('active', !running);
  }
};

// ── Slider listeners ──────────────────────────────────────────────────────────
document.getElementById('nb-spin').addEventListener('input', function () {
  spinFrac = parseFloat(this.value);
  document.getElementById('nb-spin-val').textContent = Math.round(spinFrac * 100) + '%';
});

document.getElementById('nb-lambda').addEventListener('input', function () {
  lambda = parseFloat(this.value);
  document.getElementById('nb-lambda-val').textContent = lambda.toFixed(2);
});

document.getElementById('nb-theta').addEventListener('input', function () {
  theta = parseFloat(this.value);
  document.getElementById('nb-theta-val').textContent = theta.toFixed(2);
});

})();

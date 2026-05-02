(function () {
'use strict';

/* ── Palette helpers ── */
const _cs  = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();
const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

/* ══════════════════════════════════════════════════════════════════
   ANT COLONY SIMULATION
   ──────────────────────────────────────────────────────────────────
   Pheromone field: 2D grid (GW × GH cells), scalar concentration.
   Each frame: diffuse (5-point Laplacian), evaporate.

   Ants: N particles, each with (x, y, angle, state).
   States:
     'search'   — biased random walk toward pheromone gradient,
                  deposit pheromone while returning.
     'return'   — head toward nest, deposit pheromone.
   On touching food: switch to 'return', pick up food.
   On reaching nest: switch to 'search', deposit food.

   Barriers: boolean mask on same grid. Ants reflect off them.
   Food: scalar amount per grid cell.
══════════════════════════════════════════════════════════════════ */

/* ── Grid dimensions ── */
const GW = 256, GH = 256;   // pheromone/barrier grid cells

/* ── Simulation parameters (defaults) ── */
let N_ANTS        = 200;
let EVAP_RATE     = 0.002;   // fraction lost per frame
let DIFF_RATE     = 0.15;    // diffusion coefficient per frame
let PHER_STRENGTH = 8.0;     // pheromone deposited per step
let stepsPerFrame = 2;

/* ── Ant movement ── */
const ANT_SPEED    = 0.8;    // cells per step
const SENSE_DIST   = 3.0;    // how far ahead ant sniffs (cells)
const SENSE_ANGLE  = 0.4;    // half-angle of sensor spread (rad)
let TURN_NOISE     = 0.5;    // random turn noise (rad)
const TURN_SPEED   = 0.35;   // max steering turn per step (rad)
const NEST_RADIUS  = 6;      // cells

/* ── State arrays ── */
let pheromone   = new Float32Array(GW * GH);
let nestPot     = new Float32Array(GW * GH);  // precomputed nest potential V(x,y)
let nestGradX   = new Float32Array(GW * GH);  // ∂V/∂x
let nestGradY   = new Float32Array(GW * GH);  // ∂V/∂y
let barriers    = new Uint8Array(GW * GH);
let food        = new Float32Array(GW * GH);
let ants        = [];
let nestX       = GW / 2, nestY = GH / 2;
let totalFood   = 0;       // food delivered to nest

const NEST_LAMBDA = GW / 5;   // decay length scale (~50 cells)
const NEST_WEIGHT = 3.0;       // relative weight of nest gradient vs pheromone

/* ── Build nest potential field (call after nestX/nestY changes) ── */
function buildNestPotential() {
  for (let y = 0; y < GH; y++) {
    for (let x = 0; x < GW; x++) {
      const dx = x - nestX, dy = y - nestY;
      const r = Math.sqrt(dx*dx + dy*dy);
      nestPot[y * GW + x] = Math.exp(-r / NEST_LAMBDA);
    }
  }
  // Finite-difference gradient (central differences, clamp at edges)
  for (let y = 0; y < GH; y++) {
    for (let x = 0; x < GW; x++) {
      const xp = x < GW-1 ? nestPot[y*GW + x+1] : nestPot[y*GW + x];
      const xm = x > 0    ? nestPot[y*GW + x-1] : nestPot[y*GW + x];
      const yp = y < GH-1 ? nestPot[(y+1)*GW + x] : nestPot[y*GW + x];
      const ym = y > 0    ? nestPot[(y-1)*GW + x] : nestPot[y*GW + x];
      nestGradX[y*GW + x] = (xp - xm) * 0.5;
      nestGradY[y*GW + x] = (yp - ym) * 0.5;
    }
  }
}

/* ── Texture canvases ── */
let nestTexCanvas = null;
let wallTexCanvas = null;   // tiled roughness stamp

/* ── Canvas / rendering ── */
let canvas, ctx, S;
let offscreen = null, offCtx = null;
let imgData = null, buf = null;

/* ── Interaction mode ── */
// 'none' | 'barrier' | 'food' | 'nest'
let interactMode = 'none';
let isDrawing = false;
const BRUSH_R = 4;          // brush radius in grid cells

/* ── Running flag ── */
let running = false, frameId = null;

/* ── Noise helper (simple value noise for textures) ── */
function valueNoise(x, y, seed) {
  const n = Math.sin(x * 127.1 + y * 311.7 + seed * 74.3) * 43758.5453;
  return n - Math.floor(n);
}

/* ── Build nest texture (rough circle) ── */
function buildNestTex() {
  nestTexCanvas = document.createElement('canvas');
  nestTexCanvas.width = nestTexCanvas.height = NEST_RADIUS * 2 + 4;
  const c = nestTexCanvas.getContext('2d');
  const R = NEST_RADIUS + 1;
  const [br, bg, bb] = _rgb('--bg-mid');
  const [tr, tg, tb] = _rgb('--border-mid');
  for (let py = 0; py < nestTexCanvas.height; py++) {
    for (let px = 0; px < nestTexCanvas.width; px++) {
      const dx = px - R, dy = py - R;
      const dist = Math.sqrt(dx*dx + dy*dy);
      const edge = dist - (NEST_RADIUS - 0.5);
      if (edge > 1.5) continue;
      const noise = valueNoise(px * 0.7, py * 0.7, 3) * 0.18 - 0.09;
      const t = Math.max(0, Math.min(1, 1 - edge + noise));
      const bright = 0.85 + valueNoise(px * 1.3, py * 1.3, 7) * 0.15;
      c.fillStyle = `rgba(${Math.round(br*bright)},${Math.round(bg*bright)},${Math.round(bb*bright)},${t})`;
      c.fillRect(px, py, 1, 1);
    }
  }
  // ring highlight
  c.strokeStyle = _rgba('--border-mid', 0.6);
  c.lineWidth = 1.5;
  c.beginPath();
  c.arc(R, R, NEST_RADIUS - 1, 0, Math.PI * 2);
  c.stroke();
}

/* ── Init ── */
function initSim() {
  pheromone.fill(0);
  barriers.fill(0);
  food.fill(0);
  totalFood = 0;
  nestX = GW / 2; nestY = GH / 2;
  buildNestPotential();
  spawnAnts();
  buildNestTex();
}

function spawnAnts() {
  ants = [];
  for (let i = 0; i < N_ANTS; i++) {
    const angle = Math.random() * Math.PI * 2;
    ants.push({
      x:     nestX + Math.cos(angle) * 2,
      y:     nestY + Math.sin(angle) * 2,
      angle: angle,
      state: 'search',   // 'search' or 'return'
    });
  }
}

/* ── Pheromone diffusion + evaporation ── */
const pherTmp = new Float32Array(GW * GH);
function updatePheromone() {
  const dr = DIFF_RATE;
  const ev = 1 - EVAP_RATE;
  for (let y = 0; y < GH; y++) {
    for (let x = 0; x < GW; x++) {
      const i = y * GW + x;
      if (barriers[i]) { pheromone[i] = 0; continue; }
      const n = (y > 0      && !barriers[(y-1)*GW+x]) ? pheromone[(y-1)*GW+x] : pheromone[i];
      const s = (y < GH-1   && !barriers[(y+1)*GW+x]) ? pheromone[(y+1)*GW+x] : pheromone[i];
      const ww= (x > 0      && !barriers[y*GW+x-1])   ? pheromone[y*GW+x-1]   : pheromone[i];
      const e = (x < GW-1   && !barriers[y*GW+x+1])   ? pheromone[y*GW+x+1]   : pheromone[i];
      pherTmp[i] = ev * (pheromone[i] + dr * (n + s + ww + e - 4 * pheromone[i]));
    }
  }
  pheromone.set(pherTmp);
}

/* ── Sense pheromone at a point ── */
function samplePher(x, y) {
  const ix = Math.round(x), iy = Math.round(y);
  if (ix < 0 || ix >= GW || iy < 0 || iy >= GH) return 0;
  if (barriers[iy * GW + ix]) return 0;
  return pheromone[iy * GW + ix];
}

/* ── Pheromone gradient at (x,y) via central differences ── */
function pherGrad(x, y) {
  const ix = Math.round(x), iy = Math.round(y);
  const xp = Math.min(GW-1, ix+1), xm = Math.max(0, ix-1);
  const yp = Math.min(GH-1, iy+1), ym = Math.max(0, iy-1);
  const gx = (samplePher(xp, iy) - samplePher(xm, iy)) * 0.5;
  const gy = (samplePher(ix, yp) - samplePher(ix, ym)) * 0.5;
  return [gx, gy];
}

/* ── Nest gradient at (x,y) ── */
function nestGrad(x, y) {
  const ix = Math.min(GW-1, Math.max(0, Math.round(x)));
  const iy = Math.min(GH-1, Math.max(0, Math.round(y)));
  return [nestGradX[iy*GW+ix], nestGradY[iy*GW+ix]];
}

/* ── Sample food at grid cell ── */
function sampleFood(x, y) {
  const ix = Math.round(x), iy = Math.round(y);
  if (ix < 0 || ix >= GW || iy < 0 || iy >= GH) return 0;
  return food[iy * GW + ix];
}

/* ── Check barrier ── */
function isBarrier(x, y) {
  const ix = Math.round(x), iy = Math.round(y);
  if (ix < 0 || ix >= GW || iy < 0 || iy >= GH) return true;
  return barriers[iy * GW + ix] > 0;
}

/* ── Deposit pheromone ── */
function depositPher(x, y, amount) {
  const ix = Math.round(x), iy = Math.round(y);
  if (ix < 0 || ix >= GW || iy < 0 || iy >= GH) return;
  if (barriers[iy * GW + ix]) return;
  pheromone[iy * GW + ix] = Math.min(255, pheromone[iy * GW + ix] + amount);
}

/* ── Ant step ── */
function stepAnt(a) {
  const dx = nestX - a.x, dy = nestY - a.y;
  const dNest = Math.sqrt(dx*dx + dy*dy);

  // Compute combined gradient bias
  const [pgx, pgy] = pherGrad(a.x, a.y);
  const [ngx, ngy] = nestGrad(a.x, a.y);

  let biasX, biasY;
  if (a.state === 'return') {
    // Climb nest potential (toward nest) + follow pheromone
    biasX = pgx + NEST_WEIGHT * 10 * ngx;
    biasY = pgy + NEST_WEIGHT * 10 * ngy;
    depositPher(a.x, a.y, PHER_STRENGTH);
    if (dNest < NEST_RADIUS) {
      a.state = 'search';
      totalFood++;
    }
  } else {
    // Descend nest potential (away from nest) + follow pheromone
    biasX = pgx - NEST_WEIGHT * ngx;
    biasY = pgy - NEST_WEIGHT * ngy;

    // Check for food nearby
    const fx = sampleFood(a.x + Math.cos(a.angle) * 1.5, a.y + Math.sin(a.angle) * 1.5);
    if (fx > 0) {
      const fix = Math.round(a.x + Math.cos(a.angle) * 1.5);
      const fiy = Math.round(a.y + Math.sin(a.angle) * 1.5);
      if (fix >= 0 && fix < GW && fiy >= 0 && fiy < GH && food[fiy*GW+fix] > 0) {
        food[fiy*GW+fix] = Math.max(0, food[fiy*GW+fix] - 1);
        a.state = 'return';
        a.angle += Math.PI;
      }
    }
  }

  // Steer toward bias direction if it has meaningful magnitude
  const bMag = Math.sqrt(biasX*biasX + biasY*biasY);
  if (bMag > 1e-4) {
    const biasAngle = Math.atan2(biasY, biasX);
    let da = biasAngle - a.angle;
    while (da >  Math.PI) da -= Math.PI * 2;
    while (da < -Math.PI) da += Math.PI * 2;
    a.angle += Math.max(-TURN_SPEED, Math.min(TURN_SPEED, da));
  }

  // Add random noise
  a.angle += (Math.random() - 0.5) * TURN_NOISE;

  // Move
  const nx = a.x + Math.cos(a.angle) * ANT_SPEED;
  const ny = a.y + Math.sin(a.angle) * ANT_SPEED;

  if (isBarrier(nx, ny)) {
    a.angle += Math.PI * (0.5 + Math.random() * 0.5);
  } else {
    a.x = Math.max(0, Math.min(GW - 1, nx));
    a.y = Math.max(0, Math.min(GH - 1, ny));
  }
}

/* ── Simulation step ── */
function simStep() {
  updatePheromone();
  for (let i = 0; i < ants.length; i++) stepAnt(ants[i]);
}

/* ── Color LUTs ── */
let lutPher = null, lutFood = null;

function buildLUTs() {
  const [pr, pg, pb] = _rgb('--amber');
  const [fr, fg, fb] = _rgb('--green-light');
  const [br, bg, bb] = _rgb('--bg-deep');
  const [wr, wg, wb] = _rgb('--bg-mid');
  const [brr, brg, brb] = _rgb('--border-dark');

  lutPher = new Uint8Array(256 * 4);
  for (let i = 0; i < 256; i++) {
    const t = i > 0 ? Math.min(1, Math.log(1 + i) / Math.log(32)) : 0;
    lutPher[i*4+0] = Math.round(br + (pr - br) * t);
    lutPher[i*4+1] = Math.round(bg + (pg - bg) * t);
    lutPher[i*4+2] = Math.round(bb + (pb - bb) * t);
    lutPher[i*4+3] = 255;
  }

  lutFood = new Uint8Array(256 * 4);
  for (let i = 0; i < 256; i++) {
    const t = Math.min(1, i / 30);
    lutFood[i*4+0] = Math.round(fr * t);
    lutFood[i*4+1] = Math.round(fg * t);
    lutFood[i*4+2] = Math.round(fb * t);
    lutFood[i*4+3] = i > 0 ? 220 : 0;
  }
}

/* ── Render ── */
function render() {
  if (!canvas || !ctx) return;

  const scaleX = S / GW, scaleY = S / GH;

  // Build pixel buffer
  const data = buf;
  const [wr, wg, wb] = _rgb('--bg-mid');
  const [brr, brg, brb] = _rgb('--border-dark');

  for (let y = 0; y < GH; y++) {
    for (let x = 0; x < GW; x++) {
      const i = y * GW + x;
      const base = i * 4;
      if (barriers[i]) {
        // Rough wall texture
        const noise = valueNoise(x * 0.5, y * 0.5, 1);
        const bright = 0.8 + noise * 0.2;
        data[base]   = Math.round(wr * bright);
        data[base+1] = Math.round(wg * bright);
        data[base+2] = Math.round(wb * bright);
        data[base+3] = 255;
      } else {
        const p = Math.min(255, pheromone[i]) | 0;
        data[base]   = lutPher[p*4];
        data[base+1] = lutPher[p*4+1];
        data[base+2] = lutPher[p*4+2];
        data[base+3] = 255;
        // Overlay food
        if (food[i] > 0) {
          const f = Math.min(255, food[i] * 5) | 0;
          const ta = lutFood[f*4+3] / 255;
          data[base]   = Math.round(data[base]   * (1-ta) + lutFood[f*4]   * ta);
          data[base+1] = Math.round(data[base+1] * (1-ta) + lutFood[f*4+1] * ta);
          data[base+2] = Math.round(data[base+2] * (1-ta) + lutFood[f*4+2] * ta);
        }
      }
    }
  }

  // Draw boundaries as barriers on edges
  for (let x = 0; x < GW; x++) {
    const noise0 = valueNoise(x * 0.5, 0, 2);
    const noise1 = valueNoise(x * 0.5, GH-1, 2);
    const bright0 = 0.8 + noise0 * 0.2, bright1 = 0.8 + noise1 * 0.2;
    const b0 = 0 * GW + x, b1 = (GH-1) * GW + x;
    buf[b0*4]   = Math.round(wr*bright0); buf[b0*4+1] = Math.round(wg*bright0); buf[b0*4+2] = Math.round(wb*bright0); buf[b0*4+3] = 255;
    buf[b1*4]   = Math.round(wr*bright1); buf[b1*4+1] = Math.round(wg*bright1); buf[b1*4+2] = Math.round(wb*bright1); buf[b1*4+3] = 255;
  }
  for (let y = 0; y < GH; y++) {
    const noise0 = valueNoise(0, y * 0.5, 2);
    const noise1 = valueNoise(GW-1, y * 0.5, 2);
    const bright0 = 0.8 + noise0 * 0.2, bright1 = 0.8 + noise1 * 0.2;
    const b0 = y * GW + 0, b1 = y * GW + (GW-1);
    buf[b0*4]   = Math.round(wr*bright0); buf[b0*4+1] = Math.round(wg*bright0); buf[b0*4+2] = Math.round(wb*bright0); buf[b0*4+3] = 255;
    buf[b1*4]   = Math.round(wr*bright1); buf[b1*4+1] = Math.round(wg*bright1); buf[b1*4+2] = Math.round(wb*bright1); buf[b1*4+3] = 255;
  }

  // Blit pixel buffer to offscreen, then scale to main canvas
  offCtx.putImageData(imgData, 0, 0);

  // Draw nest texture onto offscreen (grid space)
  if (nestTexCanvas) {
    const r = NEST_RADIUS + 2;
    offCtx.drawImage(nestTexCanvas, nestX - r, nestY - r);
  }

  // Scale offscreen to main canvas
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(offscreen, 0, 0, S, S);

  // Draw ants in screen space
  for (let i = 0; i < ants.length; i++) {
    const a = ants[i];
    const sx = a.x * scaleX, sy = a.y * scaleY;
    const len = 3.5;
    const col = a.state === 'return' ? _rgba('--teal-light', 0.9) : _rgba('--pink-light', 0.85);
    ctx.save();
    ctx.translate(sx, sy);
    ctx.rotate(a.angle);
    ctx.fillStyle = col;
    ctx.beginPath();
    ctx.ellipse(0, 0, len, 1.2, 0, 0, Math.PI * 2);
    ctx.fill();
    ctx.restore();
  }

  // Draw nest circle in screen space
  const nsx = nestX * scaleX, nsy = nestY * scaleY;
  const nsr = NEST_RADIUS * scaleX;
  ctx.save();
  ctx.strokeStyle = _rgba('--border-mid', 0.9);
  ctx.lineWidth = 1.5;
  // Rough circle: draw as many small arcs with jitter
  ctx.beginPath();
  for (let i = 0; i <= 60; i++) {
    const angle = (i / 60) * Math.PI * 2;
    const jitter = 1 + (valueNoise(Math.cos(angle) * 5, Math.sin(angle) * 5, 5) - 0.5) * 0.12;
    const rx = nsx + Math.cos(angle) * nsr * jitter;
    const ry = nsy + Math.sin(angle) * nsr * jitter;
    if (i === 0) ctx.moveTo(rx, ry); else ctx.lineTo(rx, ry);
  }
  ctx.closePath();
  ctx.strokeStyle = _rgba('--border-mid', 0.7);
  ctx.lineWidth = 1.5;
  ctx.stroke();
  // Filled with noisy color
  ctx.fillStyle = _rgba('--bg-mid', 0.55);
  ctx.fill();
  ctx.restore();

  // Food count readout
  ctx.fillStyle = _rgba('--text-dim', 0.8);
  ctx.font = `${Math.round(S * 0.03)}px monospace`;
  ctx.fillText(`Food: ${totalFood}`, 8, S - 8);
}

/* ── Animation loop ── */
function loop() {
  if (running) {
    for (let i = 0; i < stepsPerFrame; i++) simStep();
  }
  render();
  frameId = requestAnimationFrame(loop);
}

/* ── Interaction: convert canvas coords to grid coords ── */
function canvasToGrid(cx, cy) {
  return [
    Math.round(cx / (S / GW)),
    Math.round(cy / (S / GH)),
  ];
}

function paintAt(gx, gy) {
  for (let dy = -BRUSH_R; dy <= BRUSH_R; dy++) {
    for (let dx = -BRUSH_R; dx <= BRUSH_R; dx++) {
      if (dx*dx + dy*dy > BRUSH_R*BRUSH_R) continue;
      const ix = gx + dx, iy = gy + dy;
      if (ix < 1 || ix >= GW-1 || iy < 1 || iy >= GH-1) continue;
      const idx = iy * GW + ix;
      if (interactMode === 'barrier') {
        barriers[idx] = 1;
        pheromone[idx] = 0;
      } else if (interactMode === 'food') {
        food[idx] = Math.min(255, food[idx] + 20);
      }
    }
  }
}

function onPointerDown(e) {
  if (interactMode === 'none') return;
  e.preventDefault();
  isDrawing = true;
  const rect = canvas.getBoundingClientRect();
  const cx = (e.clientX - rect.left) * (S / rect.width);
  const cy = (e.clientY - rect.top)  * (S / rect.height);
  const [gx, gy] = canvasToGrid(cx, cy);
  if (interactMode === 'nest') {
    nestX = Math.max(NEST_RADIUS+1, Math.min(GW-NEST_RADIUS-1, gx));
    nestY = Math.max(NEST_RADIUS+1, Math.min(GH-NEST_RADIUS-1, gy));
    buildNestPotential();
    buildNestTex();
    spawnAnts();
    setInteractMode('none');
    return;
  }
  paintAt(gx, gy);
}

function onPointerMove(e) {
  if (!isDrawing || interactMode === 'none' || interactMode === 'nest') return;
  e.preventDefault();
  const rect = canvas.getBoundingClientRect();
  const cx = (e.clientX - rect.left) * (S / rect.width);
  const cy = (e.clientY - rect.top)  * (S / rect.height);
  const [gx, gy] = canvasToGrid(cx, cy);
  paintAt(gx, gy);
}

function onPointerUp(e) {
  isDrawing = false;
}

function attachEvents() {
  canvas.addEventListener('pointerdown', onPointerDown);
  canvas.addEventListener('pointermove', onPointerMove);
  canvas.addEventListener('pointerup',   onPointerUp);
  canvas.addEventListener('pointerleave',onPointerUp);
}

function detachEvents() {
  if (!canvas) return;
  canvas.removeEventListener('pointerdown', onPointerDown);
  canvas.removeEventListener('pointermove', onPointerMove);
  canvas.removeEventListener('pointerup',   onPointerUp);
  canvas.removeEventListener('pointerleave',onPointerUp);
}

/* ── Mode button helpers ── */
function setInteractMode(mode) {
  interactMode = mode;
  const modes = ['barrier', 'food', 'nest'];
  modes.forEach(m => {
    const btn = document.getElementById(`ants-btn-${m}`);
    if (!btn) return;
    if (m === mode) btn.classList.add('active');
    else btn.classList.remove('active');
  });
  // cursor
  if (canvas) canvas.style.cursor = mode === 'none' ? 'default' : 'crosshair';
}

/* ── Shell wiring ── */
const shell = new AppletShell({
  id:    'ants',
  title: 'Ant Colony',
  gap:   0,

  ctrlHTML: `
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Environment</div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" id="ants-btn-barrier" onclick="antsSetMode('barrier')">Draw Barriers</button>
      </div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" id="ants-btn-food" onclick="antsSetMode('food')">Drop Food</button>
      </div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" id="ants-btn-nest" onclick="antsSetMode('nest')">Move Nest</button>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Simulation</div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" id="ants-btn-run" onclick="antsToggleRun()">Run</button>
        <button class="applet-shell-btn" onclick="antsReset()">Reset</button>
      </div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" onclick="antsResetColony()">Reset Colony</button>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Colony Size</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Few</span>
        <input type="range" id="ants-sl-n" min="50" max="600" step="10" value="200">
        <span class="applet-shell-side">Many</span>
        <span class="applet-shell-val" id="ants-sl-n-val">200</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Evaporation</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Slow</span>
        <input type="range" id="ants-sl-evap" min="0.0005" max="0.01" step="0.0005" value="0.002">
        <span class="applet-shell-side">Fast</span>
        <span class="applet-shell-val" id="ants-sl-evap-val">0.0020</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Diffusion</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Sharp</span>
        <input type="range" id="ants-sl-diff" min="0.02" max="0.4" step="0.01" value="0.15">
        <span class="applet-shell-side">Spread</span>
        <span class="applet-shell-val" id="ants-sl-diff-val">0.15</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Pheromone Strength</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Weak</span>
        <input type="range" id="ants-sl-pher" min="1" max="30" step="0.5" value="8">
        <span class="applet-shell-side">Strong</span>
        <span class="applet-shell-val" id="ants-sl-pher-val">8.0</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Random Walk Noise</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Low</span>
        <input type="range" id="ants-sl-noise" min="0.05" max="2.0" step="0.05" value="0.5">
        <span class="applet-shell-side">High</span>
        <span class="applet-shell-val" id="ants-sl-noise-val">0.50</span>
      </div>
    </div>
  `,

  onOpen: function ({ canvas: c, S: s }) {
    canvas = c; S = s;
    ctx = canvas.getContext('2d');
    offscreen = document.createElement('canvas');
    offscreen.width = GW; offscreen.height = GH;
    offCtx = offscreen.getContext('2d');
    imgData = offCtx.createImageData(GW, GH);
    buf = imgData.data;
    buildLUTs();
    initSim();
    attachEvents();
    running = false;
    const btn = document.getElementById('ants-btn-run');
    if (btn) { btn.textContent = 'Run'; btn.classList.remove('active'); }
    if (!frameId) loop();
  },

  onClose: function () {
    running = false;
    detachEvents();
    if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    setInteractMode('none');
  },

  onResize: function ({ canvas: c, S: s }) {
    canvas = c; S = s;
    ctx = canvas.getContext('2d');
  },
});

/* ── Global entry points ── */
window.antsOpen  = () => shell.open();
window.antsClose = () => shell.close();

window.antsResetColony = function () {
  spawnAnts();
  pheromone.fill(0);
  totalFood = 0;
};

window.antsReset = function () {
  initSim();
  running = false;
  setInteractMode('none');
  const btn = document.getElementById('ants-btn-run');
  if (btn) { btn.textContent = 'Run'; btn.classList.remove('active'); }
};

window.antsToggleRun = function () {
  running = !running;
  const btn = document.getElementById('ants-btn-run');
  if (btn) {
    btn.textContent = running ? 'Pause' : 'Run';
    btn.classList.toggle('active', running);
  }
};

window.antsSetMode = function (mode) {
  setInteractMode(interactMode === mode ? 'none' : mode);
};

/* ── Slider listeners ── */
document.getElementById('ants-sl-n').addEventListener('input', function () {
  N_ANTS = parseInt(this.value);
  document.getElementById('ants-sl-n-val').textContent = N_ANTS;
  spawnAnts();
});
document.getElementById('ants-sl-evap').addEventListener('input', function () {
  EVAP_RATE = parseFloat(this.value);
  document.getElementById('ants-sl-evap-val').textContent = EVAP_RATE.toFixed(4);
});
document.getElementById('ants-sl-diff').addEventListener('input', function () {
  DIFF_RATE = parseFloat(this.value);
  document.getElementById('ants-sl-diff-val').textContent = DIFF_RATE.toFixed(2);
});
document.getElementById('ants-sl-pher').addEventListener('input', function () {
  PHER_STRENGTH = parseFloat(this.value);
  document.getElementById('ants-sl-pher-val').textContent = PHER_STRENGTH.toFixed(1);
});
document.getElementById('ants-sl-noise').addEventListener('input', function () {
  TURN_NOISE = parseFloat(this.value);
  document.getElementById('ants-sl-noise-val').textContent = TURN_NOISE.toFixed(2);
});

})();

(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  /* ── Inject CSS ── */
  (function () {
    if (document.getElementById('bm-styles')) return;
    const s = document.createElement('style');
    s.id = 'bm-styles';
    s.textContent = `
      #bm-ctrl-panel {
        display: flex;
        flex-direction: column;
        overflow: hidden;
      }
      #bm-sliders {
        flex: 1;
        overflow-y: auto;
        min-height: 0;
      }
      #bm-hist-section {
        flex: 0 0 33%;
        min-height: 0;
        display: flex;
        flex-direction: column;
        padding: 8px 12px 10px;
        border-top: 1px solid var(--border-dark);
      }
      #bm-hist-section .applet-shell-ctrl-title {
        margin-bottom: 4px;
        flex-shrink: 0;
      }
      #bm-hist-canvas {
        flex: 1;
        min-height: 0;
        width: 100%;
        display: block;
      }
    `;
    document.head.appendChild(s);
  })();

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');
  const [_PLR, _PLG, _PLB] = _rgb('--pink-light');
  const [_TDR, _TDG, _TDB] = _rgb('--teal-dark');
  const [_AMR, _AMG, _AMB] = _rgb('--amber');

  /* ── Simulation parameters ── */
  let N   = 500;   // small particles
  let T   = 1.0;
  let rs  = 5;     // small particle radius
  let L   = 500;   // box side (updated on open/resize)

  /* Big particle */
  const RB_FRAC = 0.06;   // big radius as fraction of L
  let Rb, Mb;
  let bx, by, bvx, bvy;
  let trail = [];
  const TRAIL_MAX = 600;

  /* Small particles */
  let px, py, vx, vy;
  let ms;   // mass of small particle

  let running = false, frameId = null;
  let canvas, ctx, hctx;
  let gravity = 0;

  /* ── Cell list ── */
  let cellW, nCX, nCY, cellHead, cellNext;
  const CELL_SLOTS = 2;   // cell width = CELL_SLOTS * max_diameter

  function buildCells() {
    const diam = Math.max(2 * rs, 2 * Rb, 1);
    cellW = diam * CELL_SLOTS;
    nCX   = Math.max(Math.floor(L / cellW), 1);
    nCY   = Math.max(Math.floor(L / cellW), 1);
    const nc = nCX * nCY;
    cellHead = new Int32Array(nc).fill(-1);
    // slots: 0..N-1 are small, N is big particle
    cellNext = new Int32Array(N + 1).fill(-1);
    const cW = L / nCX, cH = L / nCY;
    for (let i = 0; i < N; i++) {
      const cx = Math.min((px[i] / cW) | 0, nCX - 1);
      const cy = Math.min((py[i] / cH) | 0, nCY - 1);
      const c  = cy * nCX + cx;
      cellNext[i] = cellHead[c]; cellHead[c] = i;
    }
    // big particle index = N
    const cx = Math.min((bx / cW) | 0, nCX - 1);
    const cy = Math.min((by / cH) | 0, nCY - 1);
    const c  = cy * nCX + cx;
    cellNext[N] = cellHead[c]; cellHead[c] = N;
    buildCells._cW = cW; buildCells._cH = cH;
  }

  /* ── Random helpers ── */
  function gaussRand() {
    let u; do { u = Math.random(); } while (u === 0);
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * Math.random());
  }

  /* ── Init ── */
  function initSim() {
    Rb = Math.round(RB_FRAC * L);
    ms = rs * rs;
    Mb = Rb * Rb * 0.5;   // moderately heavier than small particles

    /* Place big particle in centre */
    bx = L / 2; by = L / 2; bvx = 0; bvy = 0;
    trail = [];

    /* Place small particles on a grid, avoiding the big one */
    px  = new Float64Array(N); py  = new Float64Array(N);
    vx  = new Float64Array(N); vy  = new Float64Array(N);
    const margin = rs + 1;
    const innerW = L - 2 * margin;
    const innerH = L - 2 * margin;
    const cols = Math.max(Math.ceil(Math.sqrt(N * innerW / Math.max(innerH, 1))), 1);
    const rows = Math.ceil(N / cols);
    const dx = innerW / cols, dy = innerH / rows;
    let placed = 0;
    outer: for (let r = 0; r < rows && placed < N; r++) {
      for (let c = 0; c < cols && placed < N; c++) {
        const x0 = margin + (c + 0.5) * dx;
        const y0 = margin + (r + 0.5) * dy;
        /* skip if overlaps big particle */
        const ddx = x0 - bx, ddy = y0 - by;
        if (ddx*ddx + ddy*ddy < (rs + Rb + 2) * (rs + Rb + 2)) continue;
        /* skip if overlaps already-placed small particle */
        let ok = true;
        for (let j = 0; j < placed; j++) {
          const ex = x0 - px[j], ey = y0 - py[j];
          if (ex*ex + ey*ey < (2*rs)*(2*rs)) { ok = false; break; }
        }
        if (!ok) continue;
        px[placed] = x0; py[placed] = y0;
        vx[placed] = gaussRand() * Math.sqrt(T);
        vy[placed] = gaussRand() * Math.sqrt(T);
        placed++;
      }
    }
    N = placed;

    /* Zero net momentum, rescale to T */
    let svx = 0, svy = 0;
    for (let i = 0; i < N; i++) { svx += vx[i]; svy += vy[i]; }
    svx /= N; svy /= N;
    for (let i = 0; i < N; i++) { vx[i] -= svx; vy[i] -= svy; }
    let ke = 0;
    for (let i = 0; i < N; i++) ke += vx[i]*vx[i] + vy[i]*vy[i];
    ke /= (2 * N);
    const sc = Math.sqrt(T / Math.max(ke, 1e-12));
    for (let i = 0; i < N; i++) { vx[i] *= sc; vy[i] *= sc; }
  }

  /* ── Exact collision between small particle i and big particle ── */
  function collideBig(i) {
    const ddx = bx - px[i], ddy = by - py[i];
    const d2  = ddx*ddx + ddy*ddy;
    const md  = rs + Rb;
    if (d2 >= md * md || d2 < 1e-12) return;

    /* Quadratic backtrack to exact contact */
    const sdx = vx[i] - bvx, sdy = vy[i] - bvy;
    const sx  = px[i] - bx,  sy  = py[i] - by;
    const a   = sdx*sdx + sdy*sdy;
    if (a < 1e-14) return;
    const b   = 2 * (sx*sdx + sy*sdy);
    const c   = sx*sx + sy*sy - md*md;
    const disc = b*b - 4*a*c;
    if (disc < 0) return;
    let t = (-b - Math.sqrt(disc)) / (2*a);
    const t2 = (-b + Math.sqrt(disc)) / (2*a);
    if (Math.abs(t) > Math.abs(t2)) t = t2;

    /* Backtrack small particle (big is handled after all pairs) */
    px[i] += t * vx[i];
    py[i] += t * vy[i];

    /* Contact normal (from big centre to small centre) */
    const nx = (px[i] - bx), ny = (py[i] - by);
    const nd = Math.sqrt(nx*nx + ny*ny);
    const nxn = nx / nd, nyn = ny / nd;

    /* COM frame elastic collision */
    const totM  = ms + Mb;
    const comdx = (ms * vx[i] + Mb * bvx) / totM;
    const comdy = (ms * vy[i] + Mb * bvy) / totM;
    const pn    = (vx[i] - comdx) * nxn + (vy[i] - comdy) * nyn;
    const dpx   = 2 * nxn * pn, dpy = 2 * nyn * pn;
    vx[i]  -= dpx;        vy[i]  -= dpy;
    bvx    += dpx * ms/Mb; bvy   += dpy * ms/Mb;

    /* Advance forward */
    if (t < 0) { px[i] -= t * vx[i]; py[i] -= t * vy[i]; }
    clampSmall(i);
  }

  /* ── Exact collision between two small particles ── */
  function collideSmall(i, j) {
    const ddx = px[j] - px[i], ddy = py[j] - py[i];
    const d2  = ddx*ddx + ddy*ddy;
    const md  = 2 * rs;
    if (d2 >= md * md || d2 < 1e-12) return;

    const sdx = vx[i] - vx[j], sdy = vy[i] - vy[j];
    const sx  = px[i] - px[j], sy  = py[i] - py[j];
    const a   = sdx*sdx + sdy*sdy;
    if (a < 1e-14) return;
    const b   = 2 * (sx*sdx + sy*sdy);
    const c   = sx*sx + sy*sy - md*md;
    const disc = b*b - 4*a*c;
    if (disc < 0) return;
    let t = (-b - Math.sqrt(disc)) / (2*a);
    const t2 = (-b + Math.sqrt(disc)) / (2*a);
    if (Math.abs(t) > Math.abs(t2)) t = t2;

    px[i] += t * vx[i]; py[i] += t * vy[i];
    px[j] += t * vx[j]; py[j] += t * vy[j];

    const nx = px[i] - px[j], ny = py[i] - py[j];
    const nd = Math.sqrt(nx*nx + ny*ny);
    const nxn = nx / nd, nyn = ny / nd;

    /* Equal-mass COM: just exchange normal components */
    const dvn = (vx[i] - vx[j]) * nxn + (vy[i] - vy[j]) * nyn;
    vx[i] -= dvn * nxn; vy[i] -= dvn * nyn;
    vx[j] += dvn * nxn; vy[j] += dvn * nyn;

    if (t < 0) {
      px[i] -= t * vx[i]; py[i] -= t * vy[i];
      px[j] -= t * vx[j]; py[j] -= t * vy[j];
    }
    clampSmall(i); clampSmall(j);
  }

  function clampSmall(i) {
    if (px[i] - rs < 0)   { px[i] = rs;     vx[i] =  Math.abs(vx[i]); }
    if (px[i] + rs > L)   { px[i] = L - rs; vx[i] = -Math.abs(vx[i]); }
    if (py[i] - rs < 0)   { py[i] = rs;     vy[i] =  Math.abs(vy[i]); }
    if (py[i] + rs > L)   { py[i] = L - rs; vy[i] = -Math.abs(vy[i]); }
  }

  function clampBig() {
    if (bx - Rb < 0)   { bx = Rb;     bvx =  Math.abs(bvx); }
    if (bx + Rb > L)   { bx = L - Rb; bvx = -Math.abs(bvx); }
    if (by - Rb < 0)   { by = Rb;     bvy =  Math.abs(bvy); }
    if (by + Rb > L)   { by = L - Rb; bvy = -Math.abs(bvy); }
  }

  /* ── Step ── */
  function step(dt) {
    /* Gravity & integrate */
    for (let i = 0; i < N; i++) {
      vy[i] += gravity * dt;
      px[i] += vx[i] * dt;
      py[i] += vy[i] * dt;
    }
    bvy += gravity * dt;
    bx  += bvx * dt;
    by  += bvy * dt;

    /* Wall bounces */
    for (let i = 0; i < N; i++) clampSmall(i);
    clampBig();

    /* Build cell list and resolve collisions */
    buildCells();
    const cW = buildCells._cW, cH = buildCells._cH;

    for (let cy = 0; cy < nCY; cy++) {
      for (let cx = 0; cx < nCX; cx++) {
        const c = cy * nCX + cx;

        /* Neighbour pairs (half-shell) */
        for (let dcy = 0; dcy <= 1; dcy++) {
          for (let dcx = (dcy === 0 ? 1 : -1); dcx <= 1; dcx++) {
            const nx2 = cx + dcx, ny2 = cy + dcy;
            if (nx2 < 0 || nx2 >= nCX || ny2 < 0 || ny2 >= nCY) continue;
            const c2 = ny2 * nCX + nx2;
            for (let i = cellHead[c]; i !== -1; i = cellNext[i]) {
              for (let j = cellHead[c2]; j !== -1; j = cellNext[j]) {
                if (i === j) continue;
                if (i === N) { if (j < N) collideBig(j); }
                else if (j === N) collideBig(i);
                else collideSmall(i, j);
              }
            }
          }
        }

        /* Within-cell pairs */
        for (let i = cellHead[c]; i !== -1; i = cellNext[i]) {
          for (let j = cellNext[i]; j !== -1; j = cellNext[j]) {
            if (i === N) { if (j < N) collideBig(j); }
            else if (j === N) collideBig(i);
            else collideSmall(i, j);
          }
        }
      }
    }

    /* Record trail */
    trail.push(bx, by);
    if (trail.length > TRAIL_MAX * 2) trail.splice(0, 2);
  }

  /* ── Neon sprites (chain/three-body aesthetic): white-hot core through
     the light palette color, additive halo in the dark color ── */
  const N_SPRITES  = 32;
  const GLOW_SCALE = 2.6;
  let sprites = null, glowSprites = null, spriteRs = -1;
  let bigSprite = null, bigGlow = null, bigSpriteR = -1;
  const _siScratch = new Int32Array(600);   // matches bm-nslider max

  function buildOrb(rad, lr, lg, lb) {
    const c = document.createElement('canvas');
    c.width = c.height = Math.ceil(2 * rad) + 2;
    const sctx = c.getContext('2d');
    const cc = c.width / 2;
    const grad = sctx.createRadialGradient(cc, cc, 0, cc, cc, rad);
    grad.addColorStop(0.0,  `rgba(255,255,255,1.0)`);
    grad.addColorStop(0.35, `rgba(255,255,255,0.95)`);
    grad.addColorStop(0.65, `rgba(${lr},${lg},${lb},0.9)`);
    grad.addColorStop(1.0,  `rgba(${lr},${lg},${lb},0)`);
    sctx.beginPath(); sctx.arc(cc, cc, rad, 0, Math.PI * 2);
    sctx.fillStyle = grad; sctx.fill();
    return c;
  }

  function buildHalo(rad, dr, dg, db, peak) {
    const gR = rad * GLOW_SCALE;
    const c = document.createElement('canvas');
    c.width = c.height = Math.ceil(2 * gR) + 2;
    const gctx = c.getContext('2d');
    const cc = c.width / 2;
    const gg = gctx.createRadialGradient(cc, cc, 0, cc, cc, gR);
    gg.addColorStop(0.0, `rgba(${dr},${dg},${db},${peak})`);
    gg.addColorStop(0.4, `rgba(${dr},${dg},${db},${peak * 0.36})`);
    gg.addColorStop(1.0, `rgba(${dr},${dg},${db},0)`);
    gctx.beginPath(); gctx.arc(cc, cc, gR, 0, Math.PI * 2);
    gctx.fillStyle = gg; gctx.fill();
    return c;
  }

  function buildSprites() {
    sprites     = new Array(N_SPRITES);
    glowSprites = new Array(N_SPRITES);
    spriteRs = rs;
    for (let s = 0; s < N_SPRITES; s++) {
      const hot = s / (N_SPRITES - 1);
      const lr = Math.round(_TLR + (_PLR - _TLR) * hot);
      const lg = Math.round(_TLG + (_PLG - _TLG) * hot);
      const lb = Math.round(_TLB + (_PLB - _TLB) * hot);
      const dr = Math.round(_TDR + (_PDR - _TDR) * hot);
      const dg = Math.round(_TDG + (_PDG - _TDG) * hot);
      const db = Math.round(_TDB + (_PDB - _TDB) * hot);
      sprites[s]     = buildOrb(rs, lr, lg, lb);
      glowSprites[s] = buildHalo(rs, dr, dg, db, 0.55);
    }
  }

  /* ── Render ── */
  function render() {
    ctx.fillStyle = _c('--black');
    ctx.fillRect(0, 0, L, L);

    /* Trail — amber neon */
    if (trail.length >= 4) {
      ctx.save();
      ctx.shadowColor = _rgba('--amber', 0.8);
      ctx.shadowBlur  = 8;
      ctx.beginPath();
      ctx.moveTo(trail[0], trail[1]);
      for (let i = 2; i < trail.length; i += 2) ctx.lineTo(trail[i], trail[i+1]);
      ctx.strokeStyle = _rgba('--amber', 0.45);
      ctx.lineWidth = 1.5;
      ctx.stroke();
      ctx.restore();
    }

    /* Small particles — colour by speed, glow pass then body pass */
    let meanSpd = 0;
    for (let i = 0; i < N; i++) meanSpd += Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
    meanSpd /= Math.max(N, 1);
    const spdRef = Math.max(meanSpd * 2, 0.001);

    if (spriteRs !== rs) buildSprites();
    for (let i = 0; i < N; i++) {
      const spd = Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
      const hot = Math.min(spd / spdRef, 1);
      _siScratch[i] = (hot * (N_SPRITES - 1) + 0.5) | 0;
    }
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    const gHalf = glowSprites[0].width / 2;
    for (let i = 0; i < N; i++) {
      ctx.drawImage(glowSprites[_siScratch[i]], px[i] - gHalf, py[i] - gHalf);
    }
    ctx.restore();
    const bHalf = sprites[0].width / 2;
    for (let i = 0; i < N; i++) {
      ctx.drawImage(sprites[_siScratch[i]], px[i] - bHalf, py[i] - bHalf);
    }

    /* Big particle — white-hot orb with amber halo */
    if (bigSpriteR !== Rb) {
      bigSprite  = buildOrb(Rb, _AMR, _AMG, _AMB);
      bigGlow    = buildHalo(Rb, _AMR, _AMG, _AMB, 0.45);
      bigSpriteR = Rb;
    }
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    ctx.drawImage(bigGlow, bx - bigGlow.width / 2, by - bigGlow.width / 2);
    ctx.restore();
    ctx.drawImage(bigSprite, bx - bigSprite.width / 2, by - bigSprite.width / 2);
  }

  /* ── Histogram (Maxwell–Boltzmann of small particles) ── */
  const N_BINS = 30;
  let smoothBins = new Float64Array(N_BINS);
  let histYMax   = 0;
  const HIST_ALPHA = 0.08, YMAX_ALPHA = 0.03;

  function renderHistogram() {
    const hc = document.getElementById('bm-hist-canvas');
    const W  = hc.clientWidth  || 200;
    const H  = hc.clientHeight || 100;
    if (hc.width !== W || hc.height !== H) { hc.width = W; hc.height = H; }

    /* Compute actual mean KE to set axis range and fit MB curve */
    let meanKE = 0;
    for (let i = 0; i < N; i++) meanKE += vx[i]*vx[i] + vy[i]*vy[i];
    meanKE /= (2 * N);                 // <½mv²> per particle (m=1)
    const Teff = Math.max(meanKE, 1e-6);  // effective temperature from data
    const vref = Math.sqrt(Teff);
    const vmax = 6 * vref + 0.5;
    const dv   = vmax / N_BINS;
    const raw  = new Float64Array(N_BINS);
    for (let i = 0; i < N; i++) {
      const spd = Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
      const b   = Math.min((spd / dv) | 0, N_BINS - 1);
      raw[b]++;
    }
    let tot = 0;
    for (let b = 0; b < N_BINS; b++) tot += raw[b];
    if (tot > 0) for (let b = 0; b < N_BINS; b++) raw[b] /= tot;
    for (let b = 0; b < N_BINS; b++)
      smoothBins[b] += HIST_ALPHA * (raw[b] - smoothBins[b]);

    /* Maxwell–Boltzmann curve fitted to Teff (2D: P(v) ∝ v·exp(-v²/2Teff)) */
    const mb = new Float64Array(N_BINS);
    let mbSum = 0;
    for (let b = 0; b < N_BINS; b++) {
      const v = (b + 0.5) * dv;
      mb[b]  = (v / Teff) * Math.exp(-v * v / (2 * Teff));
      mbSum += mb[b];
    }
    if (mbSum > 0) for (let b = 0; b < N_BINS; b++) mb[b] /= mbSum;
    const mbPeak = Math.max(...mb, 1e-9);
    histYMax += YMAX_ALPHA * (mbPeak * 1.2 - histYMax);
    const ymax = Math.max(histYMax, 1e-9);

    const PL = 4, PR = 4, PT = 6, PB = 4;
    const pw = W - PL - PR, ph = H - PT - PB;
    const bw = pw / N_BINS;
    hctx.clearRect(0, 0, W, H);
    /* Bars: light-palette gradient (matches particle bodies) with a soft
       dark-palette glow underneath (matches their halos) */
    hctx.save();
    hctx.shadowBlur = 6;
    for (let b = 0; b < N_BINS; b++) {
      const bh = Math.min(smoothBins[b] / ymax, 1) * ph;
      const v   = (b + 0.5) * dv;
      const hot = Math.min(v / (Math.sqrt(Teff) * 2), 1);
      const lr  = Math.round(_TLR + (_PLR - _TLR) * hot);
      const lg  = Math.round(_TLG + (_PLG - _TLG) * hot);
      const lb  = Math.round(_TLB + (_PLB - _TLB) * hot);
      const dr  = Math.round(_TDR + (_PDR - _TDR) * hot);
      const dg  = Math.round(_TDG + (_PDG - _TDG) * hot);
      const db  = Math.round(_TDB + (_PDB - _TDB) * hot);
      hctx.shadowColor = `rgb(${dr},${dg},${db})`;
      hctx.fillStyle   = `rgba(${lr},${lg},${lb},0.75)`;
      hctx.fillRect(PL + b * bw, PT + ph - bh, bw - 1, bh);
    }
    hctx.restore();
    /* MB curve: neon stroke — dark glow under light line, like the rods */
    hctx.save();
    hctx.shadowColor = _c('--pink-dark');
    hctx.shadowBlur  = 10;
    hctx.strokeStyle = _c('--pink-light');
    hctx.lineWidth   = 1.5;
    hctx.beginPath();
    for (let b = 0; b < N_BINS; b++) {
      const x = PL + (b + 0.5) * bw;
      const y = PT + ph - Math.min(mb[b] / ymax, 1) * ph;
      b === 0 ? hctx.moveTo(x, y) : hctx.lineTo(x, y);
    }
    hctx.stroke();
    hctx.restore();
  }

  /* ── Loop ── */
  const DT_BASE = 0.8, SUBSTEPS = 8;
  const dt = DT_BASE / SUBSTEPS;

  function loop() {
    if (running) {
      for (let s = 0; s < SUBSTEPS; s++) step(dt);
      render();
      renderHistogram();
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Reset ── */
  function doReset() {
    N       = parseInt(document.getElementById('bm-nslider').value);
    T       = parseFloat(document.getElementById('bm-tslider').value);
    rs      = parseInt(document.getElementById('bm-rslider').value);
    gravity = 0;
    document.getElementById('bm-gslider').value = 0;
    smoothBins = new Float64Array(N_BINS);
    histYMax   = 0;
    initSim();
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'bm',
    title: 'Brownian Motion &mdash; Hard Spheres',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="bmReset()">Reset</button><button class="applet-shell-header-btn" id="bm-pause-btn" onclick="bmTogglePause()">Pause</button>`,


    ctrlHTML: `
      <div id="bm-sliders">
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Particles</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Few</span>
            <input type="range" id="bm-nslider" min="50" max="600" step="10" value="500">
            <span class="applet-shell-side">Many</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Temperature</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Cold</span>
            <input type="range" id="bm-tslider" min="0.1" max="5.0" step="0.05" value="1.0">
            <span class="applet-shell-side">Hot</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Small Radius</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Small</span>
            <input type="range" id="bm-rslider" min="2" max="8" step="1" value="5">
            <span class="applet-shell-side">Large</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Gravity</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Off</span>
            <input type="range" id="bm-gslider" min="0" max="1" step="0.01" value="0">
            <span class="applet-shell-side">Strong</span>
          </div>
        </div>
      </div>
      <div id="bm-hist-section">
        <div class="applet-shell-ctrl-title">Speed distribution</div>
        <canvas id="bm-hist-canvas"></canvas>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      const hc = document.getElementById('bm-hist-canvas');
      hctx = hc.getContext('2d');
      L = S;
      setTimeout(() => {
        hc.width  = hc.clientWidth  || 200;
        hc.height = hc.clientHeight || 120;
      }, 80);
      doReset();
      running = true;
      const pb = document.getElementById('bm-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) frameId = requestAnimationFrame(loop);
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('bm-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ S }) {
      L = S;
      doReset();
    },
  });

  window.bmOpen  = () => shell.open();
  window.bmClose = () => shell.close();
  window.bmReset = doReset;

  window.bmTogglePause = function () {
    running = !running;
    const pb = document.getElementById('bm-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  document.getElementById('bm-rslider').addEventListener('input', function () {
    const newRs = parseInt(this.value);
    if (newRs === rs) return;
    const msOld = ms;
    rs = newRs;
    ms = rs * rs;
    /* rescale speeds so per-particle KE (temperature) is preserved */
    if (vx && msOld > 0) {
      const sc = Math.sqrt(msOld / ms);
      for (let i = 0; i < N; i++) { vx[i] *= sc; vy[i] *= sc; }
    }
  });

  document.getElementById('bm-tslider').addEventListener('input', function () {
    T = parseFloat(this.value);
    if (!vx) return;
    let ke = 0;
    for (let i = 0; i < N; i++) ke += vx[i]*vx[i] + vy[i]*vy[i];
    ke /= (2 * N);
    const sc = Math.sqrt(T / Math.max(ke, 1e-12));
    for (let i = 0; i < N; i++) { vx[i] *= sc; vy[i] *= sc; }
    smoothBins = new Float64Array(N_BINS); histYMax = 0;
  });

  document.getElementById('bm-gslider').addEventListener('input', function () {
    gravity = parseFloat(this.value) * 0.10;
  });

})();

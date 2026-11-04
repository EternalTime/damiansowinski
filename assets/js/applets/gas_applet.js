(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  /* ── Inject CSS ── */
  (function () {
    if (document.getElementById('gas-styles')) return;
    const s = document.createElement('style');
    s.id = 'gas-styles';
    s.textContent = `
      #gas-ctrl-panel, #gas-asm-ctrl {
        display: flex;
        flex-direction: column;
        overflow: hidden;
      }
      #gas-scrollable {
        flex: 1;
        overflow-y: auto;
        min-height: 0;
      }
      #gas-hist-section {
        flex: 0 0 33%;
        min-height: 0;
        display: flex;
        flex-direction: column;
        padding: 6px 12px 8px;
        border-top: 1px solid var(--border-dark);
      }
      #gas-hist-section .applet-shell-ctrl-title {
        flex-shrink: 0;
        margin-bottom: 4px;
      }
      #gas-hist-canvas {
        flex: 1;
        min-height: 0;
        display: block;
        width: 100%;
      }
    `;
    document.head.appendChild(s);
  })();

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');
  const [_PLR, _PLG, _PLB] = _rgb('--pink-light');
  const [_TDR, _TDG, _TDB] = _rgb('--teal-dark');

  /* ── Simulation state ── */
  const N = 300;
  let T  = 1.0;
  let R  = 10;
  let L  = 500;

  let px, py, vx, vy;
  let running = false, frameId = null;
  let wasRunning = false;   // sim state stashed while the docs panel is open
  let canvas, ctx, hctx;

  let wallX      = 0;
  let pistonIn   = false;
  let jammed     = false;
  const WALL_SPEED = 0.06;   // per base substep → wall velocity ≈ 0.5, subsonic at T=1 (quasi-static)
  let gravity = 0;

  /* The piston advances while the gas can make room for it and waits when it cannot. Disks
     closer than one base piston step count as in contact, and the gas is jammed when contacts
     chain from the piston face to the far wall: the next step would have to push that whole
     chain into overlap. */
  const CONTACT = WALL_SPEED;
  const _root = new Int32Array(N), _atPiston = new Uint8Array(N);
  function _find(i) {
    while (_root[i] !== i) { _root[i] = _root[_root[i]]; i = _root[i]; }
    return i;
  }

  let cellW, nCellsX, nCellsY;
  let cellHead, cellNext;

  function buildCellList() {
    const boxW = L - wallX;
    cellW      = Math.max(2 * R + CONTACT, 1);
    nCellsX    = Math.max(Math.floor(boxW / cellW), 1);
    nCellsY    = Math.max(Math.floor(L    / cellW), 1);
    const cW   = boxW / nCellsX;
    const cH   = L    / nCellsY;
    const nCells = nCellsX * nCellsY;
    if (!cellHead || cellHead.length < nCells) cellHead = new Int32Array(nCells);
    cellHead.fill(-1, 0, nCells);
    if (!cellNext) cellNext = new Int32Array(N);
    for (let i = 0; i < N; i++) {
      const cx = Math.min(Math.max(((px[i] - wallX) / cW) | 0, 0), nCellsX - 1);
      const cy = Math.min(Math.max((py[i]           / cH) | 0, 0), nCellsY - 1);
      const c  = cy * nCellsX + cx;
      cellNext[i] = cellHead[c];
      cellHead[c] = i;
    }
    buildCellList._cW = cW;
    buildCellList._cH = cH;
  }

  function gaussRand() {
    let u;
    do { u = Math.random(); } while (u === 0);
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * Math.random());
  }

  function initParticles() {
    px = new Float64Array(N); py = new Float64Array(N);
    vx = new Float64Array(N); vy = new Float64Array(N);
    const leftEdge = wallX + R + 1;
    const marginY  = R + 1;
    const innerW   = Math.max(L - leftEdge - R - 1, 1);
    const innerH   = Math.max(L - 2 * marginY, 1);
    const cols = Math.max(Math.ceil(Math.sqrt(N * innerW / Math.max(innerH, 1))), 1);
    const rows = Math.ceil(N / cols);
    const dx = innerW / cols, dy = innerH / rows;
    /* Every disk is placed. The square grid is used where it holds them apart; where it
       cannot, the widest triangular lattice that takes them all, the densest packing there
       is; and where even that cannot hold them apart they start overlapping, and the
       collisions push them apart as far as the box allows. */
    let placed = 0;
    if (Math.min(dx, dy) >= 2 * R) {
      for (let r = 0; r < rows && placed < N; r++) {
        for (let c = 0; c < cols && placed < N; c++) {
          px[placed] = leftEdge + (c + 0.5) * dx;
          py[placed] = marginY  + (r + 0.5) * dy;
          placed++;
        }
      }
    } else {
      const fit = a => {
        const rowsT = Math.floor(innerH / (a * Math.sqrt(3) / 2)) + 1;
        const even  = Math.floor(innerW / a) + 1;
        const odd   = innerW >= a / 2 ? Math.floor((innerW - a / 2) / a) + 1 : 0;
        return { rowsT, even, odd, count: Math.ceil(rowsT / 2) * even + Math.floor(rowsT / 2) * odd };
      };
      let lo = 0, hi = Math.max(innerW, innerH) + 1;
      for (let k = 0; k < 60; k++) { const mid = (lo + hi) / 2; if (fit(mid).count >= N) lo = mid; else hi = mid; }
      const a = lo, { rowsT, even, odd } = fit(a), h = a * Math.sqrt(3) / 2;
      for (let r = 0; r < rowsT && placed < N; r++) {
        const shift = r % 2 ? a / 2 : 0;
        for (let c = 0; c < (r % 2 ? odd : even) && placed < N; c++) {
          px[placed] = leftEdge + shift + c * a;
          py[placed] = marginY  + r * h;
          placed++;
        }
      }
      /* centre the lattice in the box, so its slack is shared by the four walls */
      let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
      for (let i = 0; i < N; i++) {
        x0 = Math.min(x0, px[i]); x1 = Math.max(x1, px[i]);
        y0 = Math.min(y0, py[i]); y1 = Math.max(y1, py[i]);
      }
      const sx = leftEdge + (innerW - (x1 - x0)) / 2 - x0, sy = marginY + (innerH - (y1 - y0)) / 2 - y0;
      for (let i = 0; i < N; i++) { px[i] += sx; py[i] += sy; }
    }
    for (let i = 0; i < N; i++) {
      vx[i] = gaussRand() * Math.sqrt(T);
      vy[i] = gaussRand() * Math.sqrt(T);
    }
    let svx = 0, svy = 0;
    for (let i = 0; i < N; i++) { svx += vx[i]; svy += vy[i]; }
    svx /= N; svy /= N;
    for (let i = 0; i < N; i++) { vx[i] -= svx; vy[i] -= svy; }
    let ke = 0;
    for (let i = 0; i < N; i++) ke += vx[i]*vx[i] + vy[i]*vy[i];
    ke /= (2 * N);
    const sc = Math.sqrt(T / Math.max(ke, 1e-12));
    for (let i = 0; i < N; i++) { vx[i] *= sc; vy[i] *= sc; }
    histT = T;
  }

  function collidePair(i, j, d2min, d2touch) {
    const ddx = px[j] - px[i], ddy = py[j] - py[i];
    const d2  = ddx*ddx + ddy*ddy;
    if (d2 < d2touch) { const a = _find(i), b = _find(j); if (a !== b) _root[a] = b; }
    if (d2 < d2min && d2 > 1e-12) {
      const d  = Math.sqrt(d2);
      const nx = ddx / d, ny = ddy / d;
      const ov = 2 * R - d;
      px[i] -= nx * ov * 0.5; py[i] -= ny * ov * 0.5;
      px[j] += nx * ov * 0.5; py[j] += ny * ov * 0.5;
      const dvx = vx[j] - vx[i], dvy = vy[j] - vy[i];
      const vn  = dvx * nx + dvy * ny;
      if (vn < 0) {
        vx[i] += vn * nx; vy[i] += vn * ny;
        vx[j] -= vn * nx; vy[j] -= vn * ny;
      }
    }
  }

  function step(dt) {
    if (gravity > 0) for (let i = 0; i < N; i++) vy[i] += gravity * dt;
    for (let i = 0; i < N; i++) { px[i] += vx[i] * dt; py[i] += vy[i] * dt; }
    /* Pushed in, the piston takes one step unless the gas is jammed, and never leaves less
       than a disk's width of box; pulled out, it withdraws to the left wall. */
    const wallFrom = wallX;
    if (pistonIn && !jammed) wallX = Math.max(Math.min(wallX + WALL_VELOCITY * dt, L - 2 * R - 2), wallFrom);
    else if (!pistonIn)      wallX = Math.max(wallX - WALL_VELOCITY * dt, 0);
    const v_wall = (wallX - wallFrom) / dt;
    for (let i = 0; i < N; i++) {
      if (px[i] - R < wallX) {
        px[i] = wallX + R;
        if (vx[i] < v_wall) vx[i] = 2 * v_wall - vx[i];  // reflect only if approaching
      }
      if (px[i] + R > L) { px[i] = L - R; vx[i] = -Math.abs(vx[i]); }
      if (py[i] - R < 0) { py[i] = R;     vy[i] =  Math.abs(vy[i]); }
      if (py[i] + R > L) { py[i] = L - R; vy[i] = -Math.abs(vy[i]); }
    }
    buildCellList();
    const d2min   = (2 * R) * (2 * R);
    const d2touch = pistonIn ? (2 * R + CONTACT) * (2 * R + CONTACT) : 0;
    if (pistonIn) for (let i = 0; i < N; i++) _root[i] = i;
    for (let cy = 0; cy < nCellsY; cy++) {
      for (let cx = 0; cx < nCellsX; cx++) {
        const c = cy * nCellsX + cx;
        for (let dcy = 0; dcy <= 1; dcy++) {
          for (let dcx = (dcy === 0 ? 1 : -1); dcx <= 1; dcx++) {
            const nx2 = cx + dcx, ny2 = cy + dcy;
            if (nx2 < 0 || nx2 >= nCellsX || ny2 < 0 || ny2 >= nCellsY) continue;
            const c2 = ny2 * nCellsX + nx2;
            for (let i = cellHead[c];  i !== -1; i = cellNext[i])
              for (let j = cellHead[c2]; j !== -1; j = cellNext[j])
                collidePair(i, j, d2min, d2touch);
          }
        }
        for (let i = cellHead[c]; i !== -1; i = cellNext[i])
          for (let j = cellNext[i]; j !== -1; j = cellNext[j])
            collidePair(i, j, d2min, d2touch);
      }
    }
    jammed = false;
    if (pistonIn) {
      for (let i = 0; i < N; i++) _atPiston[i] = 0;
      for (let i = 0; i < N; i++) if (px[i] - R - wallX < CONTACT) _atPiston[_find(i)] = 1;
      for (let i = 0; i < N && !jammed; i++) if (L - px[i] - R < CONTACT && _atPiston[_find(i)]) jammed = true;
    }
  }

  /* Pre-rendered particle sprites, bucketed by temperature.
     Neon aesthetic: white-hot core through the light palette color,
     additive halo in the paired dark color. */
  const N_SPRITES  = 32;
  const GLOW_SCALE = 2.6;
  let sprites = null, glowSprites = null, spriteR = -1;
  const _siScratch = new Int32Array(N);

  function buildSprites() {
    sprites     = new Array(N_SPRITES);
    glowSprites = new Array(N_SPRITES);
    spriteR = R;
    for (let s = 0; s < N_SPRITES; s++) {
      const hot = s / (N_SPRITES - 1);
      const lr = Math.round(_TLR + (_PLR - _TLR) * hot);
      const lg = Math.round(_TLG + (_PLG - _TLG) * hot);
      const lb = Math.round(_TLB + (_PLB - _TLB) * hot);
      const dr = Math.round(_TDR + (_PDR - _TDR) * hot);
      const dg = Math.round(_TDG + (_PDG - _TDG) * hot);
      const db = Math.round(_TDB + (_PDB - _TDB) * hot);
      /* body: white core → light color → transparent edge */
      const c = document.createElement('canvas');
      c.width = c.height = 2 * R + 2;
      const sctx = c.getContext('2d');
      const cx = R + 1, cy = R + 1;
      const grad = sctx.createRadialGradient(cx, cy, 0, cx, cy, R);
      grad.addColorStop(0.00, `rgba(255,255,255,1.0)`);
      grad.addColorStop(0.35, `rgba(255,255,255,0.95)`);
      grad.addColorStop(0.65, `rgba(${lr},${lg},${lb},0.9)`);
      grad.addColorStop(1.00, `rgba(${lr},${lg},${lb},0)`);
      sctx.beginPath();
      sctx.arc(cx, cy, R, 0, Math.PI * 2);
      sctx.fillStyle = grad;
      sctx.fill();
      sprites[s] = c;
      /* halo */
      const gR = R * GLOW_SCALE;
      const gc = document.createElement('canvas');
      gc.width = gc.height = Math.ceil(2 * gR) + 2;
      const gctx = gc.getContext('2d');
      const gcx = gc.width / 2, gcy = gc.height / 2;
      const gg = gctx.createRadialGradient(gcx, gcy, 0, gcx, gcy, gR);
      gg.addColorStop(0.0, `rgba(${dr},${dg},${db},0.55)`);
      gg.addColorStop(0.4, `rgba(${dr},${dg},${db},0.20)`);
      gg.addColorStop(1.0, `rgba(${dr},${dg},${db},0)`);
      gctx.beginPath();
      gctx.arc(gcx, gcy, gR, 0, Math.PI * 2);
      gctx.fillStyle = gg;
      gctx.fill();
      glowSprites[s] = gc;
    }
  }

  function render() {
    ctx.fillStyle = _c('--black');
    ctx.fillRect(0, 0, L, L);
    if (wallX > 0) {
      ctx.fillStyle = _rgba('--teal-light', 0.25);
      ctx.fillRect(0, 0, wallX, L);
      ctx.save();
      ctx.shadowColor = _c('--teal-light');
      ctx.shadowBlur  = 14;
      ctx.fillStyle = _c('--teal-light');
      ctx.fillRect(wallX - 2, 0, 3, L);
      ctx.restore();
    }
    let meanSpd = 0;
    for (let i = 0; i < N; i++) meanSpd += Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
    meanSpd /= N;
    const spdRef = Math.max(meanSpd * 2, 0.001);
    if (spriteR !== R) buildSprites();
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
    for (let i = 0; i < N; i++) {
      ctx.drawImage(sprites[_siScratch[i]], px[i] - R - 1, py[i] - R - 1);
    }
  }

  const N_BINS = 30;
  let smoothBins = new Float64Array(N_BINS);
  let histYMax   = 0;
  let histT      = 1;   // the gas's own temperature, smoothed, which sets the speed axis
  const HIST_ALPHA = 0.10, YMAX_ALPHA = 0.04;

  function renderHistogram() {
    const hc = document.getElementById('gas-hist-canvas');
    const W  = hc.clientWidth  || 200;
    const H  = hc.clientHeight || 120;
    if (hc.width !== W || hc.height !== H) { hc.width = W; hc.height = H; }
    /* The speed axis follows the temperature of the gas being drawn, smoothed as the height
       is, so the bars stay framed however the piston or gravity has heated or cooled it;
       the curve stays at the slider temperature, so a gas hotter than the slider shows as
       bars lying to the right of it. */
    let ke = 0;
    for (let i = 0; i < N; i++) ke += vx[i]*vx[i] + vy[i]*vy[i];
    histT += YMAX_ALPHA * (ke / (2 * N) - histT);
    const vref = Math.sqrt(histT);
    const vmax = 6 * vref + 0.5;
    const dv   = vmax / N_BINS;
    const raw  = new Float64Array(N_BINS);
    for (let i = 0; i < N; i++) {
      const spd = Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
      const b   = (spd / dv) | 0;
      if (b < N_BINS) raw[b]++;   // out-of-range values fall in unprinted tail bins
    }
    let tot = 0;
    for (let b = 0; b < N_BINS; b++) tot += raw[b];
    if (tot > 0) for (let b = 0; b < N_BINS; b++) raw[b] /= tot;
    for (let b = 0; b < N_BINS; b++)
      smoothBins[b] += HIST_ALPHA * (raw[b] - smoothBins[b]);
    /* the curve at the slider temperature, and the height set by the same law at the
       gas's own temperature, which is what the bars follow */
    const mbAt = v => (v / T) * Math.exp(-v * v / (2 * T));
    const gas = new Float64Array(N_BINS);
    let mbSum = 0, gasSum = 0;
    for (let b = 0; b < N_BINS; b++) {
      const v = (b + 0.5) * dv;
      gas[b]  = (v / histT) * Math.exp(-v * v / (2 * histT));
      mbSum  += mbAt(v);
      gasSum += gas[b];
    }
    const gasPeak = gasSum > 0 ? Math.max(...gas) / gasSum : 1e-9;
    histYMax += YMAX_ALPHA * (gasPeak * 1.20 - histYMax);
    const ymax = Math.max(histYMax, 1e-9);
    const PL = 4, PR = 4, PT = 6, PB = 4;
    const pw = W - PL - PR, ph = H - PT - PB;
    const bw = pw / N_BINS;
    hctx.clearRect(0, 0, W, H);
    /* Bars: light-palette gradient with dark-palette glow (matches particles) */
    hctx.save();
    hctx.shadowBlur = 6;
    for (let b = 0; b < N_BINS; b++) {
      const bh  = Math.min(smoothBins[b] / ymax, 1) * ph;
      const v   = (b + 0.5) * dv;
      const hot = Math.min(v / (Math.sqrt(T) * 2), 1);
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
    hctx.save();
    hctx.shadowColor = _c('--pink-dark');
    hctx.shadowBlur  = 10;
    hctx.strokeStyle = _c('--pink-light');
    hctx.lineWidth   = 1.5;
    /* The curve keeps the bars' scale, the share of the gas a bin centred on each speed
       would hold, but is sampled once per device pixel across the panel so it reads as a
       curve rather than a run of segments between bin centres. */
    const steps = Math.ceil(pw * (window.devicePixelRatio || 1));
    hctx.beginPath();
    for (let k = 0; k <= steps; k++) {
      const x = PL + pw * k / steps;
      const p = mbSum > 0 ? mbAt(vmax * k / steps) / mbSum : 0;
      const y = PT + ph - Math.min(p / ymax, 1) * ph;
      k === 0 ? hctx.moveTo(x, y) : hctx.lineTo(x, y);
    }
    hctx.stroke();
    hctx.restore();
  }

  /* Each frame advances the same time in as many substeps as keep the fastest disk to a
     twentieth of its radius per substep, so a gas the piston has heated is followed as
     closely as a cold one: four substeps at the least, sixty-four at the most. */
  const DT_BASE = 0.5, SUBSTEPS = 4, MAX_SUBSTEPS = 64, STEP_REACH = 0.05;
  const WALL_VELOCITY = WALL_SPEED * SUBSTEPS / DT_BASE;

  function loop() {
    if (running) {
      let v2max = 0;
      for (let i = 0; i < N; i++) v2max = Math.max(v2max, vx[i]*vx[i] + vy[i]*vy[i]);
      const n = Math.min(Math.max(Math.ceil(Math.sqrt(v2max) * DT_BASE / (STEP_REACH * R)), SUBSTEPS), MAX_SUBSTEPS);
      for (let s = 0; s < n; s++) step(DT_BASE / n);
      render();
      renderHistogram();
    }
    frameId = requestAnimationFrame(loop);
  }

  function withdrawPiston() {
    wallX      = 0;
    pistonIn   = false;
    jammed     = false;
    const pb = document.getElementById('gas-piston-btn');
    if (pb) { pb.textContent = 'Push In'; pb.classList.remove('active'); }
  }

  function doReset() {
    T          = parseFloat(document.getElementById('gas-tslider').value);
    R          = parseInt(document.getElementById('gas-rslider').value);
    gravity    = 0;
    withdrawPiston();
    document.getElementById('gas-gslider').value = 0;
    smoothBins = new Float64Array(N_BINS);
    histYMax   = 0;
    initParticles();
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'gas',
    title: 'Hard-Sphere Gas',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="gasReset()">Reset</button><button class="applet-shell-header-btn" id="gas-pause-btn" onclick="gasTogglePause()">Pause</button>`,


    docs: {
      whatis: `In 1738 Daniel Bernoulli proposed, in his Hydrodynamica, that a gas is a swarm of tiny bodies in ceaseless motion, and that pressure is nothing but their drumming on the container walls [bernoulli1738]. The idea lay almost forgotten for most of a century - heat was supposed to be a fluid called caloric, not motion - until Herapath, Joule, Krönig and Clausius revived it; James Clerk Maxwell then derived the equilibrium velocity distribution of such a swarm in 1860 [maxwell1860], and Ludwig Boltzmann showed in 1872 that collisions drive any initial distribution toward it, provided the velocities of colliding molecules are uncorrelated [boltzmann1872]. In two dimensions the speed distribution reads
$$P(v) = \\frac{m v}{k_B T} \\, e^{-m v^2 / 2 k_B T},$$
the two-dimensional form of the law Maxwell found in three.¶The hard-sphere gas strips the picture to its bones: disks that fly freely and collide elastically, nothing else. No attractions, no internal structure, no adjustable potential. Remarkably, this bare minimum suffices for equilibration, for pressure, for the Maxwell-Boltzmann distribution; as Alder and Wainwright discovered in 1957, in one of the first molecular dynamics simulations ever run, it even suffices for a freezing transition at high density [alder1957], which in 1962 they found for disks too [alder1962]. Entropy alone can crystallize a gas.¶The walls make the thermodynamics tangible. Push the piston inward and the moving wall returns each particle slightly faster than it arrived: microscopic work heating the gas, compression made mechanical. Switch on gravity and the box becomes a miniature atmosphere, its density thinning with height, barometrically where the gas is dilute.`,

      howto: `Disks are colored by speed, teal for slow through pink for fast. The teal wall on the left is the piston face; the panel below histograms the speeds against the two-dimensional Maxwell-Boltzmann distribution (pink curve) at the slider temperature.¶Temperature rescales all velocities immediately. Push In drives the piston in until the disks jam and can make no more room, and the histogram drifts hot as the gas is compressed; Pull Out withdraws it. Gravity pulls the disks down into a profile that thins with height.¶Radius takes effect on Reset, which re-grids the gas with fresh Maxwell-Boltzmann velocities and withdraws the piston. Crank Radius high, push the piston in, and look for crystalline patches: Alder and Wainwright's transition, live. Pause freezes the dynamics.`,

      references: ['bernoulli1738', 'maxwell1860', 'boltzmann1872', 'alder1957', 'alder1962'],
    },

    onDocsOpen:  function () { wasRunning = running; running = false; },
    onDocsClose: function () { running = wasRunning; },

    ctrlHTML: `
      <div id="gas-scrollable">
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Temperature</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Cold</span>
            <input type="range" id="gas-tslider" min="0.1" max="5.0" step="0.05" value="1.0">
            <span class="applet-shell-side">Hot</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Radius</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Small</span>
            <input type="range" id="gas-rslider" min="5" max="12" step="1" value="10">
            <span class="applet-shell-side">Large</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Piston</div>
          <div class="applet-shell-btn-row">
            <button class="applet-shell-btn" id="gas-piston-btn" onclick="gasTogglePiston()">Push In</button>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Gravity</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Off</span>
            <input type="range" id="gas-gslider" min="0" max="1" step="0.01" value="0">
            <span class="applet-shell-side">Strong</span>
          </div>
        </div>
      </div>
      <div id="gas-hist-section">
        <div class="applet-shell-ctrl-title">Speed distribution</div>
        <canvas id="gas-hist-canvas"></canvas>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      hctx   = document.getElementById('gas-hist-canvas').getContext('2d');
      L = S;
      doReset();
      running = true;
      const pb = document.getElementById('gas-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) frameId = requestAnimationFrame(loop);
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('gas-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ S }) {
      L = S;
      withdrawPiston();
      initParticles();
    },
  });

  window.gasOpen  = () => shell.open();
  window.gasClose = () => shell.close();
  window.gasReset = doReset;
  window.gasTogglePause = function () {
    running = !running;
    const pb = document.getElementById('gas-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  document.getElementById('gas-tslider').addEventListener('input', function () {
    T = parseFloat(this.value);
    if (!px) return;
    let ke = 0;
    for (let i = 0; i < N; i++) ke += vx[i]*vx[i] + vy[i]*vy[i];
    ke /= (2 * N);
    const sc = Math.sqrt(T / Math.max(ke, 1e-12));
    for (let i = 0; i < N; i++) { vx[i] *= sc; vy[i] *= sc; }
    histT      = T;
    smoothBins = new Float64Array(N_BINS);
    histYMax   = 0;
  });

  document.getElementById('gas-gslider').addEventListener('input', function () {
    gravity = parseFloat(this.value) * 0.12;
  });

  /* Push In drives the piston until the gas jams; Pull Out withdraws it. */
  window.gasTogglePiston = function () {
    pistonIn   = !pistonIn;
    jammed     = false;
    const pb = document.getElementById('gas-piston-btn');
    if (pb) {
      pb.textContent = pistonIn ? 'Pull Out' : 'Push In';
      pb.classList.toggle('active', pistonIn);
    }
  };

})();

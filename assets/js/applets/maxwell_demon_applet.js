(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');

  /* ── Simulation parameters ── */
  const N_EACH   = 50;
  let R          = 8;
  const T_INIT   = 1.0;
  const DT_BASE  = 0.5;
  const SUBSTEPS = 4;
  const dt       = DT_BASE / SUBSTEPS;

  /* ── Door ── */
  let doorOpen = false;
  let doorFrac = 0;
  const DOOR_SPEED = 0.04;

  /* ── Canvas / box geometry ── */
  let canvas, ctx;
  let BW, BH, midX;

  /* ── Particle arrays ── */
  let N, px, py, vx, vy, side;
  let T0 = 1.0;

  /* ── Running state ── */
  let running = false, frameId = null;

  /* ── Gaussian RNG ── */
  function gaussRand() {
    let u;
    do { u = Math.random(); } while (u === 0);
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * Math.random());
  }

  /* ── Initialise particles ── */
  function initParticles() {
    N    = 2 * N_EACH;
    px   = new Float64Array(N);
    py   = new Float64Array(N);
    vx   = new Float64Array(N);
    vy   = new Float64Array(N);
    side = new Int8Array(N);

    for (let side_idx = 0; side_idx < 2; side_idx++) {
      const xMin = side_idx === 0 ? R + 1        : midX + R + 5;
      const xMax = side_idx === 0 ? midX - R - 5 : BW - R - 1;
      const cols = Math.max(Math.ceil(Math.sqrt(N_EACH * (xMax - xMin) / Math.max(BH - 2*(R+1), 1))), 1);
      const rows = Math.ceil(N_EACH / cols);
      const dx   = (xMax - xMin) / cols;
      const dy   = (BH - 2*(R+1)) / rows;
      let placed = 0;
      const base = side_idx * N_EACH;
      outer:
      for (let r = 0; r < rows && placed < N_EACH; r++) {
        for (let c = 0; c < cols && placed < N_EACH; c++) {
          const x0 = xMin + (c + 0.5) * dx;
          const y0 = (R+1) + (r + 0.5) * dy;
          for (let j = 0; j < placed; j++) {
            const ddx = x0 - px[base+j], ddy = y0 - py[base+j];
            if (ddx*ddx + ddy*ddy < (2*R)*(2*R)) continue outer;
          }
          px[base+placed]   = x0;
          py[base+placed]   = y0;
          vx[base+placed]   = gaussRand() * Math.sqrt(T_INIT);
          vy[base+placed]   = gaussRand() * Math.sqrt(T_INIT);
          side[base+placed] = side_idx;
          placed++;
        }
      }
    }

    for (let side_idx = 0; side_idx < 2; side_idx++) {
      const base = side_idx * N_EACH;
      let svx = 0, svy = 0;
      for (let i = base; i < base + N_EACH; i++) { svx += vx[i]; svy += vy[i]; }
      svx /= N_EACH; svy /= N_EACH;
      for (let i = base; i < base + N_EACH; i++) { vx[i] -= svx; vy[i] -= svy; }
      let ke = 0;
      for (let i = base; i < base + N_EACH; i++) ke += vx[i]*vx[i] + vy[i]*vy[i];
      ke /= (2 * N_EACH);
      const sc = Math.sqrt(T_INIT / Math.max(ke, 1e-12));
      for (let i = base; i < base + N_EACH; i++) { vx[i] *= sc; vy[i] *= sc; }
    }
    T0 = T_INIT;
  }

  /* ── Physics step ── */
  function step(dt) {
    for (let i = 0; i < N; i++) { px[i] += vx[i] * dt; py[i] += vy[i] * dt; }

    const gapHalf    = doorFrac * 3 * R;
    const doorCY     = BH * 0.5;
    const gapTop     = doorCY - gapHalf;
    const gapBot     = doorCY + gapHalf;
    const wallExclude = R * 1.5;

    for (let i = 0; i < N; i++) {
      if (px[i] - R < 0   && vx[i] < 0) { px[i] = R;      vx[i] = -vx[i]; }
      if (px[i] + R > BW  && vx[i] > 0) { px[i] = BW - R; vx[i] = -vx[i]; }
      if (py[i] - R < 0   && vy[i] < 0) { py[i] = R;      vy[i] = -vy[i]; }
      if (py[i] + R > BH  && vy[i] > 0) { py[i] = BH - R; vy[i] = -vy[i]; }

      const inGap = py[i] > gapTop && py[i] < gapBot;
      if (!inGap) {
        if (side[i] === 0 && px[i] > midX - wallExclude) {
          px[i] = midX - wallExclude; if (vx[i] > 0) vx[i] = -vx[i];
        } else if (side[i] === 1 && px[i] < midX + wallExclude) {
          px[i] = midX + wallExclude; if (vx[i] < 0) vx[i] = -vx[i];
        }
      } else {
        if (side[i] === 0 && px[i] > midX) side[i] = 1;
        else if (side[i] === 1 && px[i] < midX) side[i] = 0;
      }
    }

    const d2min = (2 * R) * (2 * R);
    for (let i = 0; i < N - 1; i++) {
      for (let j = i + 1; j < N; j++) {
        const ddx = px[j] - px[i], ddy = py[j] - py[i];
        const d2  = ddx*ddx + ddy*ddy;
        if (d2 < d2min && d2 > 1e-12) {
          const d = Math.sqrt(d2);
          const nx = ddx/d, ny = ddy/d;
          const dvx = vx[j]-vx[i], dvy = vy[j]-vy[i];
          const vn  = dvx*nx + dvy*ny;
          if (vn < 0) {
            vx[i] += vn*nx; vy[i] += vn*ny;
            vx[j] -= vn*nx; vy[j] -= vn*ny;
            const ov = 2*R - d;
            px[i] -= nx*ov*0.5; py[i] -= ny*ov*0.5;
            px[j] += nx*ov*0.5; py[j] += ny*ov*0.5;
            const we = wallExclude;
            if (side[i]===0 && px[i]>midX-we) px[i]=midX-we;
            if (side[i]===1 && px[i]<midX+we) px[i]=midX+we;
            if (side[j]===0 && px[j]>midX-we) px[j]=midX-we;
            if (side[j]===1 && px[j]<midX+we) px[j]=midX+we;
          }
        }
      }
    }
  }

  /* ── Temperature calculation ── */
  function sideTemp(s) {
    const gapHalf = doorFrac * 3 * R;
    const gapTop  = BH*0.5 - gapHalf, gapBot = BH*0.5 + gapHalf;
    let ke = 0, count = 0;
    for (let i = 0; i < N; i++) {
      if (side[i] !== s) continue;
      if (py[i] > gapTop && py[i] < gapBot) continue;
      ke += vx[i]*vx[i] + vy[i]*vy[i]; count++;
    }
    return count > 0 ? ke / (2*count) : 0;
  }

  /* ── Render ── */
  function render() {
    ctx.fillStyle = _c('--bg-deep');
    ctx.fillRect(0, 0, BW, BH);

    const gapHalf = doorFrac * 3 * R;
    const doorCY  = BH * 0.5;
    const gapTop  = doorCY - gapHalf;
    const gapBot  = doorCY + gapHalf;
    const wallW   = R;

    ctx.fillStyle = _c('--teal-light');
    if (gapTop > 0)  ctx.fillRect(midX - wallW*0.5, 0,      wallW, gapTop);
    if (gapBot < BH) ctx.fillRect(midX - wallW*0.5, gapBot, wallW, BH - gapBot);

    if (doorFrac === 0) {
      ctx.fillStyle = _c('--pink-mid');
      ctx.fillRect(midX - wallW*0.5, doorCY - 1, wallW, 2);
    }

    let meanSpd = 0;
    for (let i = 0; i < N; i++) meanSpd += Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
    meanSpd /= N;
    const spdRef = Math.max(meanSpd * 2, 0.001);
    const LX = -0.45, LY = -0.45;

    for (let i = 0; i < N; i++) {
      const spd = Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
      const hot = Math.min(spd / spdRef, 1);
      const r = Math.round(_TLR + (_PDR - _TLR) * hot);
      const g = Math.round(_TLG + (_PDG - _TLG) * hot);
      const b = Math.round(_TLB + (_PDB - _TLB) * hot);
      const hx = px[i] + LX*R*0.45, hy = py[i] + LY*R*0.45;
      const grad = ctx.createRadialGradient(hx, hy, 0, px[i], py[i], R);
      grad.addColorStop(0.0,  `rgba(255,255,255,0.85)`);
      grad.addColorStop(0.25, `rgba(${r},${g},${b},1.0)`);
      grad.addColorStop(0.75, `rgba(${Math.round(r*0.45)},${Math.round(g*0.45)},${Math.round(b*0.45)},1.0)`);
      grad.addColorStop(1.0,  `rgba(${Math.round(r*0.25)},${Math.round(g*0.25)},${Math.round(b*0.25)},1.0)`);
      ctx.beginPath();
      ctx.arc(px[i], py[i], R, 0, Math.PI*2);
      ctx.fillStyle = grad;
      ctx.fill();
    }

    /* ── Temperature readout in ctrl panel ── */
    const TL   = sideTemp(0), TR = sideTemp(1);
    const norm = Math.max(T0, 1e-12);
    const gapHalf2 = doorFrac * 3 * R;
    const gapTop2  = BH*0.5 - gapHalf2, gapBot2 = BH*0.5 + gapHalf2;
    function sideColor(s) {
      let spd = 0, count = 0;
      for (let i = 0; i < N; i++) {
        if (side[i] !== s) continue;
        if (py[i] > gapTop2 && py[i] < gapBot2) continue;
        spd += Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]); count++;
      }
      spd = count > 0 ? spd/count : 0;
      const hot = Math.min(spd / spdRef, 1);
      const r = Math.round(_TLR + (_PDR-_TLR)*hot);
      const g = Math.round(_TLG + (_PDG-_TLG)*hot);
      const b = Math.round(_TLB + (_PDB-_TLB)*hot);
      return `rgb(${r},${g},${b})`;
    }
    const elL = document.getElementById('demon-temp-left');
    const elR = document.getElementById('demon-temp-right');
    if (elL) { elL.textContent = 'T = ' + (TL/norm).toFixed(2); elL.style.color = sideColor(0); }
    if (elR) { elR.textContent = 'T = ' + (TR/norm).toFixed(2); elR.style.color = sideColor(1); }

    /* Hired stamp */
    if (Math.abs(TL - TR) / norm >= 0.5) {
      ctx.save();
      ctx.font = `bold ${Math.round(BH*0.12)}px 'EB Garamond', Georgia, serif`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.globalAlpha = 0.55;
      ctx.strokeStyle = _c('--pink-dark');
      ctx.lineWidth = 3;
      ctx.strokeText('HIRED', BW/2, BH/2);
      ctx.fillStyle = _c('--pink-dark');
      ctx.fillText('HIRED', BW/2, BH/2);
      ctx.restore();
    }
  }

  /* ── Animation loop ── */
  function loop() {
    if (running) {
      const target = doorOpen ? 1 : 0;
      if (Math.abs(doorFrac - target) > DOOR_SPEED) doorFrac += Math.sign(target - doorFrac) * DOOR_SPEED;
      else doorFrac = target;
      for (let s = 0; s < SUBSTEPS; s++) step(dt);
      render();
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:     'demon',
    title:  'Maxwell\'s Demon',
    gap:    0,
    layout: 'stacked',

    headerBtns: `<button class="applet-shell-header-btn" onclick="demonReset()">Restart</button>`,


    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn" id="demon-door-btn" onclick="demonToggleDoor()">Open Door</button>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, W, H, S }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      BW = W || S;
      BH = H || S;
      canvas.width  = BW;
      canvas.height = BH;
      midX = Math.floor(BW / 2);
      R    = Math.max(1, Math.floor(BH / 40));

      doorOpen = false; doorFrac = 0;
      const db = document.getElementById('demon-door-btn');
      if (db) { db.textContent = 'Open Door'; db.classList.remove('active'); }

      initParticles();
      running = true;
      if (!frameId) frameId = requestAnimationFrame(loop);
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    },

    onResize: function ({ W, H, S }) {
      BW = W || S; BH = H || S;
      canvas.width  = BW; canvas.height = BH;
      midX = Math.floor(BW / 2);
      R    = Math.max(1, Math.floor(BH / 40));
      initParticles();
    },
  });

  window.demonOpen  = () => shell.open();
  window.demonClose = () => shell.close();

  window.demonReset = function () {
    doorOpen = false; doorFrac = 0;
    const db = document.getElementById('demon-door-btn');
    if (db) { db.textContent = 'Open Door'; db.classList.remove('active'); }
    initParticles();
    if (!running) render();
  };

  window.demonToggleDoor = function () {
    doorOpen = !doorOpen;
    const btn = document.getElementById('demon-door-btn');
    if (btn) {
      btn.textContent = doorOpen ? 'Close Door' : 'Open Door';
      btn.classList.toggle('active', doorOpen);
    }
  };


})();

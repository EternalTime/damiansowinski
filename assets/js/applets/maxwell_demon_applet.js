(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');

  const GAP  = 0;
  const PAD  = 20;
  const PHI  = 1.6180339887;

  /* ── Simulation parameters ── */
  const N_EACH   = 50;        // particles per side
  let R          = 8;         // particle radius px — updated in layout()
  const T_INIT   = 1.0;       // initial temperature (sets speed scale)
  const DT_BASE  = 0.5;
  const SUBSTEPS = 4;
  const dt       = DT_BASE / SUBSTEPS;

  /* ── Door ── */
  // doorOpen: logical state
  // doorFrac: 0 = fully closed, 1 = fully open (animated)
  // doorGap:  pixel gap actually enforced during collision (doorFrac * doorH)
  let doorOpen   = false;
  let doorFrac   = 0;          // 0–1
  const DOOR_SPEED = 0.04;    // fraction per frame (~0.7s full travel)

  /* ── Canvas / box geometry (set in layout()) ── */
  let canvas, ctx;
  let BW, BH;    // box width (full), box height
  let midX;      // x of centre wall

  /* ── Particle arrays ── */
  let N;
  let px, py, vx, vy;
  let side;   // Int8Array: 0 = left, 1 = right — authoritative, updated only on gap crossing

  /* ── Temperature tracking ── */
  let T0 = 1.0;  // initial mean KE per particle (set after init, used for normalisation)

  /* ── Running state ── */
  let running = false;
  let frameId = null;

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
    side = new Int8Array(N);   // 0 = left, 1 = right

    // Place N_EACH on the left half, N_EACH on the right half
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

    // Zero net momentum per side, then rescale each side to exactly T_INIT
    // so both sides start at identical temperatures
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

    // T0 is now exactly T_INIT by construction
    T0 = T_INIT;
  }

  /* ── Physics step ── */
  function step(dt) {
    // Advect
    for (let i = 0; i < N; i++) {
      px[i] += vx[i] * dt;
      py[i] += vy[i] * dt;
    }

    // Wall / door collision
    const gapHalf    = doorFrac * 3 * R;
    const doorCY     = BH * 0.5;
    const gapTop     = doorCY - gapHalf;
    const gapBot     = doorCY + gapHalf;
    const wallExclude = R * 1.5;  // particle centres kept this far from midX

    for (let i = 0; i < N; i++) {
      // Outer walls — only reflect if moving into the wall
      if (px[i] - R < 0   && vx[i] < 0) { px[i] = R;      vx[i] = -vx[i]; }
      if (px[i] + R > BW  && vx[i] > 0) { px[i] = BW - R; vx[i] = -vx[i]; }
      if (py[i] - R < 0   && vy[i] < 0) { py[i] = R;      vy[i] = -vy[i]; }
      if (py[i] + R > BH  && vy[i] > 0) { py[i] = BH - R; vy[i] = -vy[i]; }

      // Centre wall: use authoritative side[] to decide reflection, never infer from px
      // Keep particle centres at least 1.5R from midX so edges can't overlap across the wall
      const inGap = py[i] > gapTop && py[i] < gapBot;
      if (!inGap) {
        if (side[i] === 0 && px[i] > midX - wallExclude) {
          px[i] = midX - wallExclude;
          if (vx[i] > 0) vx[i] = -vx[i];
        } else if (side[i] === 1 && px[i] < midX + wallExclude) {
          px[i] = midX + wallExclude;
          if (vx[i] < 0) vx[i] = -vx[i];
        }
      } else {
        // Particle is in the gap — update side if it has crossed
        if (side[i] === 0 && px[i] > midX) side[i] = 1;
        else if (side[i] === 1 && px[i] < midX) side[i] = 0;
      }
    }

    // Particle-particle collisions — purely velocity exchange, no KE injection
    const d2min = (2 * R) * (2 * R);
    for (let i = 0; i < N - 1; i++) {
      for (let j = i + 1; j < N; j++) {
        const ddx = px[j] - px[i], ddy = py[j] - py[i];
        const d2  = ddx*ddx + ddy*ddy;
        if (d2 < d2min && d2 > 1e-12) {
          const d   = Math.sqrt(d2);
          const nx  = ddx / d, ny = ddy / d;
          // Only resolve if approaching
          const dvx = vx[j] - vx[i], dvy = vy[j] - vy[i];
          const vn  = dvx * nx + dvy * ny;
          if (vn < 0) {
            // Elastic impulse (equal masses): exchange normal components
            vx[i] += vn * nx;  vy[i] += vn * ny;
            vx[j] -= vn * nx;  vy[j] -= vn * ny;
            // Correct overlap after velocity update so next step starts clean
            const ov = 2 * R - d;
            px[i] -= nx * ov * 0.5;  py[i] -= ny * ov * 0.5;
            px[j] += nx * ov * 0.5;  py[j] += ny * ov * 0.5;
            // Clamp back to authoritative side — prevents overlap nudge from
            // pushing a particle through the centre wall
            if (side[i] === 0 && px[i] > midX - wallExclude) px[i] = midX - wallExclude;
            if (side[i] === 1 && px[i] < midX + wallExclude) px[i] = midX + wallExclude;
            if (side[j] === 0 && px[j] > midX - wallExclude) px[j] = midX - wallExclude;
            if (side[j] === 1 && px[j] < midX + wallExclude) px[j] = midX + wallExclude;
          }
        }
      }
    }
  }

  /* ── Temperature calculation — excludes particles in the gap zone ── */
  function sideTemp(s) {
    const gapHalf = doorFrac * 3 * R;
    const gapTop  = BH * 0.5 - gapHalf;
    const gapBot  = BH * 0.5 + gapHalf;
    let ke = 0, count = 0;
    for (let i = 0; i < N; i++) {
      if (side[i] !== s) continue;
      if (py[i] > gapTop && py[i] < gapBot) continue;  // transiting the gap
      ke += vx[i]*vx[i] + vy[i]*vy[i];
      count++;
    }
    return count > 0 ? ke / (2 * count) : 0;
  }

  /* ── Render ── */
  function render() {
    ctx.fillStyle = _c('--bg-deep');
    ctx.fillRect(0, 0, BW, BH);

    // Centre wall
    const doorH   = BH;
    const gapHalf = doorFrac * 3 * R;   // open gap = 10% of panel height
    const doorCY  = BH * 0.5;
    const gapTop  = doorCY - gapHalf;
    const gapBot  = doorCY + gapHalf;
    const wallW   = R;

    ctx.fillStyle = _c('--teal-light');
    // Wall above gap
    if (gapTop > 0) {
      ctx.fillRect(midX - wallW * 0.5, 0, wallW, gapTop);
    }
    // Wall below gap
    if (gapBot < BH) {
      ctx.fillRect(midX - wallW * 0.5, gapBot, wallW, BH - gapBot);
    }

    // Pink seam — visible only when door is fully closed, marks the door position
    if (doorFrac === 0) {
      ctx.fillStyle = _c('--pink-mid');
      ctx.fillRect(midX - wallW * 0.5, doorCY - 1, wallW, 2);
    }

    // Particles coloured by speed relative to global mean
    let meanSpd = 0;
    for (let i = 0; i < N; i++) meanSpd += Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
    meanSpd /= N;
    const spdRef = Math.max(meanSpd * 2, 0.001);

    // Light source direction (upper-left), offset for specular highlight
    const LX = -0.45, LY = -0.45;   // normalised light direction
    for (let i = 0; i < N; i++) {
      const spd = Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
      const hot = Math.min(spd / spdRef, 1);
      const r = Math.round(_TLR + (_PDR - _TLR) * hot);
      const g = Math.round(_TLG + (_PDG - _TLG) * hot);
      const b = Math.round(_TLB + (_PDB - _TLB) * hot);

      // Radial gradient: highlight offset toward light, base colour at edge
      const hx = px[i] + LX * R * 0.45;
      const hy = py[i] + LY * R * 0.45;
      const grad = ctx.createRadialGradient(hx, hy, 0, px[i], py[i], R);
      grad.addColorStop(0.0, `rgba(255,255,255,0.85)`);         // specular peak
      grad.addColorStop(0.25, `rgba(${r},${g},${b},1.0)`);      // base colour
      grad.addColorStop(0.75, `rgba(${Math.round(r*0.45)},${Math.round(g*0.45)},${Math.round(b*0.45)},1.0)`); // shadow
      grad.addColorStop(1.0, `rgba(${Math.round(r*0.25)},${Math.round(g*0.25)},${Math.round(b*0.25)},1.0)`); // deep shadow

      ctx.beginPath();
      ctx.arc(px[i], py[i], R, 0, Math.PI * 2);
      ctx.fillStyle = grad;
      ctx.fill();
    }
  }

  /* ── Update temperature display ── */
  function updateTemps() {
    const TL = sideTemp(0);
    const TR = sideTemp(1);
    const norm = Math.max(T0, 1e-12);

    // Use same colour mapping as particles: cool rgb(37,217,200) → hot rgb(249,38,114)
    // hot = sideSpeed / (2 * globalMeanSpeed), clamped to [0,1]
    const gapHalf = doorFrac * 3 * R;
    const gapTop  = BH * 0.5 - gapHalf;
    const gapBot  = BH * 0.5 + gapHalf;
    let meanSpd = 0;
    for (let i = 0; i < N; i++) meanSpd += Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
    meanSpd /= N;
    const spdRef = Math.max(meanSpd * 2, 0.001);

    function sideColor(s) {
      let spd = 0, count = 0;
      for (let i = 0; i < N; i++) {
        if (side[i] !== s) continue;
        if (py[i] > gapTop && py[i] < gapBot) continue;
        spd += Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
        count++;
      }
      spd = count > 0 ? spd / count : 0;
      const hot = Math.min(spd / spdRef, 1);
      const r = Math.round(_TLR + (_PDR - _TLR) * hot);
      const g = Math.round(_TLG + (_PDG - _TLG) * hot);
      const b = Math.round(_TLB + (_PDB - _TLB) * hot);
      return `rgb(${r},${g},${b})`;
    }

    const colorL = sideColor(0);
    const colorR = sideColor(1);

    const labelL = document.getElementById('demon-label-left');
    const valueL = document.getElementById('demon-temp-left');
    const labelR = document.getElementById('demon-label-right');
    const valueR = document.getElementById('demon-temp-right');

    labelL.style.color = colorL;  valueL.style.color = colorL;
    labelR.style.color = colorR;  valueR.style.color = colorR;

    valueL.textContent = (TL / norm).toFixed(2);
    valueR.textContent = (TR / norm).toFixed(2);

    // Show stamp when temperature difference reaches 0.5 (normalised)
    if (Math.abs(TL - TR) / norm >= 0.5) {
      document.getElementById('demon-hired-stamp').style.display = 'block';
    }
  }

  /* ── Animation loop ── */
  function loop() {
    if (running) {
      // Animate door
      const target = doorOpen ? 1 : 0;
      if (Math.abs(doorFrac - target) > DOOR_SPEED) {
        doorFrac += Math.sign(target - doorFrac) * DOOR_SPEED;
      } else {
        doorFrac = target;
      }

      for (let s = 0; s < SUBSTEPS; s++) step(dt);
      render();
      updateTemps();
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Layout — matches gas_applet sizing ── */
  function layout() {
    const vw = window.innerWidth;
    const vh = window.innerHeight;

    // Gas applet's S is its sim-square side; its total width = S*(1+1/PHI).
    // Use that same total width as our simW, then derive simH at 2:1.
    const sFromW  = (vw - 2 * PAD - GAP) / (1 + 1 / PHI);
    const sFromH  = (vh - 2 * PAD) / 1.1;
    const S       = Math.floor(Math.min(sFromW, sFromH));
    const totalW  = Math.floor(S * (1 + 1 / PHI));  // same footprint as gas applet

    const simW   = Math.min(totalW, vw - 2 * PAD);
    const hdrH   = Math.floor(S * 0.1);
    const simH   = Math.floor(simW * 0.5);           // 2:1 wide rectangle
    const ctrlH  = hdrH;                             // same height as header strip
    const totalH = hdrH + simH + ctrlH;
    const left   = Math.floor((vw - simW) / 2);
    const top    = Math.floor((vh - totalH) / 2);

    const el = document.getElementById('demon-overlay');
    el.style.setProperty('--demon-left',     left           + 'px');
    el.style.setProperty('--demon-top-hdr',  top            + 'px');
    el.style.setProperty('--demon-top-body', (top + hdrH)   + 'px');
    el.style.setProperty('--demon-W',        simW           + 'px');
    el.style.setProperty('--demon-H-hdr',    hdrH           + 'px');
    el.style.setProperty('--demon-H-sim',    simH           + 'px');
    el.style.setProperty('--demon-H-ctrl',   ctrlH          + 'px');

    if (canvas) {
      canvas.width  = simW;
      canvas.height = simH;
    }
    BW   = simW;
    BH   = simH;
    midX = Math.floor(BW / 2);
    R    = Math.max(1, Math.floor(BH / 40));  // diameter = BH/20
  }

  /* ── Open / close ── */
  window.demonOpen = function () {
    canvas = document.getElementById('demon-canvas');
    ctx    = canvas.getContext('2d');
    layout();

    doorOpen = false;
    doorFrac = 0;
    document.getElementById('demon-door-btn').textContent = 'Open';
    document.getElementById('demon-door-btn').classList.remove('demon-active');

    initParticles();

    document.getElementById('demon-overlay').classList.add('demon-open');
    requestAnimationFrame(() => {
      document.getElementById('demon-header').classList.add('demon-open');
      document.getElementById('demon-sim-panel').classList.add('demon-open');
      document.getElementById('demon-ctrl-panel').classList.add('demon-open');
    });

    running = true;
    if (!frameId) frameId = requestAnimationFrame(loop);
  };

  window.demonClose = function () {
    running = false;
    ['demon-header', 'demon-sim-panel', 'demon-ctrl-panel'].forEach(id => {
      const el = document.getElementById(id);
      el.classList.remove('demon-open');
      el.classList.add('demon-closing');
    });
    setTimeout(() => {
      document.getElementById('demon-overlay').classList.remove('demon-open');
      ['demon-header', 'demon-sim-panel', 'demon-ctrl-panel'].forEach(id => {
        document.getElementById(id).classList.remove('demon-closing');
      });
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    }, 550);
  };

  window.demonReset = function () {
    doorOpen = false;
    doorFrac = 0;
    document.getElementById('demon-door-btn').textContent = 'Open';
    document.getElementById('demon-door-btn').classList.remove('demon-active');
    document.getElementById('demon-hired-stamp').style.display = 'none';
    initParticles();
  };

  window.demonToggleDoor = function () {
    doorOpen = !doorOpen;
    const btn = document.getElementById('demon-door-btn');
    if (doorOpen) {
      btn.textContent = 'Close';
      btn.classList.add('demon-active');
    } else {
      btn.textContent = 'Open';
      btn.classList.remove('demon-active');
    }
  };

  window.addEventListener('resize', () => {
    if (running) { layout(); initParticles(); }
  });

})();

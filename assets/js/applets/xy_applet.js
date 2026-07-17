(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  /* ── XY model parameters ── */
  let T  = 0.5;
  let L  = 48;
  let stepsPerFrame = 5;

  let spins = null;
  let cosS = null, sinS = null, wrapT = null;
  let canvasEl, ctx;
  let W = 1, H = 1;
  let running = false, frameId = null;
  let wasRunning = false;   // sim state stashed while the docs panel is open

  function idx(x, y) { return ((y + L) % L) * L + ((x + L) % L); }

  function reset() {
    spins = new Float32Array(L * L);
    for (let i = 0; i < L * L; i++) spins[i] = Math.random() * 2 * Math.PI;
    cosS  = new Float32Array(L * L);
    sinS  = new Float32Array(L * L);
    wrapT = new Int32Array(L + 2);
    for (let i = -1; i <= L; i++) wrapT[i + 1] = ((i % L) + L) % L;
  }

  function sweep() {
    const TWO_PI = 2 * Math.PI;
    const invT = 1.0 / T;
    const N = L * L;
    for (let k = 0; k < N; k++) {
      const i = (Math.random() * N) | 0;
      const x = i % L;
      const y = (i / L) | 0;
      const old_th = spins[i];
      const new_th = Math.random() * TWO_PI;
      let fx = 0, fy = 0;
      const nn = [idx(x+1,y), idx(x-1,y), idx(x,y+1), idx(x,y-1)];
      for (let n = 0; n < 4; n++) {
        fx += Math.cos(spins[nn[n]]);
        fy += Math.sin(spins[nn[n]]);
      }
      const dE = -(Math.cos(new_th) * fx + Math.sin(new_th) * fy)
                 +(Math.cos(old_th) * fx + Math.sin(old_th) * fy);
      if (dE <= 0 || Math.random() < Math.exp(-dE * invT)) spins[i] = new_th;
    }
  }

  function wrap(dth) {
    const TWO_PI = 2 * Math.PI;
    dth = dth % TWO_PI;
    if (dth >  Math.PI) dth -= TWO_PI;
    if (dth <= -Math.PI) dth += TWO_PI;
    return dth;
  }

  function vorticityAt(x, y) {
    const s00 = spins[idx(x,   y  )];
    const s10 = spins[idx(x+1, y  )];
    const s11 = spins[idx(x+1, y+1)];
    const s01 = spins[idx(x,   y+1)];
    const circulation = wrap(s10 - s00) + wrap(s11 - s10)
                      + wrap(s01 - s11) + wrap(s00 - s01);
    return Math.round(circulation / (2 * Math.PI));
  }

  const BG = _rgb('--bg-dark');
  function mute(rgb, f) {
    return [
      Math.round(rgb[0] + f * (BG[0] - rgb[0])),
      Math.round(rgb[1] + f * (BG[1] - rgb[1])),
      Math.round(rgb[2] + f * (BG[2] - rgb[2]))
    ];
  }
  const PALETTE = [
    mute(_rgb('--teal-dark'),  0.0),
    mute(_rgb('--teal-light'), 0.0),
    mute(_rgb('--cyan'),       0.0),
    mute(_rgb('--pink-dark'),  0.0),
    mute(_rgb('--pink-light'), 0.0),
    mute(_rgb('--teal-dark'),  0.0),
  ];
  const NSTOPS = PALETTE.length - 1;

  /* Angle → color LUT (replaces per-pixel angleToRGB allocation) */
  const ALUT_SIZE = 1024;               // power of two
  const alut = new Uint8Array(ALUT_SIZE * 3);
  (function buildALUT() {
    for (let i = 0; i < ALUT_SIZE; i++) {
      const seg = (i / ALUT_SIZE) * NSTOPS;
      const s = Math.floor(seg), t = seg - s;
      const c0 = PALETTE[s], c1 = PALETTE[Math.min(s + 1, NSTOPS)];
      alut[i * 3]     = Math.round(c0[0] + t * (c1[0] - c0[0]));
      alut[i * 3 + 1] = Math.round(c0[1] + t * (c1[1] - c0[1]));
      alut[i * 3 + 2] = Math.round(c0[2] + t * (c1[2] - c0[2]));
    }
  })();
  const INV_TWO_PI = 1 / (2 * Math.PI);

  let imgData = null, pix = null;

  function render() {
    const cellW = W / L, cellH = H / L;
    const arrowLen = Math.min(cellW, cellH) * 0.38;
    if (!imgData || imgData.width !== W || imgData.height !== H) {
      imgData = ctx.createImageData(W, H);
      pix = imgData.data;
      for (let i = 3; i < pix.length; i += 4) pix[i] = 255;
    }
    const NN = L * L;
    for (let i = 0; i < NN; i++) { cosS[i] = Math.cos(spins[i]); sinS[i] = Math.sin(spins[i]); }
    let off = 0;
    for (let py = 0; py < H; py++) {
      const gy = py / cellH - 0.5;
      const iy = Math.floor(gy), ty = gy - iy;
      const r0 = wrapT[iy + 1] * L, r1 = wrapT[iy + 2] * L;
      for (let px = 0; px < W; px++, off += 4) {
        const gx = px / cellW - 0.5;
        const ix = Math.floor(gx), tx = gx - ix;
        const c0 = wrapT[ix + 1], c1 = wrapT[ix + 2];
        const w00 = (1-tx)*(1-ty), w10 = tx*(1-ty), w01 = (1-tx)*ty, w11 = tx*ty;
        const bx = cosS[r0+c0]*w00 + cosS[r0+c1]*w10 + cosS[r1+c0]*w01 + cosS[r1+c1]*w11;
        const by = sinS[r0+c0]*w00 + sinS[r0+c1]*w10 + sinS[r1+c0]*w01 + sinS[r1+c1]*w11;
        let u = Math.atan2(by, bx) * INV_TWO_PI;
        if (u < 0) u += 1;
        const li = ((u * ALUT_SIZE) | 0) % ALUT_SIZE * 3;
        pix[off] = alut[li]; pix[off+1] = alut[li+1]; pix[off+2] = alut[li+2];
      }
    }
    ctx.putImageData(imgData, 0, 0);
    const vr = Math.min(cellW, cellH) * 0.42;
    const vPosRGB = mute(_rgb('--pink-dark'), 0.45), vNegRGB = mute(_rgb('--teal-light'), 0.45);
    const vPos = 'rgb(' + vPosRGB.join(',') + ')', vNeg = 'rgb(' + vNegRGB.join(',') + ')';
    ctx.strokeStyle = _c('--black');
    ctx.lineWidth = Math.max(1.5, vr * 0.25);
    for (let y = 0; y < L; y++) {
      for (let x = 0; x < L; x++) {
        const v = vorticityAt(x, y);
        if (v === 0) continue;
        ctx.beginPath();
        ctx.arc((x + 1) * cellW, (y + 1) * cellH, vr, 0, 2 * Math.PI);
        ctx.fillStyle = v > 0 ? vPos : vNeg;
        ctx.fill();
        ctx.stroke();
      }
    }
    const lw = Math.max(1.2, cellW * 0.13);
    ctx.strokeStyle = 'rgba(0,0,0,0.82)';
    ctx.lineWidth = lw;
    ctx.lineCap = 'round';
    for (let y = 0; y < L; y++) {
      for (let x = 0; x < L; x++) {
        const k = y * L + x;
        const th = spins[k];
        const cx = (x + 0.5) * cellW, cy = (y + 0.5) * cellH;
        const cosTh = cosS[k], sinTh = sinS[k];
        const ex = cx + cosTh * arrowLen, ey = cy + sinTh * arrowLen;
        const sx = cx - cosTh * arrowLen, sy = cy - sinTh * arrowLen;
        ctx.beginPath(); ctx.moveTo(sx, sy); ctx.lineTo(ex, ey); ctx.stroke();
        const headLen = arrowLen * 0.45, headAng = 0.45;
        const a1 = th + Math.PI - headAng, a2 = th + Math.PI + headAng;
        ctx.beginPath();
        ctx.moveTo(ex, ey); ctx.lineTo(ex + Math.cos(a1) * headLen, ey + Math.sin(a1) * headLen);
        ctx.moveTo(ex, ey); ctx.lineTo(ex + Math.cos(a2) * headLen, ey + Math.sin(a2) * headLen);
        ctx.stroke();
      }
    }
  }

  function loop() {
    if (running) {
      for (let s = 0; s < stepsPerFrame; s++) sweep();
      render();
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'xy',
    title: 'XY Model',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="xyReset()">Reset</button><button class="applet-shell-header-btn" id="xy-pause-btn" onclick="xyTogglePause()">Pause</button>`,


    docs: {
      whatis: `Continuous symmetries are fragile things in two dimensions. Mermin and Wagner proved as much in 1966, with a theorem carrying the air of a final word: no two-dimensional system with a continuous symmetry can order at any finite temperature [mermin1966]. The XY model is the simplest system in the theorem's crosshairs. It frees the Ising spin from its axis, so that each site of an $L \\times L$ lattice carries a unit vector free to point in any direction in the plane — an angle $\\theta_i \\in [0, 2\\pi)$ — with neighbors preferring to align:
$$E = -J \\sum_{\\langle ij \\rangle} \\cos(\\theta_i - \\theta_j).$$
Rotating every spin by the same angle costs nothing, so long-wavelength spin waves are cheap to excite, and it is exactly these that destroy spontaneous magnetization. By the logic that served so well for the Ising model, nothing should happen here at all. And yet something does, among the most subtle transitions in statistical physics.¶The resolution came in the early 1970s from Vadim Berezinskii and, independently, from Michael Kosterlitz and David Thouless [berezinskii1971, kosterlitz1973]. The key actors are not spin waves but vortices: point defects around which the spin direction winds by a full $\\pm 2\\pi$. At low temperature, vortices and antivortices bind into tight pairs whose far fields cancel, and the spins settle into quasi-long-range order: correlations decaying not exponentially but as a power law, as if the system were critical at every temperature. Above
$$k_B T_{BKT} \\approx 0.893\\,J$$
the pairs unbind [kosterlitz1974, hasenbusch2005], free vortices proliferate into a neutral plasma, and order — even the quasi kind — dissolves. This topological phase transition, invisible to any local order parameter, earned Kosterlitz and Thouless the 2016 Nobel Prize.¶The simulation evolves by single-spin updates: a random site is offered a fresh random angle, and the move is accepted with probability $\\min(1, e^{-\\Delta E/k_B T})$, driving the lattice toward the Boltzmann distribution.`,

      howto: `The canvas paints the local spin direction through the palette's color wheel (the field is smoothly interpolated between sites), with black arrows showing the individual spins. The circulation around each plaquette is monitored: pink dots mark vortices ($+2\\pi$ winding), teal dots mark antivortices ($-2\\pi$). Opposite charges attract: watch bound pairs orbit and annihilate.¶Temperature sweeps across the transition at $T_{BKT} \\approx 0.89\\,J/k_B$: run cold to see spin waves ripple over a nearly ordered lattice with only fleeting bound pairs, then heat past the transition to watch free vortices proliferate into a plasma. Quench from hot to cold and vortex–antivortex pairs will hunt each other down and annihilate as the lattice coarsens.¶Grid Size $L$ sets the lattice from $16^2$ to $96^2$ sites (changing it resets the spins), and Speed sets sweeps per frame. Reset re-randomizes the lattice; Pause freezes the dynamics for a closer look.`,

      references: ['mermin1966', 'berezinskii1971', 'kosterlitz1973', 'kosterlitz1974', 'hasenbusch2005'],
    },

    onDocsOpen:  function () { wasRunning = running; running = false; },
    onDocsClose: function () { running = wasRunning; },

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Temperature T</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Cold</span>
          <input type="range" id="xy-temp" min="0.1" max="3.0" step="0.05" value="0.5">
          <span class="applet-shell-side">Hot</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Grid Size L</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Small</span>
          <input type="range" id="xy-size" min="16" max="96" step="8" value="48">
          <span class="applet-shell-side">Large</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Speed</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Slow</span>
          <input type="range" id="xy-speed" min="1" max="20" step="1" value="5">
          <span class="applet-shell-side">Fast</span>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvasEl = c;
      ctx = canvasEl.getContext('2d');
      W = S; H = S;
      reset();
      running = true;
      const pb = document.getElementById('xy-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) loop();
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('xy-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ S }) {
      W = S; H = S;
    },
  });

  window.xyOpen  = () => shell.open();
  window.xyClose = () => shell.close();
  window.xyReset = function () { reset(); };
  window.xyTogglePause = function () {
    running = !running;
    const pb = document.getElementById('xy-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  document.getElementById('xy-temp').addEventListener('input', function () {
    T = parseFloat(this.value);
  });
  document.getElementById('xy-size').addEventListener('input', function () {
    L = parseInt(this.value); reset();
  });
  document.getElementById('xy-speed').addEventListener('input', function () {
    stepsPerFrame = parseInt(this.value);
  });

})();

(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  /* ── Inject CSS ── */
  (function () {
    if (document.getElementById('gl-styles')) return;
    const s = document.createElement('style');
    s.id = 'gl-styles';
    s.textContent = `
      #gl-ctrl-panel {
        display: flex;
        flex-direction: column;
        overflow: hidden;
      }
      .gl-plot-section {
        flex: 1;
        min-height: 0;
        display: flex;
        flex-direction: column;
        padding: 8px 12px 10px;
      }
      #gl-plot {
        flex: 1;
        min-height: 0;
        width: 100%;
        display: block;
      }
    `;
    document.head.appendChild(s);
  })();

  const PAD = 20;
  const N   = 256;
  const DT  = 0.05;

  let phi     = new Float32Array(N * N);
  let phi_old = new Float32Array(N * N);
  let phi_new = new Float32Array(N * N);

  let T = 0.025, gamma = 1.0;
  let running = false, frameId = null;
  let wasRunning = false;   // sim state stashed while the docs panel is open

  let _spare = null;
  function randn() {
    if (_spare !== null) { const v = _spare; _spare = null; return v; }
    let u, v, s;
    do { u = Math.random()*2-1; v = Math.random()*2-1; s = u*u+v*v; } while (s>=1||s===0);
    const m = Math.sqrt(-2*Math.log(s)/s);
    _spare = v * m;
    return u * m;
  }

  function idx(i, j) {
    return ((i + N) % N) * N + ((j + N) % N);
  }

  /* Precomputed wrapped neighbor offsets (avoids per-cell modulo) */
  const rowOff = new Int32Array(N), rowUp = new Int32Array(N), rowDn = new Int32Array(N);
  const colL = new Int32Array(N), colR = new Int32Array(N);
  for (let i = 0; i < N; i++) {
    rowOff[i] = i * N;
    rowUp[i]  = ((i - 1 + N) % N) * N;
    rowDn[i]  = ((i + 1) % N) * N;
    colL[i]   = (i - 1 + N) % N;
    colR[i]   = (i + 1) % N;
  }

  function initRandom() {
    for (let k = 0; k < N * N; k++) {
      phi[k]     = (Math.random() - 0.5) * 0.5;
      phi_old[k] = phi[k];
    }
  }

  function step() {
    const dt2      = DT * DT;
    const noiseMag = Math.sqrt(2 * gamma * T / dt2);
    for (let i = 0; i < N; i++) {
      const r = rowOff[i], ru = rowUp[i], rd = rowDn[i];
      for (let j = 0; j < N; j++) {
        const k = r + j;
        const p = phi[k];
        const lap = phi[rd + j] + phi[ru + j]
                  + phi[r + colR[j]] + phi[r + colL[j]]
                  - 4 * p;
        const force = lap - (-p + p*p*p);
        phi_new[k] = p + (p - phi_old[k]) * (1 - gamma * DT)
                       + dt2 * force
                       + dt2 * noiseMag * randn();
      }
    }
    const tmp = phi_old; phi_old = phi; phi = phi_new; phi_new = tmp;
  }

  /* ── Palette LUT ── */
  const LUT_SIZE = 512;
  const lut = new Uint8Array(LUT_SIZE * 3);
  (function buildLUT() {
    const stops = [
      _rgb('--pink-dark'),
      _rgb('--pink-light'),
      _rgb('--cyan'),
      _rgb('--teal-light'),
      _rgb('--teal-dark'),
    ];
    for (let i = 0; i < LUT_SIZE; i++) {
      const t  = i / (LUT_SIZE - 1);
      const ft = t * (stops.length - 1);
      const lo = Math.floor(ft);
      const hi = Math.min(lo + 1, stops.length - 1);
      const f  = ft - lo;
      lut[i*3]   = Math.round(stops[lo][0] + f * (stops[hi][0] - stops[lo][0]));
      lut[i*3+1] = Math.round(stops[lo][1] + f * (stops[hi][1] - stops[lo][1]));
      lut[i*3+2] = Math.round(stops[lo][2] + f * (stops[hi][2] - stops[lo][2]));
    }
  })();

  function phiToLUT(v) {
    const clamped = Math.max(-1.5, Math.min(1.5, v));
    const i = Math.round((clamped + 1.5) / 3.0 * (LUT_SIZE - 1));
    return i * 3;
  }

  let canvas, ctx, imgData, buf;
  let plotCanvas, plotCtx;

  const N_BINS  = 80;
  const PHI_MIN = -2, PHI_MAX = 2;
  const BIN_W   = (PHI_MAX - PHI_MIN) / N_BINS;

  const N_CURVE = 200;
  const curveX  = new Float32Array(N_CURVE);
  const curveV  = new Float32Array(N_CURVE);
  (function buildCurve() {
    let vMin = Infinity, vMax = -Infinity;
    for (let i = 0; i < N_CURVE; i++) {
      const p = PHI_MIN + (PHI_MAX - PHI_MIN) * i / (N_CURVE - 1);
      const v = -p*p/2 + p*p*p*p/4;
      curveX[i] = p; curveV[i] = v;
      if (v < vMin) vMin = v;
      if (v > vMax) vMax = v;
    }
    for (let i = 0; i < N_CURVE; i++) {
      curveV[i] = (curveV[i] - vMin) / (vMax - vMin);
    }
  })();

  function lutColor(phiVal, alpha) {
    const li = phiToLUT(phiVal);
    return `rgba(${lut[li]},${lut[li+1]},${lut[li+2]},${alpha})`;
  }

  function drawPlot() {
    if (!plotCanvas) return;
    const W = plotCanvas.width, H = plotCanvas.height;
    const PAD_L = 8, PAD_R = 8, PAD_T = 8, PAD_B = 8;
    const pw = W - PAD_L - PAD_R;
    const ph = H - PAD_T - PAD_B;
    plotCtx.clearRect(0, 0, W, H);
    const counts = new Float32Array(N_BINS);
    for (let k = 0; k < N * N; k++) {
      const bin = Math.floor((phi[k] - PHI_MIN) / BIN_W);
      if (bin >= 0 && bin < N_BINS) counts[bin]++;
    }
    const total = N * N;
    let maxP = 0;
    for (let b = 0; b < N_BINS; b++) {
      counts[b] /= (total * BIN_W);
      if (counts[b] > maxP) maxP = counts[b];
    }
    const pScale = maxP > 0 ? ph / maxP : 1;
    const barW = pw / N_BINS;
    for (let b = 0; b < N_BINS; b++) {
      const phiCenter = PHI_MIN + (b + 0.5) * BIN_W;
      const bh = counts[b] * pScale;
      plotCtx.fillStyle = lutColor(phiCenter, 0.85);
      plotCtx.fillRect(PAD_L + b * barW, PAD_T + ph - bh, Math.max(barW - 0.5, 1), bh);
    }
    plotCtx.beginPath();
    plotCtx.strokeStyle = _rgba('--text-mid', 0.85);
    plotCtx.lineWidth   = 1.5;
    for (let i = 0; i < N_CURVE; i++) {
      const x = PAD_L + (curveX[i] - PHI_MIN) / (PHI_MAX - PHI_MIN) * pw;
      const y = PAD_T + ph - curveV[i] * ph;
      if (i === 0) plotCtx.moveTo(x, y);
      else         plotCtx.lineTo(x, y);
    }
    plotCtx.stroke();
  }

  let off, offCtx;

  function render() {
    let q = 0;
    for (let k = 0; k < N * N; k++, q += 4) {
      const li = phiToLUT(phi[k]);
      buf[q] = lut[li]; buf[q+1] = lut[li+1]; buf[q+2] = lut[li+2];
    }
    offCtx.putImageData(imgData, 0, 0);
    ctx.drawImage(off, 0, 0, canvas.width, canvas.height);
  }

  let _frameCount = 0;
  function loop() {
    if (running) {
      for (let s = 0; s < 4; s++) step();
      render();
      if (++_frameCount % 2 === 0) drawPlot();
    }
    frameId = requestAnimationFrame(loop);
  }

  function resizePlotCanvas() {
    if (!plotCanvas) return;
    const section = plotCanvas.closest('.gl-plot-section');
    if (!section) return;
    plotCanvas.width  = section.clientWidth  - 24;
    plotCanvas.height = section.clientHeight - 18;
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'gl',
    title: 'Ginzburg&ndash;Landau',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="glReset()">Initialise</button><button class="applet-shell-header-btn" id="gl-pause-btn" onclick="glTogglePause()">Pause</button>`,


    docs: {
      whatis: `Every phase transition, Lev Landau argued in 1937, tells the same story: near the transition the microscopic details of a system matter far less than the symmetry it breaks, and everything worth knowing is carried by an order parameter — a coarse-grained field $\\phi(x)$ measuring the local degree of order [landau1937]. In 1950 he and Vitaly Ginzburg gave the idea its canonical form in their theory of superconductivity [ginzburg1950], and the resulting free energy has since become the hydrogen atom of field theories:
$$F[\\phi] = \\int d^2x \\, \\Big[ \\tfrac{1}{2}(\\nabla \\phi)^2 - \\tfrac{1}{2}\\phi^2 + \\tfrac{1}{4}\\phi^4 \\Big].$$
The potential is a double well: two degenerate minima at $\\phi = \\pm 1$, each a broken-symmetry phase, separated by a hump at $\\phi = 0$. Gradients are penalized, so the field wants to be uniform; the wells force it to choose a side.¶The field evolves here under Langevin dynamics — inertia, friction, and thermal noise [langevin1908]:
$$\\ddot{\\phi} = \\nabla^2 \\phi + \\phi - \\phi^3 - \\gamma \\dot{\\phi} + \\xi,$$
where the noise strength is tied to the damping through the fluctuation–dissipation relation, $\\langle \\xi(t)\\,\\xi(t') \\rangle = 2\\gamma k_B T \\, \\delta(t-t')$. Friction drains energy while the noise injects it, and the two conspire to drive the field toward the Boltzmann distribution at temperature $T$.¶Start from small random fluctuations — a quench from high temperature — and the symmetric state $\\phi = 0$ sits unstably on the hump: fluctuations grow, and the system decomposes into domains of the two phases separated by thin walls. The walls then move under their own curvature, small domains evaporating into large ones, with the characteristic size coarsening as $L(t) \\sim t^{1/2}$ [allen1979, hohenberg1977]. The same defect-forming quench, run in a cosmological setting, is the Kibble mechanism [kibble1976]: the early universe, too, had no time to choose its vacuum uniformly.`,

      howto: `The canvas shows the field on a $256 \\times 256$ periodic lattice, colored by $\\phi$: pink for one well, teal for the other, cyan along the walls where $\\phi \\approx 0$. The panel below the sliders histograms the field values, each bar colored by its $\\phi$, with the double-well potential drawn on top; watch the single peak on the hump split in two as the field falls into the wells.¶Temperature sets the thermal noise. Cold, the domains are clean and the histogram peaks are sharp; at the top of the range fluctuations rattle the domains and blur the peaks. Damping $\\gamma$ interpolates between ringing, wave-like dynamics at low $\\gamma$ — quenches launch visible ripples — and overdamped relaxation at high $\\gamma$.¶Initialise re-seeds the field with small random fluctuations, quenching the system afresh; Pause freezes the dynamics.`,

      references: ['landau1937', 'ginzburg1950', 'langevin1908', 'allen1979', 'hohenberg1977', 'kibble1976'],
    },

    onDocsOpen:  function () { wasRunning = running; running = false; },
    onDocsClose: function () { running = wasRunning; },

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Temperature</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Cold</span>
          <input type="range" id="gl-temp" min="0" max="0.1" step="0.0005" value="0.025">
          <span class="applet-shell-side">Hot</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Damping &gamma;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Low</span>
          <input type="range" id="gl-gamma" min="0.1" max="5" step="0.05" value="1.0">
          <span class="applet-shell-side">High</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section gl-plot-section">
        <canvas id="gl-plot"></canvas>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      canvas.width  = S;
      canvas.height = S;
      ctx.imageSmoothingEnabled = false;
      if (!off) {
        off        = document.createElement('canvas');
        off.width  = N;
        off.height = N;
        offCtx     = off.getContext('2d');
        imgData    = offCtx.createImageData(N, N);
        buf        = imgData.data;
        for (let i = 3; i < buf.length; i += 4) buf[i] = 255;
      }

      plotCanvas = document.getElementById('gl-plot');
      plotCtx    = plotCanvas.getContext('2d');
      setTimeout(resizePlotCanvas, 80);

      initRandom();
      running = true;
      const pb = document.getElementById('gl-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) loop();
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('gl-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ canvas: c, S }) {
      canvas = c;
      canvas.width  = S;
      canvas.height = S;
      ctx.imageSmoothingEnabled = false;
      resizePlotCanvas();
    },
  });

  window.glOpen  = () => shell.open();
  window.glClose = () => shell.close();
  window.glReset = initRandom;
  window.glTogglePause = function () {
    running = !running;
    const pb = document.getElementById('gl-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  document.getElementById('gl-temp').addEventListener('input', function () {
    T = parseFloat(this.value);
  });
  document.getElementById('gl-gamma').addEventListener('input', function () {
    gamma = parseFloat(this.value);
  });

})();

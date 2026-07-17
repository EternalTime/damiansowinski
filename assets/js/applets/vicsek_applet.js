
(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');
  const [_PLR, _PLG, _PLB] = _rgb('--pink-light');

  /* ── Vicsek parameters ── */
  const R   = 25;
  let eta   = 0.2;
  let v0    = 2.0;
  let N     = 300;
  let kNN   = 7;        // neighbours for topological mode
  let mode  = 'metric'; // 'metric' | 'topological'

  /* ── Boid state arrays ── */
  const MAX = 2500;
  const px  = new Float32Array(MAX);
  const py  = new Float32Array(MAX);
  const th  = new Float32Array(MAX);

  let canvasEl, ctx;
  let W = 1, H = 1;
  let running = false, frameId = null;
  let wasRunning = false;   // sim state stashed while the docs panel is open

  function reset() {
    const sl = document.getElementById('vicsek-count');
    if (sl) N = parseInt(sl.value);
    for (let i = 0; i < N; i++) {
      px[i] = Math.random() * W;
      py[i] = Math.random() * H;
      th[i] = (Math.random() - 0.5) * 2 * Math.PI;
    }
  }

  /* ── Shared per-step scratch (no per-step allocation) ── */
  const cosT  = new Float32Array(MAX);
  const sinT  = new Float32Array(MAX);
  const newTh = new Float32Array(MAX);

  /* ── Cell list for metric neighbor search ── */
  let cellHead = null;
  const cellNext = new Int32Array(MAX);
  let nCX = 0, nCY = 0, _cw = 1, _ch = 1;

  function buildCells() {
    nCX = Math.max(1, Math.floor(W / R));
    nCY = Math.max(1, Math.floor(H / R));
    _cw = W / nCX;
    _ch = H / nCY;
    if (!cellHead || cellHead.length < nCX * nCY) cellHead = new Int32Array(nCX * nCY);
    cellHead.fill(-1, 0, nCX * nCY);
    for (let i = 0; i < N; i++) {
      const cx = Math.min((px[i] / _cw) | 0, nCX - 1);
      const cy = Math.min((py[i] / _ch) | 0, nCY - 1);
      const c  = cy * nCX + cx;
      cellNext[i] = cellHead[c];
      cellHead[c] = i;
    }
  }

  /* ── Metric step (original Vicsek) ── */
  function stepMetric() {
    const R2 = R * R;
    const halfEtaPi = eta * Math.PI;
    for (let i = 0; i < N; i++) { cosT[i] = Math.cos(th[i]); sinT[i] = Math.sin(th[i]); }
    const useCells = Math.floor(W / R) >= 3 && Math.floor(H / R) >= 3;
    if (useCells) buildCells();
    for (let i = 0; i < N; i++) {
      let sx = 0, sy = 0;
      const xi = px[i], yi = py[i];
      if (useCells) {
        const cx = Math.min((xi / _cw) | 0, nCX - 1);
        const cy = Math.min((yi / _ch) | 0, nCY - 1);
        for (let dcy = -1; dcy <= 1; dcy++) {
          const cyw = (cy + dcy + nCY) % nCY;
          for (let dcx = -1; dcx <= 1; dcx++) {
            const cxw = (cx + dcx + nCX) % nCX;
            for (let j = cellHead[cyw * nCX + cxw]; j !== -1; j = cellNext[j]) {
              let dx = px[j] - xi, dy = py[j] - yi;
              if (dx >  W / 2) dx -= W;
              if (dx < -W / 2) dx += W;
              if (dy >  H / 2) dy -= H;
              if (dy < -H / 2) dy += H;
              if (dx * dx + dy * dy <= R2) { sx += cosT[j]; sy += sinT[j]; }
            }
          }
        }
      } else {
        for (let j = 0; j < N; j++) {
          let dx = px[j] - xi, dy = py[j] - yi;
          if (dx >  W / 2) dx -= W;
          if (dx < -W / 2) dx += W;
          if (dy >  H / 2) dy -= H;
          if (dy < -H / 2) dy += H;
          if (dx * dx + dy * dy <= R2) { sx += cosT[j]; sy += sinT[j]; }
        }
      }
      newTh[i] = Math.atan2(sy, sx) + (Math.random() - 0.5) * 2 * halfEtaPi;
    }
    for (let i = 0; i < N; i++) {
      th[i]  = newTh[i];
      px[i]  = (px[i] + v0 * Math.cos(th[i]) + W) % W;
      py[i]  = (py[i] + v0 * Math.sin(th[i]) + H) % H;
    }
  }

  /* ── Topological step (k nearest neighbours) ── */
  const _dist2 = new Float32Array(MAX);
  const _idx   = new Int32Array(MAX);
  function stepTopological() {
    const halfEtaPi = eta * Math.PI;
    for (let i = 0; i < N; i++) { cosT[i] = Math.cos(th[i]); sinT[i] = Math.sin(th[i]); }
    const k = Math.min(kNN, N - 1);
    for (let i = 0; i < N; i++) {
      const xi = px[i], yi = py[i];
      /* compute distances to all others */
      for (let j = 0; j < N; j++) {
        let dx = px[j] - xi, dy = py[j] - yi;
        if (dx >  W / 2) dx -= W;
        if (dx < -W / 2) dx += W;
        if (dy >  H / 2) dy -= H;
        if (dy < -H / 2) dy += H;
        _dist2[j] = dx * dx + dy * dy;
        _idx[j]   = j;
      }
      /* partial sort: find k nearest via selection */
      for (let m = 0; m < k; m++) {
        let minD = _dist2[m], minJ = m;
        for (let j = m + 1; j < N; j++) {
          if (_dist2[j] < minD) { minD = _dist2[j]; minJ = j; }
        }
        /* swap */
        let tmp = _dist2[m]; _dist2[m] = _dist2[minJ]; _dist2[minJ] = tmp;
        let ti  = _idx[m];   _idx[m]   = _idx[minJ];   _idx[minJ]   = ti;
      }
      /* align with k nearest (index 0 is self, skip if dist2==0) */
      let sx = 0, sy = 0;
      for (let m = 0; m < k; m++) {
        const j = _idx[m];
        sx += cosT[j]; sy += sinT[j];
      }
      newTh[i] = Math.atan2(sy, sx) + (Math.random() - 0.5) * 2 * halfEtaPi;
    }
    for (let i = 0; i < N; i++) {
      th[i] = newTh[i];
      px[i] = (px[i] + v0 * Math.cos(th[i]) + W) % W;
      py[i] = (py[i] + v0 * Math.sin(th[i]) + H) % H;
    }
  }

  function step() {
    if (mode === 'topological') stepTopological();
    else stepMetric();
  }

  const BODY = 7;
  const WING = 3.5;

  /* ── Neon look: bright white boids whose additive halo color carries the
     heading information (teal-light→pink-light by direction). ── */
  const N_HUES = 32;
  let _halos = null, _white = null;
  const _hueScratch = new Int32Array(MAX);

  function buildBoidStyles() {
    _white = _c('--white');
    _halos = new Array(N_HUES);
    const hR = BODY * 1.7;
    for (let s = 0; s < N_HUES; s++) {
      const t = s / (N_HUES - 1);
      const lr = Math.round(_TLR + t * (_PLR - _TLR));
      const lg = Math.round(_TLG + t * (_PLG - _TLG));
      const lb = Math.round(_TLB + t * (_PLB - _TLB));
      const c = document.createElement('canvas');
      c.width = c.height = Math.ceil(2 * hR) + 2;
      const hctx = c.getContext('2d');
      const cc = c.width / 2;
      const gg = hctx.createRadialGradient(cc, cc, 0, cc, cc, hR);
      gg.addColorStop(0.0, `rgba(${lr},${lg},${lb},0.6)`);
      gg.addColorStop(0.4, `rgba(${lr},${lg},${lb},0.22)`);
      gg.addColorStop(1.0, `rgba(${lr},${lg},${lb},0)`);
      hctx.beginPath(); hctx.arc(cc, cc, hR, 0, Math.PI * 2);
      hctx.fillStyle = gg; hctx.fill();
      _halos[s] = c;
    }
  }

  function render() {
    ctx.fillStyle = _c('--black');
    ctx.fillRect(0, 0, W, H);
    if (!_halos) buildBoidStyles();
    for (let i = 0; i < N; i++) {
      const t = (Math.sin(th[i]) + 1) * 0.5;
      _hueScratch[i] = (t * (N_HUES - 1) + 0.5) | 0;
    }
    /* glow pass */
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    const gHalf = _halos[0].width / 2;
    for (let i = 0; i < N; i++) {
      ctx.drawImage(_halos[_hueScratch[i]], px[i] - gHalf, py[i] - gHalf);
    }
    ctx.restore();
    /* body pass: bright white boids, one batched path */
    ctx.fillStyle = _white;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = px[i], y = py[i], a = th[i];
      const cos = Math.cos(a), sin = Math.sin(a);
      const tipX = x + cos * BODY,       tipY = y + sin * BODY;
      const b1X  = x - cos * (BODY*0.5) - sin * WING;
      const b1Y  = y - sin * (BODY*0.5) + cos * WING;
      const b2X  = x - cos * (BODY*0.5) + sin * WING;
      const b2Y  = y - sin * (BODY*0.5) - cos * WING;
      ctx.moveTo(tipX, tipY);
      ctx.lineTo(b1X, b1Y);
      ctx.lineTo(b2X, b2Y);
      ctx.closePath();
    }
    ctx.fill();
  }

  function loop() {
    if (running) { step(); render(); }
    frameId = requestAnimationFrame(loop);
  }

  function syncModeButtons() {
    document.getElementById('vicsek-btn-metric').classList.toggle('active',      mode === 'metric');
    document.getElementById('vicsek-btn-topo').classList.toggle('active',        mode === 'topological');
    const kRow = document.getElementById('vicsek-k-row');
    if (kRow) kRow.style.display = mode === 'topological' ? '' : 'none';
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'vicsek',
    title: 'Vicsek Model',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="vicsekReset()">Reset</button><button class="applet-shell-header-btn" id="vicsek-pause-btn" onclick="vicsekTogglePause()">Pause</button>`,


    docs: {
      whatis: `A flock of starlings wheeling at dusk has no leader, no plan, no bird aware of more than its handful of neighbors; and yet the whole turns as one. How does order on the scale of thousands emerge from rules on the scale of one? Tamás Vicsek and his collaborators posed the question in its barest form in 1995 [vicsek1995], stripping a bird down to a point moving at fixed speed $v_0$, carrying only a heading. At each tick every particle looks within a radius $R$, averages the directions of its neighbors, adds a dash of noise, and steps:
$$\\theta_i(t+1) = \\langle \\theta_j \\rangle_{|r_{ij}| < R} + \\eta \\, \\xi_i.$$
That is the entire model: the XY model set loose to move, alignment carried from place to place by the motion it produces.¶Turn the noise down and a disordered swarm spontaneously picks a common direction; turn it up and the consensus shatters. This is a genuine phase transition to collective motion, with the flock-averaged velocity as its order parameter, and its most remarkable feature is that it survives in two dimensions at all. The Mermin–Wagner theorem forbids exactly this kind of ordering in the equilibrium XY model, but a flock is not in equilibrium — it burns energy to move — and motion itself carries alignment information faster than fluctuations can destroy it. Whether the transition is continuous or abrupt turns out to hinge on subtle details of the noise, a question that took the better part of a decade to settle [gregoire2004].¶Real starlings, it turns out, do not use a fixed radius at all. Tracking whole flocks, Ballerini and collaborators found each bird attends to a fixed number of neighbors — roughly seven — no matter how near or far [ballerini2008]. This topological rule keeps the flock cohesive under predator attack, where a metric rule would let it tear apart. Both rules live in this applet.`,

      howto: `Each boid is a white dart trailing a halo colored by its heading; the domain is periodic, so a flock leaving one edge returns on the opposite one.¶Interaction switches the neighborhood rule. Metric averages over everyone within a fixed radius; sparse regions leave a boid nearly alone. Topological averages over the $k$ nearest neighbors regardless of distance, exposing the Neighbours slider; set it near seven for the starling value.¶Noise $\\eta$ is the control parameter: sweep it up from zero and watch an ordered flock melt into a directionless swarm at the transition. Boids sets the population, Speed sets $v_0$. Reset scatters the flock with random headings; Pause freezes it.`,

      references: ['vicsek1995', 'gregoire2004', 'ballerini2008'],
    },

    onDocsOpen:  function () { wasRunning = running; running = false; },
    onDocsClose: function () { running = wasRunning; },

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Interaction</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn active" id="vicsek-btn-metric" onclick="vicsekSetMode('metric')">Metric</button>
          <button class="applet-shell-btn"        id="vicsek-btn-topo"   onclick="vicsekSetMode('topological')">Topological</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Boids</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Few</span>
          <input type="range" id="vicsek-count" min="300" max="2500" step="50" value="300">
          <span class="applet-shell-side">Many</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Noise &eta;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Low</span>
          <input type="range" id="vicsek-noise" min="0.0" max="1.0" step="0.01" value="0.2">
          <span class="applet-shell-side">High</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section" id="vicsek-k-row" style="display:none;">
        <div class="applet-shell-ctrl-title">Neighbours</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Few</span>
          <input type="range" id="vicsek-k" min="1" max="20" step="1" value="7">
          <span class="applet-shell-side">Many</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Speed</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Slow</span>
          <input type="range" id="vicsek-speed" min="0.5" max="5.0" step="0.5" value="2.0">
          <span class="applet-shell-side">Fast</span>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvasEl = c;
      ctx = canvasEl.getContext('2d');
      W = canvasEl.clientWidth  || S;
      H = canvasEl.clientHeight || S;
      canvasEl.width  = W;
      canvasEl.height = H;
      mode = 'metric';
      syncModeButtons();
      reset();
      running = true;
      const pb = document.getElementById('vicsek-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) loop();
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('vicsek-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ canvas: c, S }) {
      W = c.clientWidth  || S;
      H = c.clientHeight || S;
      c.width  = W;
      c.height = H;
      reset();
    },
  });

  window.vicsekOpen  = () => shell.open();
  window.vicsekClose = () => shell.close();
  window.vicsekReset = function () { reset(); if (!running) render(); };
  window.vicsekTogglePause = function () {
    running = !running;
    const pb = document.getElementById('vicsek-pause-btn');
    if (pb) { pb.textContent = running ? 'Pause' : 'Resume'; pb.classList.toggle('active', !running); }
  };
  window.vicsekSetMode = function (m) {
    mode = m;
    syncModeButtons();
  };

  document.getElementById('vicsek-noise').addEventListener('input', function () {
    eta = parseFloat(this.value);
  });
  document.getElementById('vicsek-count').addEventListener('input', function () {
    reset();
  });
  document.getElementById('vicsek-k').addEventListener('input', function () {
    kNN = parseInt(this.value);
  });
  document.getElementById('vicsek-speed').addEventListener('input', function () {
    v0 = parseFloat(this.value);
  });

})();

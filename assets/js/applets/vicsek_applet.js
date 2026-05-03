
(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');

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

  function reset() {
    const sl = document.getElementById('vicsek-count');
    if (sl) N = parseInt(sl.value);
    for (let i = 0; i < N; i++) {
      px[i] = Math.random() * W;
      py[i] = Math.random() * H;
      th[i] = (Math.random() - 0.5) * 2 * Math.PI;
    }
  }

  /* ── Metric step (original Vicsek) ── */
  function stepMetric() {
    const R2 = R * R;
    const halfEtaPi = eta * Math.PI;
    const newTh = new Float32Array(N);
    for (let i = 0; i < N; i++) {
      let sx = 0, sy = 0;
      const xi = px[i], yi = py[i];
      for (let j = 0; j < N; j++) {
        let dx = px[j] - xi, dy = py[j] - yi;
        if (dx >  W / 2) dx -= W;
        if (dx < -W / 2) dx += W;
        if (dy >  H / 2) dy -= H;
        if (dy < -H / 2) dy += H;
        if (dx * dx + dy * dy <= R2) { sx += Math.cos(th[j]); sy += Math.sin(th[j]); }
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
    const newTh = new Float32Array(N);
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
        sx += Math.cos(th[j]); sy += Math.sin(th[j]);
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

  function render() {
    ctx.fillStyle = _c('--bg-dark');
    ctx.fillRect(0, 0, W, H);
    for (let i = 0; i < N; i++) {
      const x = px[i], y = py[i], a = th[i];
      const cos = Math.cos(a), sin = Math.sin(a);
      const tipX = x + cos * BODY,       tipY = y + sin * BODY;
      const b1X  = x - cos * (BODY*0.5) - sin * WING;
      const b1Y  = y - sin * (BODY*0.5) + cos * WING;
      const b2X  = x - cos * (BODY*0.5) + sin * WING;
      const b2Y  = y - sin * (BODY*0.5) - cos * WING;
      const t = (Math.sin(a) + 1) * 0.5;
      const r = Math.round(_TLR + t * (_PDR - _TLR));
      const g = Math.round(_TLG + t * (_PDG - _TLG));
      const b = Math.round(_TLB + t * (_PDB - _TLB));
      ctx.beginPath();
      ctx.moveTo(tipX, tipY);
      ctx.lineTo(b1X, b1Y);
      ctx.lineTo(b2X, b2Y);
      ctx.closePath();
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      ctx.fill();
    }
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
    title: 'Vicsek Model &mdash; Flocking Dynamics',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="vicsekReset()">Reset</button><button class="applet-shell-header-btn" id="vicsek-pause-btn" onclick="vicsekTogglePause()">Pause</button>`,


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

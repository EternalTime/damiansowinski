(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');

  /* ── Simulation constants ── */
  const N  = 400;
  const TC = 2.269;

  /* ── State ── */
  const spins = new Int8Array(N * N);
  let T              = 2.3;
  let sweepsPerFrame = 10;
  let running        = false;
  let frameId        = null;
  let lattice        = 'square';
  let J              = 1;
  let h              = 0;
  let stagger        = false;

  /* exp(-dE/T) factorized: dE = 2*(m + h*s) with m = J*sum*s ∈ [-6, 6] */
  const expJ = new Float64Array(13);
  let expHp = 1, expHm = 1;

  function buildExpTable() {
    for (let m = -6; m <= 6; m++) expJ[m + 6] = Math.exp(-2 * m / T);
    expHp = Math.exp(-2 * h / T);
    expHm = Math.exp( 2 * h / T);
  }

  function hotStart() {
    for (let i = 0; i < N * N; i++) spins[i] = Math.random() < 0.5 ? 1 : -1;
  }
  function coldStart() { spins.fill(1); }

  function sweep() {
    const tri = lattice === 'triangle';
    for (let k = 0; k < N * N; k++) {
      const i   = (Math.random() * N * N) | 0;
      const row = (i / N) | 0;
      const col = i % N;
      const s   = spins[i];
      const rU  = ((row - 1 + N) % N) * N;
      const rD  = ((row + 1)     % N) * N;
      const cL  = (col - 1 + N) % N;
      const cR  = (col + 1)     % N;
      let sum =
        spins[rU + col] + spins[rD + col] +
        spins[row * N + cL] + spins[row * N + cR];
      if (tri) {
        const cX = (row % 2 === 0) ? cR : cL;
        sum += spins[rU + cX] + spins[rD + cX];
      }
      const m  = J * sum * s;
      const dE = 2 * (m + h * s);
      if (dE <= 0 || Math.random() < expJ[m + 6] * (s === 1 ? expHp : expHm)) spins[i] = -s;
    }
  }

  /* ── Rendering ── */
  let canvas, ctx, off, offCtx, imgData, buf;

  function render() {
    let q = 0, k = 0;
    for (let row = 0; row < N; row++) {
      for (let col = 0; col < N; col++, k++, q += 4) {
        const raw  = spins[k];
        const disp = stagger ? raw * (((row + col) & 1) ? -1 : 1) : raw;
        const up   = disp === 1;
        buf[q]     = up ? _TLR : _PDR;
        buf[q + 1] = up ? _TLG : _PDG;
        buf[q + 2] = up ? _TLB : _PDB;
      }
    }
    offCtx.putImageData(imgData, 0, 0);
    ctx.drawImage(off, 0, 0, canvas.width, canvas.height);
  }

  function loop() {
    if (running) {
      for (let s = 0; s < sweepsPerFrame; s++) sweep();
      render();
    }
    frameId = requestAnimationFrame(loop);
  }

  function initCanvas(S) {
    canvas = document.getElementById('ising-canvas');
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
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'ising',
    title: 'Ising Model &mdash; Glauber Dynamics',
    gap:   0,

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Lattice</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn active" id="ising-btn-square"   onclick="isingSetLattice('square')">Square</button>
          <button class="applet-shell-btn"        id="ising-btn-triangle" onclick="isingSetLattice('triangle')">Triangular</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Coupling</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn active" id="ising-btn-ferro"     onclick="isingSetCoupling('ferro')">Ferro</button>
          <button class="applet-shell-btn"        id="ising-btn-antiferro" onclick="isingSetCoupling('antiferro')">Antiferro</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Display</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn active" id="ising-btn-display" onclick="isingToggleDisplay()">Magnetization</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Initialise</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn" onclick="isingHotStart()">Hot</button>
          <button class="applet-shell-btn" onclick="isingColdStart()">Cold</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Temperature</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Cold</span>
          <input type="range" id="ising-temp" min="0.5" max="5.0" step="0.05" value="2.5">
          <span class="applet-shell-side">Hot</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Speed</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Slow</span>
          <input type="range" id="ising-speed" min="1" max="20" step="1" value="10">
          <span class="applet-shell-side">Fast</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">External Field</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">&minus;</span>
          <div class="applet-shell-slider-wrap">
            <input type="range" id="ising-field" min="-1" max="1" step="0.02" value="0">
            <div class="applet-shell-tick" style="left:50%;"></div>
          </div>
          <span class="applet-shell-side">&plus;</span>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      initCanvas(S);
      hotStart();
      buildExpTable();
      running = true;
      if (!frameId) loop();
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    },

    onResize: function ({ S }) {
      if (!canvas) return;
      canvas.width  = S;
      canvas.height = S;
      ctx.imageSmoothingEnabled = false;
    },
  });

  window.isingOpen  = () => shell.open();
  window.isingClose = () => shell.close();

  window.isingHotStart  = hotStart;
  window.isingColdStart = coldStart;

  window.isingSetLattice = function (type) {
    lattice = type;
    buildExpTable();
    document.getElementById('ising-btn-square').classList.toggle('active',   type === 'square');
    document.getElementById('ising-btn-triangle').classList.toggle('active', type === 'triangle');
  };

  window.isingSetCoupling = function (type) {
    J = type === 'ferro' ? 1 : -1;
    document.getElementById('ising-btn-ferro').classList.toggle('active',     type === 'ferro');
    document.getElementById('ising-btn-antiferro').classList.toggle('active', type === 'antiferro');
  };

  window.isingToggleDisplay = function () {
    stagger = !stagger;
    const btn = document.getElementById('ising-btn-display');
    if (stagger) {
      btn.textContent = 'Staggered Magnetization';
      btn.classList.remove('active');
      btn.classList.add('active-pink');
    } else {
      btn.textContent = 'Magnetization';
      btn.classList.remove('active-pink');
      btn.classList.add('active');
    }
  };

  document.getElementById('ising-temp').addEventListener('input', function () {
    T = parseFloat(this.value); buildExpTable();
  });
  document.getElementById('ising-speed').addEventListener('input', function () {
    sweepsPerFrame = parseInt(this.value);
  });
  document.getElementById('ising-field').addEventListener('input', function () {
    h = parseFloat(this.value);
  });

})();

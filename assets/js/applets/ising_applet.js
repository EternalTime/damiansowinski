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

  /* ── Docs figures: square vs triangular lattice ── */
  function latticeFig(tri) {
    const a = 30, N = 5, R = 5, x0 = 15, y0 = 15;
    const dy = tri ? 26 : a;
    const H  = y0 + dy * (N - 1) + 15;
    const cols = r => (tri && (r & 1)) ? N - 1 : N;
    const pos  = (r, c) => [x0 + a * c + ((tri && (r & 1)) ? a / 2 : 0), y0 + dy * r];
    const line = (p, q) =>
      `<line x1="${p[0]}" y1="${p[1]}" x2="${q[0]}" y2="${q[1]}" style="stroke: var(--border-mid); stroke-width: 1.5;"/>`;
    const nbrs = tri
      ? [[2,1],[2,3],[1,1],[1,2],[3,1],[3,2]]
      : [[1,2],[3,2],[2,1],[2,3]];

    let bonds = '', nodes = '';
    for (let r = 0; r < N; r++) {
      for (let c = 0; c < cols(r); c++) {
        if (c < cols(r) - 1) bonds += line(pos(r, c), pos(r, c + 1));
        if (r < N - 1) {
          const targets = tri ? ((r & 1) ? [c, c + 1] : [c - 1, c]) : [c];
          targets.forEach(function (c2) {
            if (c2 >= 0 && c2 < cols(r + 1)) bonds += line(pos(r, c), pos(r + 1, c2));
          });
        }
        const isCenter = r === 2 && c === 2;
        const isNbr    = nbrs.some(n => n[0] === r && n[1] === c);
        const fill     = isCenter ? 'var(--pink-dark)' : isNbr ? 'var(--teal-light)' : 'var(--text-dim)';
        const [x, y]   = pos(r, c);
        nodes += `<circle cx="${x}" cy="${y}" r="${R}" style="fill: ${fill};"/>`;
      }
    }
    const caption = tri ? 'Triangular &mdash; 6 neighbors' : 'Square &mdash; 4 neighbors';
    return `<figure><svg width="150" height="${H}" viewBox="0 0 150 ${H}" xmlns="http://www.w3.org/2000/svg">${bonds}${nodes}</svg><figcaption>${caption}</figcaption></figure>`;
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'ising',
    title: 'Ising Model',
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

    docs: {
      whatis: `In the beginning, there was a model. Wilhelm Lenz conceived it in 1920, pondering the mystery of ferromagnetism: the simplest possible magnet, a lattice of atomic moments — spins — each pointing up ($s_i = +1$) or down ($s_i = -1$), interacting only with its nearest neighbors [lenz1920]. The entire model is contained in a single energy function,
$$E = -J\\sum_{\\langle ij \\rangle} s_i s_j - h\\sum_i s_i,$$
where the first sum runs over neighboring pairs and $h$ is an external magnetic field. For ferromagnetic coupling ($J>0$) neighbors lower the energy by aligning, while thermal fluctuations at temperature $T$ fight that order by flipping spins at random. The competition between the two is the entire story.¶Lenz handed the model to his doctoral student Ernst Ising, who solved the one-dimensional chain and found no phase transition at any finite temperature [ising1925]. Ising concluded — incorrectly — that the model had nothing to say about ferromagnetism in any dimension, and left research physics altogether. Peierls proved otherwise for two dimensions in 1936 [peierls1936]; Kramers and Wannier tickled the location of the critical point out of the model with a duality argument in 1941 [kramers1941]; and in 1944 Onsager solved the two-dimensional model exactly [onsager1944], a tour-de-force of mathematical magic confirming the transition at
$$k_B T_c = \\frac{2J}{\\ln(1+\\sqrt{2})} \\approx 2.269\\,J.$$
Below $T_c$ the lattice magnetizes; above, order dissolves. Exactly at $T_c$ fluctuations exist on all scales — patches within patches within patches — a self-similarity that is the fingerprint of a continuous phase transition.¶The simulation evolves via single-spin-flip dynamics [glauber1963]: a random spin is selected, the energy change $\\Delta E$ of flipping it is computed, and the flip is accepted with probability $\\min(1, e^{-\\Delta E/k_B T})$. One sweep is one such attempt per site. This Markov chain satisfies detailed balance, so the lattice relaxes toward the Boltzmann distribution at temperature $T$.¶Antiferromagnetic coupling ($J<0$) prefers neighbors to anti-align, forming a checkerboard on the square lattice, invisible to the ordinary magnetization but revealed by its staggered counterpart. On the triangular lattice the checkerboard is impossible: three mutually neighboring spins cannot all anti-align. This frustration destroys the ordered phase entirely, as Wannier showed in 1950 [wannier1950]: the triangular antiferromagnet remains disordered at every temperature, with a macroscopic ground-state entropy.`,

      howto: `The canvas shows a $400 \\times 400$ lattice with periodic boundaries: teal pixels are spins up, pink are down.¶Lattice selects square (4 neighbors) or triangular (6 neighbors) connectivity:¶FIG::lattices¶Coupling switches $J$ between ferromagnetic and antiferromagnetic. Display toggles between magnetization (raw spins) and staggered magnetization, which multiplies spins by a checkerboard sign so antiferromagnetic order appears as solid domains. Initialise resets the lattice: Hot scrambles every spin at random, Cold aligns them all.¶Temperature runs the lattice from deep in the ordered phase to well above $T_c \\approx 2.269\\,J/k_B$; cool slowly through the transition to watch domains coarsen, or park at $T_c$ to see fluctuations on every scale. Speed sets the number of sweeps per animation frame. External Field applies $h$, biasing spins up or down; near $T_c$ even a small field tilts the balance decisively.`,

      references: ['lenz1920', 'ising1925', 'peierls1936', 'kramers1941', 'onsager1944', 'wannier1950', 'glauber1963'],

      figures: { lattices: latticeFig(false) + latticeFig(true) },
    },

    onDocsOpen:  function () { running = false; },
    onDocsClose: function () { running = true; },

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

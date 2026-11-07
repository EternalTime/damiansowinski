(function () {
  'use strict';

  const TAU = 2 * Math.PI;
  const SOFTENING = 0.01;
  const PRESETS = {
    sync: { label: 'Static sync', J: 0.1, K: 1 },
    async: { label: 'Static async', J: 0.1, K: -1 },
    wave: { label: 'Static phase wave', J: 1, K: 0 },
    splintered: { label: 'Splintered phase wave', J: 1, K: -0.1 },
    active: { label: 'Active phase wave', J: 1, K: -0.75 },
  };

  function createModel(count, random = Math.random) {
    const n = count;
    const state = new Float64Array(3 * n);
    const first = new Float64Array(3 * n);
    const second = new Float64Array(3 * n);
    const trial = new Float64Array(3 * n);
    // Rejection sampling avoids a singular initial pair while keeping a random box.
    for (let i = 0; i < n; i++) {
      let x, y, clear;
      do {
        x = 2 * random() - 1;
        y = 2 * random() - 1;
        clear = true;
        for (let j = 0; j < i; j++) {
          if ((x - state[3 * j]) ** 2 + (y - state[3 * j + 1]) ** 2 < 0.04 ** 2) {
            clear = false;
            break;
          }
        }
      } while (!clear);
      state[3 * i] = x;
      state[3 * i + 1] = y;
      state[3 * i + 2] = TAU * random();
    }

    function rates(at, out, J, K) {
      out.fill(0);
      // Evaluate each unordered pair once, preserving reciprocal interactions.
      for (let i = 0; i < n; i++) {
        const a = 3 * i;
        for (let j = i + 1; j < n; j++) {
          const b = 3 * j;
          const dx = at[b] - at[a], dy = at[b + 1] - at[a + 1];
          const invR = 1 / Math.sqrt(dx * dx + dy * dy + SOFTENING ** 2);
          const delta = at[b + 2] - at[a + 2];
          const force = ((1 + J * Math.cos(delta)) * invR - invR * invR) / n;
          const phase = K * Math.sin(delta) * invR / n;
          out[a] += dx * force;
          out[a + 1] += dy * force;
          out[a + 2] += phase;
          out[b] -= dx * force;
          out[b + 1] -= dy * force;
          out[b + 2] -= phase;
        }
      }
      return out;
    }

    function advance(duration, J, K, maxStep = 0.01) {
      let remaining = duration;
      while (remaining > 1e-12) {
        rates(state, first, J, K);
        let fastest = 0;
        for (let i = 0; i < n; i++) {
          const a = 3 * i;
          fastest = Math.max(fastest, Math.hypot(first[a], first[a + 1]), Math.abs(first[a + 2]));
        }
        // Limit position and phase increments during close encounters.
        const dt = Math.min(remaining, maxStep, 0.01 / Math.max(1, fastest));
        for (let a = 0; a < state.length; a++) trial[a] = state[a] + dt * first[a];
        rates(trial, second, J, K);
        for (let a = 0; a < state.length; a++) state[a] += dt * (first[a] + second[a]) / 2;
        for (let i = 0; i < n; i++) state[3 * i + 2] = ((state[3 * i + 2] % TAU) + TAU) % TAU;
        remaining -= dt;
      }
    }

    function observe() {
      let cx = 0, cy = 0, c = 0, s = 0;
      for (let i = 0; i < n; i++) {
        cx += state[3 * i]; cy += state[3 * i + 1];
        c += Math.cos(state[3 * i + 2]); s += Math.sin(state[3 * i + 2]);
      }
      cx /= n; cy /= n;
      let cp = 0, sp = 0, cm = 0, sm = 0;
      for (let i = 0; i < n; i++) {
        const phi = Math.atan2(state[3 * i + 1] - cy, state[3 * i] - cx);
        const theta = state[3 * i + 2];
        cp += Math.cos(phi + theta); sp += Math.sin(phi + theta);
        cm += Math.cos(phi - theta); sm += Math.sin(phi - theta);
      }
      return { x: cx, y: cy, r: Math.hypot(c, s) / n,
        S: Math.max(Math.hypot(cp, sp), Math.hypot(cm, sm)) / n };
    }
    return { state, n, rates, advance, observe };
  }

  // The same numerical core runs in the browser and in deterministic checks.
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { createModel, PRESETS, SOFTENING };
    return;
  }

  let canvas, ctx, model, opener, overlay;
  let playing = false, docsOpen = false, frameId = null, lastTime = null;
  let J = 1, K = 0, count = 160, speed = 4, preset = 'wave', elapsed = 0;
  let palette = [], background, foreground, extent = 1.6;
  const byId = suffix => document.getElementById('swarmalator-' + suffix);
  const style = document.createElement('style');
  style.textContent = `
    #swarmalator-presets { flex-wrap: wrap; }
    #swarmalator-presets button { flex: 1 1 45%; }
    #swarmalator-readout { font-size: calc(13px * var(--shell-fs, 1)); line-height: 1.5; }
    #swarmalator-status { margin: 0; }
    #swarmalator-legend { display: block; width: 100%; height: 12px; margin-top: 8px; }
    #swarmalator-overlay :focus-visible, #swarmalator-asm-overlay :focus-visible {
      outline: 2px solid var(--teal-light); outline-offset: 2px;
    }
  `;
  document.head.appendChild(style);

  function slider(name, label, min, max, step, value, left, right) {
    return `<div class="applet-shell-ctrl-section">
      <label class="applet-shell-ctrl-title" for="swarmalator-${name}">${label}
        <output class="applet-shell-val" id="swarmalator-${name}-value" for="swarmalator-${name}">${value}</output>
      </label>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">${left}</span>
        <input type="range" id="swarmalator-${name}" min="${min}" max="${max}" step="${step}" value="${value}">
        <span class="applet-shell-side">${right}</span>
      </div></div>`;
  }

  function buildPalette() {
    const css = getComputedStyle(document.documentElement);
    const color = name => css.getPropertyValue(name).trim();
    background = color('--black'); foreground = color('--white');
    const stops = ['--teal-light', '--teal-dark', '--cyan', '--pink-dark', '--pink-light', '--teal-light']
      .map(name => {
        const hex = parseInt(color(name).replace('#', ''), 16);
        return [(hex >> 16) & 255, (hex >> 8) & 255, hex & 255];
      });
    palette = Array.from({ length: 128 }, (_, i) => {
      const u = i / 128 * (stops.length - 1), a = Math.floor(u), t = u - a;
      return 'rgb(' + stops[a].map((v, k) => Math.round(v + t * (stops[a + 1][k] - v))).join(',') + ')';
    });
    const legend = byId('legend');
    const lc = legend.getContext('2d');
    legend.width = 256; legend.height = 12;
    palette.forEach((c, i) => { lc.fillStyle = c; lc.fillRect(2 * i, 0, 2, 12); });
  }

  function syncControls() {
    for (const [name, value] of Object.entries({ J, K, count, speed })) {
      byId(name).value = value;
      byId(name + '-value').value = name === 'J' || name === 'K' ? value.toFixed(2) : String(value);
    }
    for (const key of Object.keys(PRESETS)) {
      const button = byId('preset-' + key);
      button.classList.toggle('active', key === preset);
      button.setAttribute('aria-pressed', String(key === preset));
    }
    byId('regime').textContent = preset ? 'Preset: ' + PRESETS[preset].label : 'Custom parameters';
  }

  function render() {
    if (!model || !ctx) return;
    const metrics = model.observe();
    let needed = 1.6;
    for (let i = 0; i < model.n; i++) {
      needed = Math.max(needed, 1.15 * Math.hypot(model.state[3 * i] - metrics.x, model.state[3 * i + 1] - metrics.y));
    }
    extent = needed;
    const W = canvas.width, H = canvas.height;
    const scale = Math.min(W, H) / (2 * extent);
    ctx.fillStyle = background;
    ctx.fillRect(0, 0, W, H);
    const radius = Math.max(2.5, Math.min(4.5, Math.min(W, H) / 100));
    for (let i = 0; i < model.n; i++) {
      const a = 3 * i, theta = model.state[a + 2];
      const x = W / 2 + (model.state[a] - metrics.x) * scale;
      const y = H / 2 - (model.state[a + 1] - metrics.y) * scale;
      ctx.fillStyle = palette[Math.floor(theta / TAU * palette.length) % palette.length];
      ctx.shadowColor = ctx.fillStyle;
      ctx.shadowBlur = radius * 2;
      ctx.beginPath(); ctx.arc(x, y, radius, 0, TAU); ctx.fill();
      ctx.shadowBlur = 0;
      // The white phase hand also carries phase information without color.
      ctx.strokeStyle = foreground; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x, y);
      ctx.lineTo(x + radius * Math.cos(theta), y - radius * Math.sin(theta)); ctx.stroke();
    }
    byId('metrics').textContent = `t = ${elapsed.toFixed(1)}   r = ${metrics.r.toFixed(2)}   S = ${metrics.S.toFixed(2)}`;
    byId('view').textContent = `View width: ${(2 * extent).toFixed(1)} units`;
  }

  function schedule() {
    if (playing && !docsOpen && frameId === null) frameId = requestAnimationFrame(loop);
  }
  function haltFrame() {
    if (frameId !== null) cancelAnimationFrame(frameId);
    frameId = null; lastTime = null;
  }
  function setPlaying(value) {
    playing = value; haltFrame();
    const button = byId('pause-btn');
    button.textContent = playing ? 'Pause' : 'Play';
    button.classList.toggle('active', !playing);
    byId('status').textContent = playing ? 'Running' : 'Paused';
    schedule();
  }
  function reset() {
    model = createModel(count); elapsed = 0; lastTime = null;
    render();
    byId('status').textContent = playing ? 'Reset. Running' : 'Reset. Paused';
  }
  function loop(now) {
    frameId = null;
    if (!playing || docsOpen) return;
    if (lastTime !== null) {
      // Drop excess wall time after a hidden tab or a slow frame, avoiding catch-up bursts.
      const dt = Math.min((now - lastTime) / 1000, 0.05) * speed;
      model.advance(dt, J, K); elapsed += dt;
    }
    lastTime = now;
    render(); schedule();
  }

  function onKeyDown(event) {
    if (event.key === 'Escape') { event.preventDefault(); shell.close(); return; }
    if (event.key !== 'Tab') return;
    const focusable = Array.from(overlay.querySelectorAll('button, input, a[href]'))
      .filter(el => !el.disabled && el.getClientRects().length && getComputedStyle(el).visibility !== 'hidden');
    const first = focusable[0], last = focusable[focusable.length - 1];
    if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last.focus(); }
    else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first.focus(); }
  }

  const shell = new AppletShell({
    id: 'swarmalator', title: 'Swarmalator', gap: 0,
    headerBtns: `<button class="applet-shell-header-btn" id="swarmalator-reset-btn">Reset</button><button class="applet-shell-header-btn" id="swarmalator-pause-btn">Pause</button>`,
    docs: {
      whatis: `A swarmalator moves through the plane while carrying an internal phase.
Its phase affects whom it approaches, and its neighbors' distances affect how its phase changes.
Kevin P. O'Keeffe, Hyunsuk Hong, and Steven H. Strogatz studied this coupling in 2017 [okeeffe2017swarm].¶We use their identical planar model, with unit attraction and repulsion, in a frame where the common velocity and natural frequency vanish.
For particle $i$, position $\\mathbf{x}_i$, phase $\\theta_i$, population $N$, and separation $r_{ij}=|\\mathbf{x}_j-\\mathbf{x}_i|$, the equations are [okeeffe2017swarm]
$$\\dot{\\mathbf{x}}_i=\\frac{1}{N}\\sum_{j\\ne i}\\left[\\frac{\\mathbf{x}_j-\\mathbf{x}_i}{r_{ij}}(1+J\\cos(\\theta_j-\\theta_i))-\\frac{\\mathbf{x}_j-\\mathbf{x}_i}{r_{ij}^2}\\right],$$
$$\\dot{\\theta}_i=\\frac{K}{N}\\sum_{j\\ne i}\\frac{\\sin(\\theta_j-\\theta_i)}{r_{ij}}.$$
Positive $J$ strengthens attraction between similar phases; negative $J$ favors opposite phases.
Positive $K$ promotes synchrony; negative $K$ separates phases [okeeffe2017swarm].¶The presets use the published $(J,K)$ pairs: static sync $(0.1,1)$, static async $(0.1,-1)$, static phase wave $(1,0)$, splintered phase wave $(1,-0.1)$, and active phase wave $(1,-0.75)$.
In a static phase wave, phases sort around an annulus; splintered waves form oscillating clusters, and active waves have circulating particles [okeeffe2017swarm].¶The readouts measure synchrony $r$ and phase-position correlation $S$:
$$r=\\left|\\frac{1}{N}\\sum_j e^{i\\theta_j}\\right|,\\qquad S=\\max_{\\pm}\\left|\\frac{1}{N}\\sum_j e^{i(\\phi_j\\pm\\theta_j)}\\right|.$$
Here $\\phi_j$ is the spatial angle about the swarm's center; $r=1$ means equal phases, while large $S$ indicates angular sorting [okeeffe2017swarm].`,
      howto: `Color and the white hand on each particle show its phase, from $0$ to $2\\pi$.
The view follows the swarm's center and scales to contain it; View width is in model units.¶A preset sets $J$ and $K$ and scatters the particles anew.
The labels name the published regimes; the pattern develops through a transient that depends on the random start and population, especially for splintered waves [okeeffe2017swarm].
Changing either coupling keeps the current swarm and selects Custom parameters.¶Particles changes the population and resets it.
Simulation speed sets model time per second; $t$ records elapsed model time.
Reset scatters positions and phases while keeping the controls and play state.
Close and reopen starts a fresh swarm with the selected controls; Pause freezes it and Play resumes it.
Opening these notes pauses evolution, and returning restores the play state.
Reduced motion starts the applet paused.¶For numerical evolution we replace each $r_{ij}$ by $\\sqrt{r_{ij}^2+0.01^2}$ and use Heun steps of at most $0.01$ model time, shortened during fast motion.
Initial positions lie in a square of side $2$, with a minimum separation of $0.04$; initial phases are random.
These are the applet's numerical choices.`,
      references: ['okeeffe2017swarm'],
    },
    ctrlHTML: `<div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title" id="swarmalator-presets-label">Preset regimes</div>
      <div class="applet-shell-btn-row" id="swarmalator-presets" role="group" aria-labelledby="swarmalator-presets-label">
      ${Object.entries(PRESETS).map(([key, p]) => `<button class="applet-shell-btn" id="swarmalator-preset-${key}" aria-pressed="false">${p.label}</button>`).join('')}
      </div></div>
      ${slider('J', 'Spatial coupling J', -0.9, 1, 0.05, 1, '-0.9', '1')}
      ${slider('K', 'Phase coupling K', -1, 1, 0.05, 0, '-1', '1')}
      ${slider('count', 'Particles', 80, 240, 20, 160, '80', '240')}
      ${slider('speed', 'Simulation speed', 0.5, 10, 0.5, 4, '0.5', '10')}
      <div class="applet-shell-ctrl-section" id="swarmalator-readout">
        <div id="swarmalator-regime"></div>
        <div id="swarmalator-metrics" aria-live="off"></div>
        <div id="swarmalator-view"></div>
        <p id="swarmalator-status" role="status"></p>
        <canvas id="swarmalator-legend" role="img" aria-label="Phase palette, 0 to 2 pi; white hands show the same phase angle"></canvas>
        <div>Phase: 0 to 2&pi;</div>
      </div>`,
    onOpen: function ({ canvas: c }) {
      opener = document.activeElement;
      canvas = c; ctx = c.getContext('2d');
      canvas.setAttribute('role', 'img');
      canvas.setAttribute('aria-label', 'Planar swarmalators. Color and white hand show phase. Readouts give synchrony r and phase-position correlation S.');
      overlay = byId('overlay') || byId('asm-overlay');
      overlay.setAttribute('role', 'dialog'); overlay.setAttribute('aria-modal', 'true');
      overlay.setAttribute('aria-label', 'Swarmalator');
      docsOpen = false; buildPalette(); syncControls(); reset();
      setPlaying(!window.matchMedia('(prefers-reduced-motion: reduce)').matches);
      byId('pause-btn').focus();
      document.addEventListener('keydown', onKeyDown);
    },
    onClose: function () {
      setPlaying(false); docsOpen = false;
      overlay.removeAttribute('aria-modal');
      document.removeEventListener('keydown', onKeyDown);
      if (opener && opener.isConnected) opener.focus();
    },
    onResize: function () { render(); },
    onDocsOpen: function () { docsOpen = true; haltFrame(); },
    onDocsClose: function () { docsOpen = false; schedule(); },
  });

  byId('reset-btn').addEventListener('click', reset);
  byId('pause-btn').addEventListener('click', () => setPlaying(!playing));
  for (const [key, p] of Object.entries(PRESETS)) {
    byId('preset-' + key).addEventListener('click', function () {
      J = p.J; K = p.K; preset = key; syncControls(); reset();
    });
  }
  for (const name of ['J', 'K', 'count', 'speed']) {
    byId(name).addEventListener('input', function () {
      const value = Number(this.value);
      if (name === 'J') { J = value; preset = null; }
      if (name === 'K') { K = value; preset = null; }
      if (name === 'speed') speed = value;
      if (name === 'count') count = value;
      syncControls();
      if (name === 'count') reset();
    });
  }
  window.swarmalatorOpen = () => shell.open();
  window.swarmalatorClose = () => shell.close();
})();

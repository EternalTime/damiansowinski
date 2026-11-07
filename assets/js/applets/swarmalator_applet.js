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
      return { x: cx / n, y: cy / n, r: Math.hypot(c, s) / n };
    }
    return { state, n, rates, advance, observe };
  }

  // The same numerical core runs in the browser and in deterministic checks.
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { createModel, PRESETS, SOFTENING };
    return;
  }

  let canvas, ctx, model;
  let running = false, frameId = null, lastTime = null;
  const shell = new AppletShell({
    id: 'swarmalator', title: 'Swarmalator', gap: 0,
    onOpen: function ({ canvas: c }) {
      canvas = c; ctx = c.getContext('2d'); model = createModel(160);
      running = true; lastTime = null;
      if (frameId === null) frameId = requestAnimationFrame(loop);
    },
    onClose: function () {
      running = false;
      cancelAnimationFrame(frameId); frameId = null;
    },
  });
  function loop(now) {
    if (running) {
      if (lastTime !== null) model.advance(Math.min((now - lastTime) / 1000, 0.05), 1, 0);
      lastTime = now;
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      for (let i = 0; i < model.n; i++) {
        const a = 3 * i;
        ctx.fillStyle = 'white';
        ctx.beginPath();
        ctx.arc(canvas.width / 2 + model.state[a] * canvas.width / 4,
          canvas.height / 2 - model.state[a + 1] * canvas.height / 4, 3, 0, TAU);
        ctx.fill();
      }
    }
    frameId = requestAnimationFrame(loop);
  }
  window.swarmalatorOpen = () => shell.open();
  window.swarmalatorClose = () => shell.close();
})();

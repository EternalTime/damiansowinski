/* ─────────────────────────────────────────────────────────────────────────────
   AppletShell  —  shared scaffold for all standard-layout physics applets.

   Usage:
     const shell = new AppletShell({
       id:       'ising',                          // prefix used for all IDs
       title:    'Ising Model \u2014 Glauber Dynamics',  // header text (HTML allowed)
       gap:      0,                                // px gap between sim and ctrl panels
       onOpen:   ({ canvas, S }) => { ... },       // called after panels animate in
       onClose:  ()              => { ... },       // called before panels animate out
       onResize: ({ canvas, S }) => { ... },       // called on window resize (optional)
     });
     window.isingOpen  = () => shell.open();
     window.isingClose = () => shell.close();

   The shell injects:
     - A single shared <style> block (first call only) with all visual constants.
     - Per-applet <style> scoped to the applet's prefix.
     - The overlay + three-panel HTML scaffold.

   Visual constants (change here → all applets update):
     Colours, border radii, transition curve, glow values,
     button styles, slider thumb, section dividers.
───────────────────────────────────────────────────────────────────────────── */

(function () {
  'use strict';

  /* ── Layout constants ──────────────────────────────────────────────────── */
  const PHI        = 1.6180339887;
  const PAD        = 20;    // viewport margin (px)
  const FS_REF     = 300;   // ctrl panel width at which font scale = 1.0
  const FS_MIN     = 0.65;  // never shrink below 65% of reference sizes
  const FS_MAX     = 1.15;  // never grow above 115%

  /* ── Shared CSS — injected once ────────────────────────────────────────── */
  const SHARED_STYLE_ID = 'applet-shell-shared-styles';

  function injectSharedStyles() {
    if (document.getElementById(SHARED_STYLE_ID)) return;
    const s = document.createElement('style');
    s.id = SHARED_STYLE_ID;
    s.textContent = `
/* ── Per-panel glow colours (fixed for all applets) ── */
.applet-shell-header  {
  border-color: var(--pink-dark);
  box-shadow: 0 -8px 30px rgba(var(--pink-dark-rgb), 0.45),
              0  8px 20px rgba(var(--pink-dark-rgb), 0.25);
}
.applet-shell-sim {
  border-color: var(--teal-dark);
  box-shadow: -12px 0 35px rgba(var(--teal-dark-rgb), 0.5),
                0 12px 30px rgba(var(--teal-dark-rgb), 0.3);
}
.applet-shell-ctrl {
  border-color: var(--cyan);
  box-shadow: 12px 0 35px rgba(var(--cyan-rgb), 0.5),
               0 12px 30px rgba(var(--cyan-rgb), 0.3);
  clip-path: inset(0px -60px -60px -60px);
}

/* ── Shared panel base ── */
.applet-shell-panel {
  position: fixed;
  z-index: 910;
  background: var(--bg-dark);
  color: var(--text-bright);
  font-family: 'EB Garamond', Georgia, serif;
  border-width: 2px;
  border-style: solid;
  border-radius: 14px;
  overflow: hidden;
  transition: transform 0.55s cubic-bezier(0.22, 1, 0.36, 1),
              opacity  0.0s ease;
  opacity: 0;
  pointer-events: none;
  box-sizing: border-box;
}
.applet-shell-panel.applet-shell-open {
  opacity: 1;
  pointer-events: auto;
}
.applet-shell-panel.applet-shell-closing {
  opacity: 1;
  pointer-events: none;
}

/* ── Header ── */
.applet-shell-header-inner {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0 18px;
  height: 100%;
}
.applet-shell-header-inner h2 {
  font-size: calc(15px * var(--shell-fs, 1));
  letter-spacing: 2px;
  color: var(--text-bright);
  font-weight: normal;
  margin: 0;
  text-transform: uppercase;
}
.applet-shell-close-btn {
  background: none;
  border: 1px solid var(--border-mid);
  color: var(--text-bright);
  font-family: 'EB Garamond', Georgia, serif;
  font-size: calc(13px * var(--shell-fs, 1));
  letter-spacing: 1px;
  padding: calc(3px * var(--shell-fs, 1)) calc(12px * var(--shell-fs, 1));
  border-radius: 4px;
  cursor: pointer;
}
.applet-shell-close-btn:hover { background: var(--bg-control); }

/* ── Canvas fills sim panel ── */
.applet-shell-canvas {
  flex: 1;
  display: block;
  width: 100%;
  min-height: 0;
  image-rendering: pixelated;
}

/* ── Control panel sections ── */
.applet-shell-ctrl-section {
  padding: calc(14px * var(--shell-fs, 1)) calc(16px * var(--shell-fs, 1)) calc(10px * var(--shell-fs, 1));
  border-bottom: 1px solid var(--border-dark);
  display: flex;
  flex-direction: column;
  gap: calc(10px * var(--shell-fs, 1));
}
.applet-shell-ctrl-section:last-child { border-bottom: none; }

.applet-shell-ctrl-title {
  font-size: calc(14px * var(--shell-fs, 1));
  text-transform: uppercase;
  letter-spacing: 1.5px;
  color: var(--text-dim);
  margin: 0 0 2px;
}

/* ── Slider rows ── */
.applet-shell-slider-row {
  display: flex;
  align-items: center;
  gap: calc(8px * var(--shell-fs, 1));
  width: 100%;
  box-sizing: border-box;
}
.applet-shell-slider-row .applet-shell-side {
  font-size: calc(14px * var(--shell-fs, 1));
  letter-spacing: 0.5px;
  color: var(--text-dim);
  white-space: nowrap;
  flex-shrink: 0;
}
.applet-shell-slider-row input[type=range] {
  flex: 1;
  min-width: 0;
  -webkit-appearance: none;
  appearance: none;
  height: 3px;
  border-radius: 2px;
  background: rgba(var(--slider-track-rgb), 0.35);
  outline: none;
  cursor: pointer;
}
.applet-shell-slider-row input[type=range]::-webkit-slider-thumb {
  -webkit-appearance: none;
  width: calc(12px * var(--shell-fs, 1));
  height: calc(12px * var(--shell-fs, 1));
  border-radius: 50%;
  background: var(--teal-light);
  cursor: pointer;
  margin-top: calc(-4.5px * var(--shell-fs, 1));
}
.applet-shell-slider-row input[type=range]::-moz-range-thumb {
  width: calc(12px * var(--shell-fs, 1));
  height: calc(12px * var(--shell-fs, 1));
  border-radius: 50%;
  background: var(--teal-light);
  border: none;
  cursor: pointer;
}
.applet-shell-val {
  font-size: calc(16px * var(--shell-fs, 1));
  color: var(--teal-light);
  min-width: calc(34px * var(--shell-fs, 1));
  text-align: right;
  letter-spacing: 0.5px;
  flex-shrink: 0;
}

/* ── Buttons ── */
.applet-shell-btn-row {
  display: flex;
  gap: calc(8px * var(--shell-fs, 1));
  flex-wrap: wrap;
}
.applet-shell-btn {
  background: var(--bg-control);
  color: var(--text-bright);
  border: 1px solid var(--border-mid);
  padding: calc(5px * var(--shell-fs, 1)) calc(14px * var(--shell-fs, 1));
  cursor: pointer;
  font-family: 'EB Garamond', Georgia, serif;
  font-size: calc(16px * var(--shell-fs, 1));
  letter-spacing: 1px;
  border-radius: 4px;
  transition: background 0.15s;
}
.applet-shell-btn:hover { background: var(--bg-hover); }
.applet-shell-btn.active {
  background: var(--bg-active-teal);
  border-color: var(--teal-light);
  color: var(--teal-light);
  box-shadow: 0 0 10px rgba(var(--teal-light-rgb), 0.6),
              0 0 20px rgba(var(--teal-light-rgb), 0.3);
}
.applet-shell-btn.active-pink {
  background: var(--bg-active-pink2);
  border-color: var(--pink-dark);
  color: var(--pink-dark);
  box-shadow: 0 0 10px rgba(var(--pink-dark-rgb), 0.6),
              0 0 20px rgba(var(--pink-dark-rgb), 0.3);
}
    `;
    document.head.appendChild(s);
  }

  /* ── Per-applet CSS — positioning via CSS custom properties ─────────────── */
  function injectAppletStyles(id) {
    const styleId = 'applet-shell-styles-' + id;
    if (document.getElementById(styleId)) return;
    const s = document.createElement('style');
    s.id = styleId;
    s.textContent = `
#${id}-overlay {
  display: none;
  position: fixed;
  inset: 0;
  z-index: 900;
  pointer-events: none;
}
#${id}-overlay.${id}-open {
  display: block;
  pointer-events: auto;
}

/* Header — slides from top */
#${id}-header {
  left:   var(--${id}-left);
  top:    var(--${id}-top-hdr);
  width:  var(--${id}-W);
  height: var(--${id}-H-hdr);
  transform: translateY(-120px);
  border-bottom-left-radius: 0;
  border-bottom-right-radius: 0;
}
#${id}-header.applet-shell-open { transform: translateY(0); }

/* Sim panel — slides from left */
#${id}-sim-panel {
  left:   var(--${id}-left);
  top:    var(--${id}-top-body);
  width:  var(--${id}-W-sim);
  height: var(--${id}-H-body);
  display: flex;
  flex-direction: column;
  transform: translateX(-110vw);
  border-top-left-radius: 0;
  border-top-right-radius: 0;
  border-bottom-right-radius: 0;
}
#${id}-sim-panel.applet-shell-open { transform: translateX(0); }

/* Ctrl panel — slides from right */
#${id}-ctrl-panel {
  top:    var(--${id}-top-body);
  left:   calc(var(--${id}-left) + var(--${id}-W-sim) + var(--${id}-gap, 0px));
  width:  var(--${id}-W-ctrl);
  height: var(--${id}-H-body);
  display: flex;
  flex-direction: column;
  transform: translateX(110vw);
  border-top-left-radius: 0;
  border-top-right-radius: 0;
  border-bottom-left-radius: 0;
  overflow-y: auto;
}
#${id}-ctrl-panel.applet-shell-open { transform: translateX(0); }
    `;
    document.head.appendChild(s);
  }

  /* ── HTML scaffold ──────────────────────────────────────────────────────── */
  function buildScaffold(id, title, ctrlHTML) {
    const div = document.createElement('div');
    div.innerHTML = `
<div id="${id}-overlay">

  <div id="${id}-header" class="applet-shell-panel applet-shell-header">
    <div class="applet-shell-header-inner">
      <h2>${title}</h2>
      <button class="applet-shell-close-btn" data-shell-close="${id}">Close</button>
    </div>
  </div>

  <div id="${id}-sim-panel" class="applet-shell-panel applet-shell-sim">
    <canvas id="${id}-canvas" class="applet-shell-canvas"></canvas>
  </div>

  <div id="${id}-ctrl-panel" class="applet-shell-panel applet-shell-ctrl">
    ${ctrlHTML}
  </div>

</div>
    `.trim();
    return div.firstChild;
  }

  /* ── Layout computation ─────────────────────────────────────────────────── */
  function computeLayout(gap) {
    const vw = window.innerWidth;
    const vh = window.innerHeight;
    const sFromW = (vw - 2 * PAD - gap) / (1 + 1 / PHI);
    const sFromH = (vh - 2 * PAD) / 1.1;
    const S      = Math.floor(Math.min(sFromW, sFromH));
    const ctrlW  = Math.floor(S / PHI);
    const hdrH   = Math.floor(S * 0.1);
    const totalW = S + gap + ctrlW;
    const left   = Math.floor((vw - totalW) / 2);
    const top    = Math.floor((vh - (hdrH + S)) / 2);
    return { S, ctrlW, hdrH, totalW, left, top };
  }

  function applyLayout(id, gap) {
    const { S, ctrlW, hdrH, totalW, left, top } = computeLayout(gap);
    const fs = Math.min(FS_MAX, Math.max(FS_MIN, ctrlW / FS_REF));
    const el = document.getElementById(id + '-overlay');
    el.style.setProperty('--' + id + '-left',     left           + 'px');
    el.style.setProperty('--' + id + '-top-hdr',  top            + 'px');
    el.style.setProperty('--' + id + '-top-body', (top + hdrH)   + 'px');
    el.style.setProperty('--' + id + '-W',        totalW         + 'px');
    el.style.setProperty('--' + id + '-W-sim',    S              + 'px');
    el.style.setProperty('--' + id + '-W-ctrl',   ctrlW          + 'px');
    el.style.setProperty('--' + id + '-H-hdr',    hdrH           + 'px');
    el.style.setProperty('--' + id + '-H-body',   S              + 'px');
    el.style.setProperty('--' + id + '-gap',      gap            + 'px');
    el.style.setProperty('--shell-fs',             fs.toFixed(4));
    return S;
  }

  /* ── AppletShell constructor ─────────────────────────────────────────────── */
  function AppletShell(cfg) {
    const id       = cfg.id;
    const title    = cfg.title;
    const gap      = cfg.gap || 0;
    const onOpen   = cfg.onOpen   || function () {};
    const onClose  = cfg.onClose  || function () {};
    const onResize = cfg.onResize || null;
    const ctrlHTML = cfg.ctrlHTML || '';

    // Inject styles
    injectSharedStyles();
    injectAppletStyles(id);

    // Build and insert HTML scaffold
    const scaffold = buildScaffold(id, title, ctrlHTML);
    document.body.appendChild(scaffold);

    // Wire up close button
    scaffold.querySelector('[data-shell-close]').addEventListener('click', function () {
      self.close();
    });

    const panelIds = [id + '-header', id + '-sim-panel', id + '-ctrl-panel'];

    function layout() {
      const S      = applyLayout(id, gap);
      const canvas = document.getElementById(id + '-canvas');
      if (canvas) {
        canvas.width  = S;
        canvas.height = S;
      }
      return S;
    }

    const self = {
      open: function () {
        const S      = layout();
        const canvas = document.getElementById(id + '-canvas');

        document.getElementById(id + '-overlay').classList.add(id + '-open');
        requestAnimationFrame(function () {
          panelIds.forEach(function (pid) {
            document.getElementById(pid).classList.add('applet-shell-open');
          });
        });

        onOpen({ canvas: canvas, S: S });
      },

      close: function () {
        onClose();

        panelIds.forEach(function (pid) {
          const el = document.getElementById(pid);
          el.classList.remove('applet-shell-open');
          el.classList.add('applet-shell-closing');
        });

        setTimeout(function () {
          document.getElementById(id + '-overlay').classList.remove(id + '-open');
          panelIds.forEach(function (pid) {
            document.getElementById(pid).classList.remove('applet-shell-closing');
          });
        }, 550);
      },
    };

    // Resize handler
    window.addEventListener('resize', function () {
      const overlay = document.getElementById(id + '-overlay');
      if (!overlay.classList.contains(id + '-open')) return;
      const S      = layout();
      const canvas = document.getElementById(id + '-canvas');
      if (onResize) onResize({ canvas: canvas, S: S });
    });

    return self;
  }

  window.AppletShellDesktop = AppletShell;
  window.AppletShell = AppletShell;

})();

/* ─────────────────────────────────────────────────────────────────────────────
   AppletShellMobile  —  full-screen mobile scaffold for physics applets.

   Layout (portrait and landscape):
     ┌─────────────────────┐
     │       HEADER        │  48px, slides from top
     ├─────────────────────┤
     │                     │
     │    SIM (square)     │  min(vw, vh*0.55), slides from left
     │                     │
     ├─────────────────────┤
     │   CONTROLS (scroll) │  remaining height, slides from bottom
     └─────────────────────┘

   Interface: identical to AppletShell — id, title, ctrlHTML,
   onOpen, onClose, onResize. Individual applet files unchanged.
───────────────────────────────────────────────────────────────────────────── */

(function () {
  'use strict';

  const HDR_H    = 48;    // px — header height
  const SIM_FRAC = 0.55;  // fraction of vh used to cap sim size
  const FS_BASE  = 375;   // reference vw for font scaling
  const FS_MIN   = 0.85;  // minimum font scale
  const FS_MAX   = 1.30;  // maximum font scale

  /* ── Shared CSS — injected once ────────────────────────────────────────── */
  const SHARED_STYLE_ID = 'applet-shell-mobile-shared-styles';

  function injectSharedStyles() {
    if (document.getElementById(SHARED_STYLE_ID)) return;
    const s = document.createElement('style');
    s.id = SHARED_STYLE_ID;
    s.textContent = `
/* ── Mobile panel glows ── */
.asm-header {
  border-color: var(--pink-dark);
  box-shadow: 0 4px 20px rgba(var(--pink-dark-rgb), 0.45);
}
.asm-sim {
  border-color: var(--teal-dark);
  box-shadow: 0 4px 20px rgba(var(--teal-dark-rgb), 0.4);
}
.asm-ctrl {
  border-color: var(--cyan);
  box-shadow: 0 -4px 20px rgba(var(--cyan-rgb), 0.4);
}

/* ── Base panel ── */
.asm-panel {
  position: fixed;
  z-index: 910;
  background: var(--bg-dark);
  color: var(--text-bright);
  font-family: 'EB Garamond', Georgia, serif;
  border-width: 2px;
  border-style: solid;
  border-radius: 0;
  box-sizing: border-box;
  opacity: 0;
  pointer-events: none;
  transition: transform 0.45s cubic-bezier(0.22, 1, 0.36, 1);
}
.asm-panel.asm-open {
  opacity: 1;
  pointer-events: auto;
}
.asm-panel.asm-closing {
  opacity: 1;
  pointer-events: none;
}

/* ── Header inner ── */
.asm-header-inner {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0 16px;
  height: 100%;
}
.asm-header-inner h2 {
  font-size: calc(13px * var(--shell-fs, 1));
  letter-spacing: 1.5px;
  color: var(--text-bright);
  font-weight: normal;
  margin: 0;
  text-transform: uppercase;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}
.asm-close-btn {
  background: none;
  border: 1px solid var(--border-mid);
  color: var(--text-bright);
  font-family: 'EB Garamond', Georgia, serif;
  font-size: calc(13px * var(--shell-fs, 1));
  letter-spacing: 1px;
  padding: calc(4px * var(--shell-fs, 1)) calc(14px * var(--shell-fs, 1));
  border-radius: 4px;
  cursor: pointer;
  flex-shrink: 0;
  margin-left: 10px;
}
.asm-close-btn:hover { background: var(--bg-control); }

/* ── Canvas fills sim panel ── */
.asm-canvas {
  display: block;
  width: 100%;
  height: 100%;
  image-rendering: pixelated;
}

/* ── Control panel ── */
.asm-ctrl {
  overflow-y: auto;
  -webkit-overflow-scrolling: touch;
}

/* ── Control sections ── */
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
  height: 5px;
  border-radius: 3px;
  background: rgba(var(--slider-track-rgb), 0.35);
  outline: none;
  cursor: pointer;
}
.applet-shell-slider-row input[type=range]::-webkit-slider-thumb {
  -webkit-appearance: none;
  width: calc(22px * var(--shell-fs, 1));
  height: calc(22px * var(--shell-fs, 1));
  border-radius: 50%;
  background: var(--teal-light);
  cursor: pointer;
}
.applet-shell-slider-row input[type=range]::-moz-range-thumb {
  width: calc(22px * var(--shell-fs, 1));
  height: calc(22px * var(--shell-fs, 1));
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
  gap: calc(10px * var(--shell-fs, 1));
  flex-wrap: wrap;
}
.applet-shell-btn {
  background: var(--bg-control);
  color: var(--text-bright);
  border: 1px solid var(--border-mid);
  padding: calc(10px * var(--shell-fs, 1)) calc(18px * var(--shell-fs, 1));
  cursor: pointer;
  font-family: 'EB Garamond', Georgia, serif;
  font-size: calc(18px * var(--shell-fs, 1));
  letter-spacing: 1px;
  border-radius: 6px;
  transition: background 0.15s;
  touch-action: manipulation;
}
.applet-shell-btn:hover  { background: var(--bg-hover); }
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

  /* ── Per-applet CSS ─────────────────────────────────────────────────────── */
  function injectAppletStyles(id) {
    const styleId = 'asm-styles-' + id;
    if (document.getElementById(styleId)) return;
    const s = document.createElement('style');
    s.id = styleId;
    s.textContent = `
#${id}-asm-overlay {
  display: none;
  position: fixed;
  inset: 0;
  z-index: 900;
  background: var(--bg-black);
  pointer-events: none;
}
#${id}-asm-overlay.${id}-asm-open {
  display: block;
  pointer-events: auto;
}

/* Header — slides from top */
#${id}-asm-header {
  top: 0; left: 0; right: 0;
  height: var(--${id}-asm-hdr-h);
  transform: translateY(calc(-1 * var(--${id}-asm-hdr-h)));
}
#${id}-asm-header.asm-open { transform: translateY(0); }

/* Sim — slides from left */
#${id}-asm-sim {
  top: var(--${id}-asm-hdr-h);
  left: 0;
  right: 0;
  height: var(--${id}-asm-S);
  transform: translateX(-110vw);
}
#${id}-asm-sim.asm-open { transform: translateX(0); }

/* Ctrl — slides from bottom (portrait) or right (landscape) */
#${id}-asm-ctrl {
  top: calc(var(--${id}-asm-hdr-h) + var(--${id}-asm-S));
  left: 0; right: 0;
  height: var(--${id}-asm-ctrl-h);
  transform: translateY(110vh);
}
#${id}-asm-ctrl.asm-open        { transform: translateY(0); }
#${id}-asm-ctrl.asm-landscape   { transform: translateX(110vw); }
#${id}-asm-ctrl.asm-landscape.asm-open { transform: translateX(0); }
    `;
    document.head.appendChild(s);
  }

  /* ── HTML scaffold ──────────────────────────────────────────────────────── */
  function buildScaffold(id, title, ctrlHTML) {
    const div = document.createElement('div');
    div.innerHTML = `
<div id="${id}-asm-overlay">

  <div id="${id}-asm-header" class="asm-panel asm-header">
    <div class="asm-header-inner">
      <h2>${title}</h2>
      <button class="asm-close-btn" data-asm-close="${id}">Close</button>
    </div>
  </div>

  <div id="${id}-asm-sim" class="asm-panel asm-sim">
    <canvas id="${id}-canvas" class="asm-canvas"></canvas>
  </div>

  <div id="${id}-asm-ctrl" class="asm-panel asm-ctrl">
    ${ctrlHTML}
  </div>

</div>
    `.trim();
    return div.firstChild;
  }

  /* ── Layout computation ─────────────────────────────────────────────────── */
  function computeLayout() {
    const vw = window.innerWidth;
    const vh = window.innerHeight;
    const landscape = vw > vh;

    let S, ctrlH, ctrlW;
    if (landscape) {
      /* Landscape: sim square fits available height below header */
      S      = vh - HDR_H;
      ctrlH  = vh - HDR_H;
      ctrlW  = vw - S;
    } else {
      /* Portrait: sim is full width, controls fill remaining height */
      S      = vw;
      ctrlH  = Math.max(80, vh - HDR_H - S);
      ctrlW  = vw;
    }
    const fs = Math.min(FS_MAX, Math.max(FS_MIN, vw / FS_BASE));
    return { S, ctrlH, ctrlW, landscape, fs };
  }

  function applyLayout(id) {
    const { S, ctrlH, ctrlW, landscape, fs } = computeLayout();
    const el = document.getElementById(id + '-asm-overlay');
    el.style.setProperty('--' + id + '-asm-hdr-h',   HDR_H  + 'px');
    el.style.setProperty('--' + id + '-asm-S',        S      + 'px');
    el.style.setProperty('--' + id + '-asm-ctrl-h',   ctrlH  + 'px');
    el.style.setProperty('--' + id + '-asm-ctrl-w',   ctrlW  + 'px');
    el.style.setProperty('--shell-fs',                 fs.toFixed(4));

    /* Swap positioning for landscape vs portrait */
    const simEl  = document.getElementById(id + '-asm-sim');
    const ctrlEl = document.getElementById(id + '-asm-ctrl');
    if (landscape) {
      simEl.style.right   = '';
      simEl.style.width   = S + 'px';
      simEl.style.height  = S + 'px';
      ctrlEl.style.top    = HDR_H + 'px';
      ctrlEl.style.left   = S + 'px';
      ctrlEl.style.right  = '0';
      ctrlEl.style.height = ctrlH + 'px';
      ctrlEl.classList.add('asm-landscape');
    } else {
      simEl.style.right   = '0';
      simEl.style.width   = '';
      simEl.style.height  = S + 'px';
      ctrlEl.style.top    = (HDR_H + S) + 'px';
      ctrlEl.style.left   = '0';
      ctrlEl.style.right  = '0';
      ctrlEl.style.height = ctrlH + 'px';
      ctrlEl.classList.remove('asm-landscape');
    }
    return { S, ctrlW: landscape ? (window.innerWidth - S) : window.innerWidth };
  }

  /* ── AppletShellMobile constructor ──────────────────────────────────────── */
  function AppletShellMobile(cfg) {
    const id       = cfg.id;
    const title    = cfg.title;
    const onOpen   = cfg.onOpen   || function () {};
    const onClose  = cfg.onClose  || function () {};
    const onResize = cfg.onResize || null;
    const ctrlHTML = cfg.ctrlHTML || '';

    injectSharedStyles();
    injectAppletStyles(id);

    const scaffold = buildScaffold(id, title, ctrlHTML);
    document.body.appendChild(scaffold);

    scaffold.querySelector('[data-asm-close]').addEventListener('click', function () {
      self.close();
    });

    const panelIds = [
      id + '-asm-header',
      id + '-asm-sim',
      id + '-asm-ctrl',
    ];

    function layout() {
      const { S, ctrlW } = applyLayout(id);
      const canvas = document.getElementById(id + '-canvas');
      if (canvas) {
        canvas.width  = S;
        canvas.height = S;
      }
      return { S, ctrlW };
    }

    const self = {
      open: function () {
        const { S, ctrlW } = layout();
        const canvas = document.getElementById(id + '-canvas');

        /* Prevent body scroll while applet is open */
        document.body.style.overflow = 'hidden';

        document.getElementById(id + '-asm-overlay').classList.add(id + '-asm-open');
        requestAnimationFrame(function () {
          panelIds.forEach(function (pid) {
            document.getElementById(pid).classList.add('asm-open');
          });
        });

        /* Wire touch events to canvas alongside mouse events */
        if (canvas) {
          canvas.addEventListener('touchstart', function (e) {
            e.preventDefault();
            const t = e.touches[0];
            canvas.dispatchEvent(new MouseEvent('mousedown', { clientX: t.clientX, clientY: t.clientY }));
          }, { passive: false });
          canvas.addEventListener('touchmove', function (e) {
            e.preventDefault();
            const t = e.touches[0];
            canvas.dispatchEvent(new MouseEvent('mousemove', { clientX: t.clientX, clientY: t.clientY }));
          }, { passive: false });
          canvas.addEventListener('touchend', function (e) {
            e.preventDefault();
            canvas.dispatchEvent(new MouseEvent('mouseup', {}));
          }, { passive: false });
        }

        onOpen({ canvas: canvas, S: S, ctrlW: ctrlW });
      },

      close: function () {
        onClose();

        document.body.style.overflow = '';

        panelIds.forEach(function (pid) {
          const el = document.getElementById(pid);
          el.classList.remove('asm-open');
          el.classList.add('asm-closing');
        });

        setTimeout(function () {
          document.getElementById(id + '-asm-overlay').classList.remove(id + '-asm-open');
          panelIds.forEach(function (pid) {
            document.getElementById(pid).classList.remove('asm-closing');
          });
        }, 450);
      },
    };

    /* Resize / rotation handler */
    window.addEventListener('resize', function () {
      const overlay = document.getElementById(id + '-asm-overlay');
      if (!overlay.classList.contains(id + '-asm-open')) return;
      const { S, ctrlW } = layout();
      const canvas = document.getElementById(id + '-canvas');
      if (onResize) onResize({ canvas: canvas, S: S, ctrlW: ctrlW });
    });

    return self;
  }

  window.AppletShellMobile = AppletShellMobile;

  /* ── Dispatcher — runs immediately when this script is parsed ── */
  if (document.documentElement.classList.contains('is-mobile')) {
    window.AppletShell = AppletShellMobile;
  }

})();

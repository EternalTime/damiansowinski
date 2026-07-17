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

/* ── Docs panel (header stays on top; docs slides beneath it) ── */
.applet-shell-panel.applet-shell-header { z-index: 915; }
.applet-shell-docs-clip {
  position: fixed;
  z-index: 912;
  pointer-events: none;
  /* clip at the top edge only — glow may spill on the other three sides */
  clip-path: inset(0px -80px -80px -80px);
}
.applet-shell-docs-clip .applet-shell-panel {
  position: absolute;
  left: 0;
  top: 0;
  width: 100%;
  height: 100%;
}
.applet-shell-panel.applet-shell-docs {
  background: var(--bg-void);   /* match the sim canvas, not the ctrl panel */
  border-color: var(--amber);
  box-shadow: 0 12px 35px rgba(var(--amber-rgb), 0.45),
              0 -8px 20px rgba(var(--amber-rgb), 0.25);
  overflow-y: auto;
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
.applet-shell-header-actions {
  display: flex;
  align-items: center;
  gap: calc(6px * var(--shell-fs, 1));
}
.applet-shell-hdr-extra {
  display: flex;
  align-items: center;
  gap: calc(6px * var(--shell-fs, 1));
}
.applet-shell-hdr-extra:empty { display: none; }
.applet-shell-close-btn,
.applet-shell-header-btn {
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
.applet-shell-close-btn:hover,
.applet-shell-header-btn:hover { background: var(--bg-control); }
.applet-shell-header-btn.active {
  background: var(--bg-active-teal);
  border-color: var(--teal-light);
  color: var(--teal-light);
}

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
/* ── Slider tick marks ── */
.applet-shell-slider-wrap {
  position: relative;
  flex: 1;
  min-width: 0;
  display: flex;
  align-items: center;
}
.applet-shell-slider-wrap input[type=range] {
  width: 100%;
}
.applet-shell-tick {
  position: absolute;
  bottom: -6px;
  width: 2px;
  height: 5px;
  background: var(--text-dim);
  transform: translateX(-50%);
  pointer-events: none;
  border-radius: 1px;
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

/* Docs panel — slides down from under the header */
#${id}-docs-clip {
  left:   calc(var(--${id}-left) + 0.05 * var(--${id}-W));
  top:    var(--${id}-top-body);
  width:  calc(0.9 * var(--${id}-W));
  height: var(--${id}-H-body);
}
#${id}-docs-panel {
  transform: translateY(calc(-100% - 80px));
  border-top-left-radius: 0;
  border-top-right-radius: 0;
}
#${id}-docs-panel.applet-shell-open { transform: translateY(0); }

/* Docs mode — applet-specific header buttons hidden */
#${id}-overlay.${id}-docs .applet-shell-hdr-extra { display: none; }

/* Docs mode — sim slides left, ctrl slides right, both frozen */
#${id}-overlay.${id}-docs #${id}-sim-panel.applet-shell-open {
  transform: translateX(calc(0.05 * var(--${id}-W) - var(--${id}-W-sim)));
}
#${id}-overlay.${id}-docs #${id}-ctrl-panel.applet-shell-open {
  transform: translateX(calc(0.95 * var(--${id}-W) - var(--${id}-W-sim)));
}
#${id}-overlay.${id}-docs #${id}-sim-panel,
#${id}-overlay.${id}-docs #${id}-ctrl-panel { pointer-events: none; }
    `;
    document.head.appendChild(s);
  }

  function injectAppletStylesStacked(id) {
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
  top:    var(--${id}-top-sim);
  width:  var(--${id}-W);
  height: var(--${id}-H-sim);
  display: flex;
  flex-direction: column;
  transform: translateX(-110vw);
  border-radius: 0;
}
#${id}-sim-panel.applet-shell-open { transform: translateX(0); }

/* Ctrl panel — slides from bottom */
#${id}-ctrl-panel {
  left:   var(--${id}-left);
  top:    var(--${id}-top-ctrl);
  width:  var(--${id}-W);
  height: var(--${id}-H-ctrl);
  display: flex;
  flex-direction: row;
  align-items: center;
  transform: translateY(120px);
  border-top-left-radius: 0;
  border-top-right-radius: 0;
  overflow-x: auto;
  overflow-y: hidden;
}
#${id}-ctrl-panel.applet-shell-open { transform: translateY(0); }

/* Docs panel — slides down from under the header, covers sim + ctrl */
#${id}-docs-clip {
  left:   var(--${id}-left);
  top:    var(--${id}-top-sim);
  width:  var(--${id}-W);
  height: calc(var(--${id}-H-sim) + var(--${id}-H-ctrl));
}
#${id}-docs-panel {
  transform: translateY(calc(-100% - 80px));
  border-radius: 0;
}
#${id}-docs-panel.applet-shell-open { transform: translateY(0); }

/* Docs mode — freeze sim and ctrl, hide applet-specific header buttons */
#${id}-overlay.${id}-docs #${id}-sim-panel,
#${id}-overlay.${id}-docs #${id}-ctrl-panel { pointer-events: none; }
#${id}-overlay.${id}-docs .applet-shell-hdr-extra { display: none; }
    `;
    document.head.appendChild(s);
  }

  /* ── HTML scaffold ──────────────────────────────────────────────────────── */
  function buildScaffold(id, title, ctrlHTML, headerBtns, docsHTML) {
    const div = document.createElement('div');
    div.innerHTML = `
<div id="${id}-overlay">

  <div id="${id}-header" class="applet-shell-panel applet-shell-header">
    <div class="applet-shell-header-inner">
      <h2>${title}</h2>
      <div class="applet-shell-header-actions">
        <span class="applet-shell-hdr-extra">${headerBtns || ''}</span>
        <button class="applet-shell-header-btn" id="${id}-docs-btn" data-shell-docs="${id}">What is this?</button>
        <button class="applet-shell-close-btn" data-shell-close="${id}">Close</button>
      </div>
    </div>
  </div>

  <div id="${id}-sim-panel" class="applet-shell-panel applet-shell-sim">
    <canvas id="${id}-canvas" class="applet-shell-canvas"></canvas>
  </div>

  <div id="${id}-ctrl-panel" class="applet-shell-panel applet-shell-ctrl">
    ${ctrlHTML}
  </div>

  <div id="${id}-docs-clip" class="applet-shell-docs-clip">
    <div id="${id}-docs-panel" class="applet-shell-panel applet-shell-docs">
      ${docsHTML || ''}
    </div>
  </div>

</div>
    `.trim();
    return div.firstChild;
  }

  /* ── Layout computation — 'side' (default) ─────────────────────────────── */
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

  function applyLayoutSide(id, gap) {
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
    return { W: S, H: S };
  }

  /* ── Layout computation — 'stacked' ─────────────────────────────────────── */
  function computeLayoutStacked() {
    const vw = window.innerWidth;
    const vh = window.innerHeight;
    // H/W = 1.1/PHI; fit within viewport with PAD margins
    const wFromW = vw - 2 * PAD;
    const wFromH = (vh - 2 * PAD) * PHI / 1.1;
    const W      = Math.floor(Math.min(wFromW, wFromH));
    const H      = Math.floor(W * 1.1 / PHI);
    const hdrH   = Math.floor(H * 0.1);
    const ctrlH  = Math.floor(H * 0.1);
    const simH   = H - hdrH - ctrlH;
    const left   = Math.floor((vw - W) / 2);
    const top    = Math.floor((vh - H) / 2);
    return { W, H, simH, hdrH, ctrlH, left, top };
  }

  function applyLayoutStacked(id) {
    const { W, H, simH, hdrH, ctrlH, left, top } = computeLayoutStacked();
    const fs = Math.min(FS_MAX, Math.max(FS_MIN, W / (FS_REF * PHI)));
    const el = document.getElementById(id + '-overlay');
    el.style.setProperty('--' + id + '-left',      left              + 'px');
    el.style.setProperty('--' + id + '-top-hdr',   top               + 'px');
    el.style.setProperty('--' + id + '-top-sim',   (top + hdrH)      + 'px');
    el.style.setProperty('--' + id + '-top-ctrl',  (top + hdrH + simH) + 'px');
    el.style.setProperty('--' + id + '-W',         W                 + 'px');
    el.style.setProperty('--' + id + '-H-hdr',     hdrH              + 'px');
    el.style.setProperty('--' + id + '-H-sim',     simH              + 'px');
    el.style.setProperty('--' + id + '-H-ctrl',    ctrlH             + 'px');
    el.style.setProperty('--shell-fs',              fs.toFixed(4));
    return { W, H: simH };
  }

  function applyLayout(id, gap, layout) {
    if (layout === 'stacked') return applyLayoutStacked(id);
    return applyLayoutSide(id, gap);
  }

  /* ── AppletShell constructor ─────────────────────────────────────────────── */
  function AppletShell(cfg) {
    const id       = cfg.id;
    const title    = cfg.title;
    const gap      = cfg.gap || 0;
    const layoutMode = cfg.layout || 'side';
    const onOpen      = cfg.onOpen   || function () {};
    const onClose     = cfg.onClose  || function () {};
    const onResize    = cfg.onResize || null;
    const onDocsOpen  = cfg.onDocsOpen  || function () {};
    const onDocsClose = cfg.onDocsClose || function () {};
    const ctrlHTML    = cfg.ctrlHTML    || '';
    const headerBtns  = cfg.headerBtns  || '';
    const docsHTML    = cfg.docsHTML    || '';
    const docsSpec    = cfg.docs        || null;

    let docsOpen     = false;
    let docsRendered = false;

    // Inject styles
    injectSharedStyles();
    if (layoutMode === 'stacked') injectAppletStylesStacked(id);
    else injectAppletStyles(id);

    // Build and insert HTML scaffold
    const scaffold = buildScaffold(id, title, ctrlHTML, headerBtns, docsHTML);
    document.body.appendChild(scaffold);

    // Wire up close button
    scaffold.querySelector('[data-shell-close]').addEventListener('click', function () {
      self.close();
    });

    // Wire up docs toggle button
    scaffold.querySelector('[data-shell-docs]').addEventListener('click', function () {
      if (docsOpen) self.closeDocs();
      else          self.openDocs();
      this.blur();
    });

    const panelIds = [id + '-header', id + '-sim-panel', id + '-ctrl-panel'];

    function layout() {
      const { W, H } = applyLayout(id, gap, layoutMode);
      const canvas = document.getElementById(id + '-canvas');
      if (canvas) {
        canvas.width  = W;
        canvas.height = H;
      }
      return { W, H };
    }

    const self = {
      open: function () {
        const { W, H } = layout();
        const canvas   = document.getElementById(id + '-canvas');

        document.getElementById(id + '-overlay').classList.add(id + '-open');
        requestAnimationFrame(function () {
          panelIds.forEach(function (pid) {
            document.getElementById(pid).classList.add('applet-shell-open');
          });
        });

        // Pass S for side layout (W===H===S), or {W,H} for stacked.
        // For backwards compat, S = W for side layout.
        onOpen({ canvas: canvas, S: W, W: W, H: H });
      },

      close: function () {
        if (docsOpen) self.closeDocs(true);
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

      openDocs: function () {
        if (docsOpen) return;
        docsOpen = true;
        if (!docsRendered && docsSpec && window.AppletDocs) {
          window.AppletDocs.render(document.getElementById(id + '-docs-panel'), docsSpec, id);
          docsRendered = true;
        }
        document.getElementById(id + '-overlay').classList.add(id + '-docs');
        document.getElementById(id + '-docs-panel').classList.add('applet-shell-open');
        document.getElementById(id + '-docs-btn').textContent = 'Back to Applet';
        onDocsOpen();
      },

      // skipHook — used when the whole applet is closing; the applet's
      // onClose already stops the sim, so don't fire onDocsClose (resume).
      closeDocs: function (skipHook) {
        if (!docsOpen) return;
        docsOpen = false;
        document.getElementById(id + '-overlay').classList.remove(id + '-docs');
        const dp = document.getElementById(id + '-docs-panel');
        dp.classList.remove('applet-shell-open');
        dp.classList.add('applet-shell-closing');
        setTimeout(function () { dp.classList.remove('applet-shell-closing'); }, 550);
        document.getElementById(id + '-docs-btn').textContent = 'What is this?';
        if (!skipHook) onDocsClose();
      },
    };

    // Resize handler
    window.addEventListener('resize', function () {
      const overlay = document.getElementById(id + '-overlay');
      if (!overlay.classList.contains(id + '-open')) return;
      const { W, H } = layout();
      const canvas   = document.getElementById(id + '-canvas');
      if (onResize) onResize({ canvas: canvas, S: W, W: W, H: H });
    });

    return self;
  }

  window.AppletShellDesktop = AppletShell;
  window.AppletShell = AppletShell;

  /* ═══════════════════════════ AppletDocs ═══════════════════════════════
     Shared docs-panel renderer, used by both the desktop and mobile shells.

     Spec (per applet, passed as cfg.docs):
       {
         whatis:     '…¶…',        // paragraphs split on ¶; $…$ / $$…$$ math
         howto:      '…¶…',
         references: ['bibkey', …] // keys into /assets/data/references.bib
       }
     [bibkey] in prose renders as a superscript [n] linking to the entry.
     Bib parsing/formatting mirrors the MFS site (_layouts/mfs.html).
  ═══════════════════════════════════════════════════════════════════════ */

  const BIB_URL = '/assets/data/references.bib';
  const DOCS_STYLE_ID = 'applet-docs-styles';
  let _bibCache = null;

  function injectDocsStyles() {
    if (document.getElementById(DOCS_STYLE_ID)) return;
    const s = document.createElement('style');
    s.id = DOCS_STYLE_ID;
    s.textContent = `
.adoc-section {
  padding: calc(14px * var(--shell-fs, 1)) calc(20px * var(--shell-fs, 1)) calc(10px * var(--shell-fs, 1));
  border-bottom: 1px solid rgba(var(--amber-rgb), 0.18);
}
.adoc-section:last-child { border-bottom: none; }
.adoc-label {
  font-size: calc(14px * var(--shell-fs, 1));
  text-transform: uppercase;
  letter-spacing: 1.5px;
  color: rgba(var(--amber-rgb), 0.65);
  margin-bottom: calc(8px * var(--shell-fs, 1));
}
.adoc-body p {
  font-size: calc(17px * var(--shell-fs, 1));
  line-height: 1.65;
  color: var(--text-bright);
  margin: 0 0 0.9em;
}
.adoc-body p:last-child { margin-bottom: 0; }
.adoc-math {
  margin: 0.4em 0 1em;
  text-align: center;
  color: var(--text-bright);
  overflow-x: auto;
}
.adoc-figure {
  margin: 0.6em 0 1.1em;
  display: flex;
  justify-content: center;
  align-items: flex-end;
  gap: calc(24px * var(--shell-fs, 1));
  flex-wrap: wrap;
}
.adoc-figure figure { margin: 0; text-align: center; }
.adoc-figure figcaption {
  font-size: calc(14px * var(--shell-fs, 1));
  color: rgba(var(--amber-rgb), 0.65);
  letter-spacing: 0.5px;
  margin-top: 4px;
}
a.adoc-cite, a.adoc-cite:link, a.adoc-cite:visited { color: var(--amber) !important; text-decoration: none; }
a.adoc-cite:hover { color: var(--amber) !important; text-decoration: underline; }
.adoc-reference {
  font-size: calc(15px * var(--shell-fs, 1));
  color: rgba(var(--amber-rgb), 0.55);
  line-height: 1.6;
  margin-bottom: calc(10px * var(--shell-fs, 1));
  word-break: break-word;
}
.adoc-ref-num     { color: var(--amber); margin-right: 0.3em; font-size: 0.85em; }
.adoc-ref-author  { color: rgba(var(--amber-rgb), 0.85); }
.adoc-ref-year    { color: rgba(var(--amber-rgb), 0.45); }
.adoc-ref-title   { font-style: italic; color: rgba(var(--amber-rgb), 0.75); }
.adoc-ref-journal { color: rgba(var(--amber-rgb), 0.55); }
.adoc-ref-volume  { color: rgba(var(--amber-rgb), 0.45); }
a.adoc-ref-link {
  color: rgba(var(--amber-rgb), 0.75) !important;
  text-decoration: underline !important;
  text-decoration-color: rgba(var(--amber-rgb), 0.3) !important;
  transition: color 0.15s, text-shadow 0.15s;
}
a.adoc-ref-link:hover {
  color: var(--amber) !important;
  text-decoration-color: var(--amber) !important;
  text-shadow: 0 0 8px rgba(var(--amber-rgb), 0.8), 0 0 20px rgba(var(--amber-rgb), 0.5);
}
    `;
    document.head.appendChild(s);
  }

  function fetchBib(cb) {
    if (_bibCache) { cb(_bibCache); return; }
    fetch(BIB_URL)
      .then(function (r) { return r.text(); })
      .then(function (text) {
        _bibCache = {};
        const entries = text.split(/(?=@\w+\s*\{)/);
        entries.forEach(function (entry) {
          entry = entry.trim();
          if (!entry) return;
          const keyMatch = entry.match(/^@\w+\s*\{\s*([\w:]+)\s*,/);
          if (keyMatch) _bibCache[keyMatch[1]] = entry;
        });
        cb(_bibCache);
      })
      .catch(function (e) { console.error('Failed to load references.bib', e); cb({}); });
  }

  function unlatex(str) {
    return str
      .replace(/\\["]([aeiouAEIOUy])/g, function (_, c) {
        return { a:'ä',e:'ë',i:'ï',o:'ö',u:'ü',A:'Ä',E:'Ë',I:'Ï',O:'Ö',U:'Ü',y:'ÿ' }[c] || c;
      })
      .replace(/\\'([aeiouAEIOUy])/g, function (_, c) {
        return { a:'á',e:'é',i:'í',o:'ó',u:'ú',A:'Á',E:'É',I:'Í',O:'Ó',U:'Ú',y:'ý' }[c] || c;
      })
      .replace(/\\`([aeiouAEIOU])/g, function (_, c) {
        return { a:'à',e:'è',i:'ì',o:'ò',u:'ù',A:'À',E:'È',I:'Ì',O:'Ò',U:'Ù' }[c] || c;
      })
      .replace(/\\\^([aeiouAEIOU])/g, function (_, c) {
        return { a:'â',e:'ê',i:'î',o:'ô',u:'û',A:'Â',E:'Ê',I:'Î',O:'Ô',U:'Û' }[c] || c;
      })
      .replace(/\\~([nNaAoO])/g, function (_, c) {
        return { n:'ñ',N:'Ñ',a:'ã',A:'Ã',o:'õ',O:'Õ' }[c] || c;
      })
      .replace(/\\c\{?([cCsS])\}?/g, function (_, c) {
        return { c:'ç',C:'Ç',s:'ş',S:'Ş' }[c] || c;
      })
      .replace(/\\ss\b/g, 'ß')
      .replace(/\{([^}]*)\}/g, '$1')
      .trim();
  }

  function parseBibtex(entry) {
    const typeKey = entry.match(/^@(\w+)\s*\{\s*([\w:]+)\s*,/);
    const type = typeKey ? typeKey[1].toLowerCase() : 'misc';
    const fields = {};
    const fieldRe = /(\w+)\s*=\s*(?:\{([^}]*(?:\{[^}]*\}[^}]*)*)\}|"([^"]*)")/g;
    let m;
    while ((m = fieldRe.exec(entry)) !== null) {
      fields[m[1].toLowerCase()] = unlatex(m[2] !== undefined ? m[2] : m[3]);
    }
    return { type: type, fields: fields };
  }

  function formatBibtex(entry) {
    const p = parseBibtex(entry);
    const f = p.fields;
    const author  = f.author  || '';
    const title   = f.title   || '';
    const year    = f.year    || '';
    const journal = f.journal || f.booktitle || f.publisher || '';
    const volume  = f.volume  || '';
    const pages   = f.pages   || '';
    const url     = f.url     || '';
    let html = '';
    if (author)  html += '<span class="adoc-ref-author">' + author + '</span> ';
    if (year)    html += '<span class="adoc-ref-year">(' + year + ')</span>. ';
    if (title)   html += url
      ? '<span class="adoc-ref-title"><a class="adoc-ref-link" href="' + url + '" target="_blank" rel="noopener">' + title + '</a></span>. '
      : '<span class="adoc-ref-title">' + title + '</span>. ';
    if (journal) html += '<span class="adoc-ref-journal">' + journal + '</span>';
    if (volume)  html += ' <span class="adoc-ref-volume">' + volume + '</span>';
    if (pages)   html += ', ' + pages;
    return html || entry;
  }

  function renderParagraphs(txt, refs, id, figures) {
    return txt.split('¶').map(function (p) {
      // Figure block: a paragraph of the form "FIG::name" pulls raw SVG/HTML
      // from the spec's figures map and renders it centered.
      if (p.trim().indexOf('FIG::') === 0) {
        const name = p.trim().slice(5).trim();
        const fig  = figures && figures[name];
        if (!fig) return '';
        return '<div class="adoc-figure">' + fig + '</div>';
      }
      const parts = p.trim().split(/(\$\$[\s\S]*?\$\$)/g);
      return parts.map(function (part) {
        if (/^\$\$[\s\S]*\$\$$/.test(part)) {
          return '<div class="adoc-math">' + part + '</div>';
        }
        if (!part.trim()) return '';
        const text = part.replace(/\s*\[([^\]]+)\]/g, function (match, inner) {
          const keys = inner.split(',').map(function (k) { return k.trim(); });
          const sups = keys.map(function (key) {
            const idx = refs.indexOf(key);
            if (idx === -1) return '';
            return '<a class="adoc-cite" href="#' + id + '-ref-' + key + '">[' + (idx + 1) + ']</a>';
          }).filter(Boolean).join('');
          return sups ? '<sup>' + sups + '</sup>' : match;  // non-bibkey brackets pass through
        });
        return '<p>' + text + '</p>';
      }).join('');
    }).join('');
  }

  function renderDocs(panel, spec, id) {
    injectDocsStyles();
    fetchBib(function (bib) {
      const refs = spec.references || [];
      let html = '';
      if (spec.whatis) {
        html += '<div class="adoc-section"><div class="adoc-label">What is this?</div>' +
                '<div class="adoc-body">' + renderParagraphs(spec.whatis, refs, id, spec.figures) + '</div></div>';
      }
      if (spec.howto) {
        html += '<div class="adoc-section"><div class="adoc-label">How to use</div>' +
                '<div class="adoc-body">' + renderParagraphs(spec.howto, refs, id, spec.figures) + '</div></div>';
      }
      if (refs.length) {
        html += '<div class="adoc-section"><div class="adoc-label">References</div>' +
          refs.map(function (key, i) {
            const entry = bib[key];
            if (!entry) return '<div class="adoc-reference" id="' + id + '-ref-' + key + '">' + key + '</div>';
            return '<div class="adoc-reference" id="' + id + '-ref-' + key + '">' +
                   '<span class="adoc-ref-num">[' + (i + 1) + ']</span> ' + formatBibtex(entry) + '</div>';
          }).join('') + '</div>';
      }
      panel.innerHTML = html;
      if (window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise([panel]);
    });
  }

  window.AppletDocs = { render: renderDocs };

})();

(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  const PINK_DARK   = _c('--pink-dark');
  const PINK_LIGHT  = _c('--pink-light');
  const TEAL_DARK   = _c('--teal-dark');
  const TEAL_LIGHT  = _c('--teal-light');

  /* ── State ── */
  let beta1 = 0.0;   // frame 1 (pink)
  let beta2 = 0.0;   // frame 2 (teal)
  let canvas, ctx, S;

  /* ── Physics helpers ── */
  // For a frame moving at beta, the ct' axis tilts by arctan(beta) from vertical
  // and the x' axis tilts by arctan(beta) from horizontal — both toward the light cone.
  // Grid lines of constant x' are parallel to ct', lines of constant t' are parallel to x'.

  function basisVectors(beta) {
    // Returns unit vectors (in diagram coords) for the ct' and x' axes.
    // ct' direction: (beta, 1) normalised in diagram space (but we don't normalise — we scale by grid spacing)
    // x'  direction: (1, beta) normalised similarly
    return {
      ect: { x: beta, y: 1  },   // direction of ct' axis in (x, ct) diagram coords
      ex:  { x: 1,    y: beta }  // direction of x'  axis
    };
  }

  /* ── Rendering ── */
  function render() {
    ctx.clearRect(0, 0, S, S);

    // Background
    ctx.fillStyle = _c('--bg-black');
    ctx.fillRect(0, 0, S, S);

    const cx = S / 2;   // diagram origin in canvas pixels
    const cy = S / 2;
    const scale = S * 0.42;  // how many pixels per unit (ct or x)

    // Helper: diagram coords (x, ct) → canvas pixels
    function toCanvas(x, ct) {
      return { px: cx + x * scale, py: cy - ct * scale };
    }

    // Helper: draw a line from diagram point a to diagram point b
    function diagLine(ax, act, bx, bct) {
      const a = toCanvas(ax, act);
      const b = toCanvas(bx, bct);
      ctx.beginPath();
      ctx.moveTo(a.px, a.py);
      ctx.lineTo(b.px, b.py);
      ctx.stroke();
    }

    // ── Light cone ──
    ctx.save();

    // Shade the four quadrants: timelike (future/past) brighter, spacelike darker
    const reach = S;  // large enough to fill canvas corners
    const future = toCanvas( 0,  0);
    const tl     = toCanvas(-1.5,  1.5);
    const tr     = toCanvas( 1.5,  1.5);
    const bl     = toCanvas(-1.5, -1.5);
    const br     = toCanvas( 1.5, -1.5);
    const top    = toCanvas( 0,  1.5);
    const bot    = toCanvas( 0, -1.5);
    const left   = toCanvas(-1.5, 0);
    const right  = toCanvas( 1.5, 0);
    const orig   = toCanvas( 0,  0);

    // Future lightcone (upper triangle)
    ctx.beginPath();
    ctx.moveTo(orig.px, orig.py);
    ctx.lineTo(tl.px, tl.py);
    ctx.lineTo(top.px, top.py);
    ctx.lineTo(tr.px, tr.py);
    ctx.closePath();
    ctx.fillStyle = _rgba('--gold', 0.1);
    ctx.fill();

    // Past lightcone (lower triangle)
    ctx.beginPath();
    ctx.moveTo(orig.px, orig.py);
    ctx.lineTo(bl.px, bl.py);
    ctx.lineTo(bot.px, bot.py);
    ctx.lineTo(br.px, br.py);
    ctx.closePath();
    ctx.fillStyle = _rgba('--gold', 0.1);
    ctx.fill();

    // Spacelike regions (left and right triangles) — dimmer
    ctx.beginPath();
    ctx.moveTo(orig.px, orig.py);
    ctx.lineTo(tl.px, tl.py);
    ctx.lineTo(left.px, left.py);
    ctx.lineTo(bl.px, bl.py);
    ctx.closePath();
    ctx.fillStyle = _rgba('--gold', 0.0);
    ctx.fill();

    ctx.beginPath();
    ctx.moveTo(orig.px, orig.py);
    ctx.lineTo(tr.px, tr.py);
    ctx.lineTo(right.px, right.py);
    ctx.lineTo(br.px, br.py);
    ctx.closePath();
    ctx.fillStyle = _rgba('--gold', 0.0);
    ctx.fill();

    // Light cone lines
    ctx.strokeStyle = _rgba('--gold', 1.0);
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 6]);
    diagLine(-1.5, -1.5,  1.5,  1.5);
    diagLine( 1.5, -1.5, -1.5,  1.5);
    ctx.setLineDash([]);
    ctx.restore();

    // ── Draw one frame's grid + worldline ──
    function drawFrame(beta, gridColor, lineColor) {
      const { ect, ex } = basisVectors(beta);
      const N      = 12;      // coarse grid lines each side
      const sub    = 8;     // fine subdivisions per coarse cell
      const Nfine  = N * sub;
      const reach  = 2.5;

      ctx.save();
      ctx.strokeStyle = gridColor;

      // Fine sub-grid (skip multiples of sub — those are the coarse lines)
      ctx.globalAlpha = 0.85;
      ctx.lineWidth = 0.9;
      for (let n = -Nfine; n <= Nfine; n++) {
        if (n % sub === 0) continue;
        const step = n / sub;
        const ox  = step * ex.x,  oct  = step * ex.y;
        diagLine(ox - reach * ect.x, oct - reach * ect.y,
                 ox + reach * ect.x, oct + reach * ect.y);
        const ox2 = step * ect.x, oct2 = step * ect.y;
        diagLine(ox2 - reach * ex.x, oct2 - reach * ex.y,
                 ox2 + reach * ex.x, oct2 + reach * ex.y);
      }

      // Coarse grid
      ctx.globalAlpha = 0.85;
      ctx.lineWidth = 0.8;
      for (let n = -N; n <= N; n++) {
        const ox  = n * ex.x,  oct  = n * ex.y;
        diagLine(ox - reach * ect.x, oct - reach * ect.y,
                 ox + reach * ect.x, oct + reach * ect.y);
        const ox2 = n * ect.x, oct2 = n * ect.y;
        diagLine(ox2 - reach * ex.x, oct2 - reach * ex.y,
                 ox2 + reach * ex.x, oct2 + reach * ex.y);
      }

      ctx.globalAlpha = 1.0;

      // Worldline: x' = 0, i.e. the ct' axis — thicker, lighter colour
      ctx.strokeStyle = lineColor;
      ctx.lineWidth = 3;
      diagLine(-reach * ect.x, -reach * ect.y,
                reach * ect.x,  reach * ect.y);

      ctx.restore();
    }

    drawFrame(beta1, PINK_DARK, PINK_LIGHT);
    drawFrame(beta2, TEAL_DARK, TEAL_LIGHT);
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'ltz',
    title: 'Lorentz Transformation &mdash; Worldlines',
    gap:   0,

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Worldline A &nbsp;<span style="color:var(--pink-dark)">&#9632;</span></div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">&minus;1</span>
          <input type="range" id="ltz-beta1" min="-0.99" max="0.99" step="0.01" value="0">
          <span class="applet-shell-side">&plus;1</span>
          <span class="applet-shell-val" id="ltz-beta1-val">0.00</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Worldline B &nbsp;<span style="color:var(--teal-dark)">&#9632;</span></div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">&minus;1</span>
          <input type="range" id="ltz-beta2" min="-0.99" max="0.99" step="0.01" value="0">
          <span class="applet-shell-side">&plus;1</span>
          <span class="applet-shell-val" id="ltz-beta2-val">0.00</span>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, S: s }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      S      = s;
      render();
    },

    onClose: function () {
      canvas = null;
      ctx    = null;
    },

    onResize: function ({ canvas: c, S: s }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      S      = s;
      render();
    },
  });

  window.ltzOpen  = () => shell.open();
  window.ltzClose = () => shell.close();

  document.getElementById('ltz-beta1').addEventListener('input', function () {
    beta1 = parseFloat(this.value);
    document.getElementById('ltz-beta1-val').textContent = beta1.toFixed(2);
    if (canvas) render();
  });

  document.getElementById('ltz-beta2').addEventListener('input', function () {
    beta2 = parseFloat(this.value);
    document.getElementById('ltz-beta2-val').textContent = beta2.toFixed(2);
    if (canvas) render();
  });

})();

#!/usr/bin/env node
/* Turn a figure of every kind on the spacetimes page by hand, in a browser, as a reader does,
   and hold what the page does to what _tools/README.md, "Turning the figure", says. It needs a
   served copy of the site and, for Chrome, nothing else; _tools/chrome.mjs starts Chrome:

     bundle exec jekyll serve
     node _tools/turn_drag.mjs http://127.0.0.1:4000

   --phone lays the page out as an iPhone held upright, 390 by 844, with touch; --width and
   --height set another screen; and --text <px> sets the browser's default font size, 16 by
   default.

   The figures are Schwarzschild's embedding diagram, a surface of revolution, the Krasnikov
   tube's, a height over a plane, and Gödel's light cones about the axis, a figure of light
   cones. For each it checks that:

   - a mouse drag a fifth of the drawing's width across turns it: the drawing changes, the reset
     button shows in the top right corner of the frame, exactly as tall as the drawing's margin
     above its box at the usual text size, every point drawn stays inside the box, every label
     shown stays inside the drawing, and nothing in the figure glows;
   - a drag a whole width up tilts it to look straight up its axis, and a drag further up
     changes nothing;
   - a double click brings back the published drawing, path for path and label for label, and
     hides the button, and so do Escape after the left arrow has turned it and the reset button;
   - on a phone, a finger drawn across turns it, one drawn up scrolls the page and leaves the
     drawing as it is, and a double tap brings back the published drawing;
   - turned, the print copy the print button sets holds the published drawing, with no reset
     button in it.

   It prints one line per check and exits non-zero if any fails or the page raised any error.
   The checks are run() and the measures in the page PAGE, both exported, so a driver for
   another browser, which gives run() the same few actions, runs the same checks. */
import { pathToFileURL } from 'node:url';

export const FIGURES = [
  { name: 'surface of revolution', metric: 'schwarzschild', kind: 'surfaces' },
  { name: 'height over a plane', metric: 'krasnikov', kind: 'surfaces' },
  { name: 'light cones', metric: 'godel', kind: 'cones' },
];

/* What the page is asked, installed before its own scripts run: every error it raises, and for
   the drawing that turns, what it shows, where, and whether all of it stays in its place. */
export const PAGE = `
  window.__turnErrors = [];
  window.addEventListener('error', function (e) {
    window.__turnErrors.push(e.message + (e.filename ? ' (' + e.filename + ':' + e.lineno + ')' : ''));
  });
  window.addEventListener('unhandledrejection', function (e) {
    window.__turnErrors.push('unhandled rejection: ' + (e.reason && e.reason.stack || e.reason));
  });
  window.__turn = {
    // The drawing of the figure that turns, in the panel or in the print copy.
    drawing: function (metric, kind, root) {
      return (root || document.getElementById('mfs-content-panel'))
        .querySelector('.mfs-cd-drawing[data-turn="' + kind + '"][data-turn-metric="' + metric + '"]');
    },
    // The part of the drawing a reader can reach, on the page's own coordinates: the drawing
    // cut to the window and to every box around it that clips what it holds, such as its own
    // frame, which scrolls sideways, and the panel, which scrolls down and at a large text
    // size starts halfway down the window, below the page's title.
    reach: function (d) {
      var box = { left: 0, top: 0, right: innerWidth, bottom: innerHeight };
      for (var e = d; e && e !== document.documentElement; e = e.parentElement) {
        var cs = getComputedStyle(e);
        if (e === d || cs.overflow !== 'visible') {
          var r = e.getBoundingClientRect();
          box = { left: Math.max(box.left, r.left), top: Math.max(box.top, r.top),
                  right: Math.min(box.right, r.right), bottom: Math.min(box.bottom, r.bottom) };
        }
        // Nothing around a box fixed to the window, as the panel is, clips it.
        if (cs.position === 'fixed') break;
      }
      return { x: (box.left + box.right) / 2, y: (box.top + box.bottom) / 2, width: d.clientWidth,
               left: box.left, right: box.right, top: box.top, bottom: box.bottom };
    },
    state: function (d) {
      var svg = d.querySelector('svg'), figure = d.closest('.mfs-cd-figure'), reset = figure.querySelector('.mfs-turn-reset');
      var box = svg.viewBox.baseVal, margin = 34, out = 0;
      // Every coordinate of every path and point, against the box inside the drawing's margin,
      // in the SVG's own units, which the page writes to a tenth.
      [].forEach.call(svg.querySelectorAll('path'), function (p) {
        var n = (p.getAttribute('d').match(/-?[\\d.]+/g) || []).map(Number);
        for (var i = 0; i + 1 < n.length; i += 2) {
          out = Math.max(out, margin - n[i], n[i] - (box.width - margin), margin - n[i + 1], n[i + 1] - (box.height - margin));
        }
      });
      [].forEach.call(svg.querySelectorAll('circle'), function (c) {
        var x = +c.getAttribute('cx'), y = +c.getAttribute('cy');
        out = Math.max(out, margin - x, x - (box.width - margin), margin - y, y - (box.height - margin));
      });
      var r = d.getBoundingClientRect(), strays = 0;
      var labels = [].map.call(d.querySelectorAll('.mfs-cd-labels .cd-at:not(.mfs-slice-label)'), function (span) {
        var q = span.getBoundingClientRect(), shown = getComputedStyle(span).visibility !== 'hidden';
        if (shown) strays = Math.max(strays, r.left - q.left, q.right - r.right, r.top - q.top, q.bottom - r.bottom);
        return [span.style.left, span.style.top, shown];
      });
      // The reset button is as tall as the drawing's margin at the usual text size, the frame's
      // width over the drawing's units, which the drawing never falls below.
      var b = reset.getBoundingClientRect(), f = figure.getBoundingClientRect(), w = d.parentElement.getBoundingClientRect().width;
      var glow = [figure].concat([].slice.call(figure.querySelectorAll('*'))).filter(function (el) {
        var cs = getComputedStyle(el);
        return (cs.textShadow && cs.textShadow !== 'none') || (cs.boxShadow && cs.boxShadow !== 'none') || /blur|drop-shadow/.test(cs.filter || '');
      }).length;
      return { svg: svg.innerHTML, labels: JSON.stringify(labels), outside: out, strays: strays, glow: glow,
               reset: !reset.hidden && getComputedStyle(reset).display !== 'none',
               corner: reset.hidden ? null : [b.top - r.top, b.right - f.right, b.height - margin * w / box.width,
                                              margin * r.width / box.width - b.height],
               scroll: (document.scrollingElement.scrollTop || 0) + document.getElementById('mfs-content-panel').scrollTop };
    }
  };`;

/* The checks, given a driver: open(metric) shows the spacetime; show(figure) brings its
   drawing into view; evaluate(expression); mouse(x, y, dx, dy), a drag with the left button;
   dblclick(x, y); click(x, y); key(name); and on a phone touch(x, y, dx, dy), a drag of one
   finger, and tap(x, y). report(ok, text) keeps each result. */
export async function run(driver, figure, report, phone) {
  const { evaluate } = driver;
  const q = `window.__turn.drawing(${JSON.stringify(figure.metric)}, ${JSON.stringify(figure.kind)})`;
  const state = () => evaluate(`window.__turn.state(${q})`);
  const reach = async () => { await driver.show(figure); return evaluate(`window.__turn.reach(${q})`); };
  const settle = () => evaluate('new Promise(function (r) { requestAnimationFrame(function () { requestAnimationFrame(r); }); })');
  const at = `${figure.metric}, ${figure.name}`;
  /* A drag of the whole (dx, dy) by `hand`, 'mouse' or 'touch', made as a reader makes a long
     one, in as many strokes as it takes, each within the part of the drawing in reach, which
     at a large text size is smaller than the drag: each stroke turns the figure on from where
     the last one left it. */
  async function drag(hand, dx, dy) {
    for (let left = dx, up = dy, strokes = 0; Math.abs(left) > 0.5 || Math.abs(up) > 0.5; strokes++) {
      const R = await reach(), wide = R.right - R.left - 12, tall = R.bottom - R.top - 12;
      if (wide < 24 || tall < 24 || strokes >= 20) {
        report(false, `${at}: a drag of ${dx.toFixed(0)}, ${dy.toFixed(0)} does not fit what is in reach, ` +
               `${(R.right - R.left).toFixed(0)} by ${(R.bottom - R.top).toFixed(0)}px`);
        return;
      }
      const sx = Math.sign(left) * Math.min(Math.abs(left), wide), sy = Math.sign(up) * Math.min(Math.abs(up), tall);
      await driver[hand](R.x - sx / 2, R.y - sy / 2, sx, sy);
      await settle();
      left -= sx;
      up -= sy;
    }
  }

  await driver.open(figure);
  if (!(await evaluate(`!!${q}`))) { report(false, `${at}: no drawing that turns`); return; }
  let R = await reach();
  const published = await state();
  report(!published.reset, `${at}: the reset button is hidden at the start`);

  function held(s, what) {
    report(s.glow === 0, `${at}: ${what}, nothing in the figure glows (${s.glow} elements)`);
    report(s.outside <= 0.06, `${at}: ${what}, every point drawn stays inside the box (${s.outside.toFixed(2)} out)`);
    report(s.strays <= 0.5, `${at}: ${what}, every label shown stays inside the drawing (${s.strays.toFixed(2)}px out)`);
  }
  held(published, 'at the start');
  function home(s, what) {
    report(s.svg === published.svg && s.labels === published.labels && !s.reset,
           `${at}: ${what} brings back the published drawing and hides the reset button`);
  }

  await drag('mouse', R.width / 5, 0);
  let s = await state();
  report(s.svg !== published.svg, `${at}: a mouse drag across turns the drawing`);
  report(s.reset, `${at}: the reset button shows once it is turned`);
  if (s.corner) {
    report(Math.abs(s.corner[0]) <= 0.5 && Math.abs(s.corner[1]) <= 0.5 && Math.abs(s.corner[2]) <= 0.5 && s.corner[3] >= -0.5,
           `${at}: the reset button stands in the top right corner, as tall as the margin above the box at the usual text size and within the margin (${s.corner.map(v => v.toFixed(2)).join(', ')})`);
  }
  held(s, 'turned');

  await drag('mouse', 0, -R.width);
  const below = await state();
  report(below.svg !== s.svg, `${at}: a drag up tilts it`);
  held(below, 'seen from straight below');
  await drag('mouse', 0, -R.width / 4);
  report((await state()).svg === below.svg, `${at}: a drag further up leaves it seen from straight below`);

  R = await reach();
  await driver.dblclick(R.x, R.y);
  await settle();
  home(await state(), 'a double click');

  R = await reach();
  await driver.click(R.x, R.y);
  await driver.key('ArrowLeft');
  await settle();
  report((await state()).svg !== published.svg, `${at}: the left arrow turns it`);
  await driver.key('Escape');
  await settle();
  home(await state(), 'Escape');

  await drag('mouse', -R.width / 7, R.width / 9);
  // The reader scrolls to the button where the drawing is taller than the screen, as at a large
  // text size, and presses it where nothing covers it.
  const button = await evaluate(`(function (b) {
    b.scrollIntoView({ block: 'nearest', inline: 'nearest' });
    return new Promise(function (r) { requestAnimationFrame(function () { requestAnimationFrame(function () {
      var q = b.getBoundingClientRect(), x = (q.left + q.right) / 2, y = (q.top + q.bottom) / 2;
      r({ x: x, y: y, top: document.elementFromPoint(x, y) === b });
    }); }); });
  })(${q}.closest('.mfs-cd-figure').querySelector('.mfs-turn-reset'))`);
  report(button.top, `${at}: nothing covers the reset button`);
  await driver.click(button.x, button.y);
  await settle();
  home(await state(), 'the reset button');

  if (phone) {
    await drag('touch', R.width / 5, 0);
    s = await state();
    report(s.svg !== published.svg, `${at}: a finger drawn across turns it`);
    held(s, 'turned by a finger');
    R = await reach();
    await driver.tap(R.x, R.y);
    await driver.tap(R.x, R.y);
    await settle();
    home(await state(), 'a double tap');
    R = await reach();
    const before = await state();
    await driver.touch(R.x, Math.min(R.bottom - 10, R.y + 60), 0, -120);
    await settle();
    s = await state();
    report(s.scroll !== before.scroll && s.svg === published.svg,
           `${at}: a finger drawn up scrolls the page (${(s.scroll - before.scroll).toFixed(0)}px) and leaves the drawing as it is`);
  }

  // Turned, the print copy still holds the published drawing.
  await drag('mouse', R.width / 5, 0);
  const printed = await evaluate(`new Promise(function (resolve) {
    window.print = function () {
      var d = window.__turn.drawing(${JSON.stringify(figure.metric)}, ${JSON.stringify(figure.kind)}, document.getElementById('mfs-print-root'));
      var reset = d && d.closest('.mfs-cd-figure').querySelector('.mfs-turn-reset');
      resolve(d ? { svg: d.querySelector('svg').innerHTML, reset: !!reset && reset.hidden } : null);
    };
    document.getElementById('mfs-print-btn').click();
  })`);
  report(!!printed && printed.svg === published.svg && printed.reset,
         `${at}: turned, the print copy holds the published drawing and a hidden reset button`);
  if (driver.printMedia) {
    await driver.printMedia(true);
    const shown = await evaluate(`getComputedStyle(window.__turn.drawing(${JSON.stringify(figure.metric)}, ${JSON.stringify(figure.kind)}, document.getElementById('mfs-print-root')).closest('.mfs-cd-figure').querySelector('.mfs-turn-reset')).display`);
    report(shown === 'none', `${at}: paper shows no reset button`);
    await driver.printMedia(false);
  }
  R = await reach();
  await driver.dblclick(R.x, R.y);
  await settle();

  for (const e of await evaluate('window.__turnErrors.splice(0)')) report(false, `${at}: page error: ${e}`);
}

/* Open a spacetime in the list and wait until it is drawn, typeset and settled; then, for a
   figure of light cones, choose the chart and the view that show it. */
export function opener(evaluate, sleep) {
  return async function open(figure) {
    for (let i = 0; i < 600; i++) {
      if (await evaluate(`!!document.querySelector('.mfs-result[data-id="${figure.metric}"]')`)) break;
      await sleep(100);
    }
    // Ready once the panel has slid in and takes the pointer again, which it does not while it
    // slides, with its mathematics set and its fonts in.
    const ready = `(function () {
      var panel = document.getElementById('mfs-content-panel'), header = panel.querySelector('.mfs-header');
      return !!header && header.dataset.turnSeen !== '1' && !/\\\\[(\\[]/.test(panel.textContent) &&
        !panel.classList.contains('mfs-fonts-wait') && /^translateX\\(0(px)?\\)$/.test(panel.style.transform) &&
        panel.style.pointerEvents !== 'none' && !panel.style.transition;
    })()`;
    await evaluate(`(function (h) { if (h) h.dataset.turnSeen = '1'; })(document.querySelector('#mfs-content-panel .mfs-header'));
                    document.querySelector('.mfs-result[data-id="${figure.metric}"]').click();`);
    for (let i = 0; i < 1200 && !(await evaluate(ready)); i++) await sleep(50);
    await evaluate('document.fonts.ready.then(function () { return new Promise(function (r) { setTimeout(r, 300); }); })');
    if (figure.kind !== 'cones') return;
    const charts = await evaluate(`document.querySelectorAll('#mfs-content-panel .mfs-charts [data-chart]').length || 1`);
    for (let c = 1; c < charts; c++) {
      if (await evaluate(`!!window.__turn.drawing(${JSON.stringify(figure.metric)}, 'cones')`)) break;
      await evaluate(`document.querySelector('#mfs-content-panel .mfs-charts [data-chart="${c}"]').click()`);
      for (let i = 0; i < 600 && !(await evaluate(`!/\\\\[(\\[]/.test(document.getElementById('mfs-content-panel').textContent)`)); i++) await sleep(50);
      await evaluate('document.fonts.ready.then(function () { return new Promise(function (r) { setTimeout(r, 300); }); })');
    }
  };
}

// Bring the figure's view into sight: its button chosen, and its drawing scrolled to the middle.
export const SHOW = `function (d) {
  var view = d.closest('.mfs-nr-view'), section = d.closest('.mfs-diagram');
  if (view.classList.contains('mfs-nr-hidden')) {
    var views = [].slice.call(section.querySelectorAll('.mfs-nr-view'));
    section.querySelectorAll('.mfs-section-head .mfs-choice')[views.indexOf(view)].click();
  }
  d.scrollIntoView({ block: 'center', inline: 'nearest' });
  return new Promise(function (r) { requestAnimationFrame(function () { requestAnimationFrame(r); }); });
}`;

async function main() {
  const { launch, sleep } = await import('./chrome.mjs');
  const argv = process.argv.slice(2);
  const option = (name, fallback) => { const i = argv.indexOf(name); return i >= 0 ? argv[i + 1] : fallback; };
  const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
  const phone = argv.includes('--phone');
  const width = +option('--width', phone ? 390 : 1440), height = +option('--height', phone ? 844 : 900);
  const text = +option('--text', 16);
  const { send, evaluate, finish, onError } = await launch({ width, height, phone, text });
  const results = [];
  let at = 'the list';
  onError(message => results.push([false, `${at}: ${message}`]));
  await send('Page.addScriptToEvaluateOnNewDocument', { source: PAGE });
  await send('Page.navigate', { url: `${base}/MFS/` });

  const mouse = (type, x, y, extra = {}) => send('Input.dispatchMouseEvent', { type, x, y, button: 'left', ...extra });
  const touches = (type, points) => send('Input.dispatchTouchEvent', { type, touchPoints: points });
  const steps = 12;
  const driver = {
    evaluate,
    open: opener(evaluate, sleep),
    show: figure => evaluate(`(${SHOW})(window.__turn.drawing(${JSON.stringify(figure.metric)}, ${JSON.stringify(figure.kind)}))`),
    async mouse(x, y, dx, dy) {
      await mouse('mouseMoved', x, y, { button: 'none' });
      await mouse('mousePressed', x, y, { clickCount: 1, buttons: 1 });
      for (let i = 1; i <= steps; i++) await mouse('mouseMoved', x + dx * i / steps, y + dy * i / steps, { buttons: 1 });
      await mouse('mouseReleased', x + dx, y + dy, { clickCount: 1 });
    },
    async click(x, y) {
      await mouse('mousePressed', x, y, { clickCount: 1, buttons: 1 });
      await mouse('mouseReleased', x, y, { clickCount: 1 });
    },
    async dblclick(x, y) {
      await driver.click(x, y);
      await mouse('mousePressed', x, y, { clickCount: 2, buttons: 1 });
      await mouse('mouseReleased', x, y, { clickCount: 2 });
    },
    async key(name) {
      const code = { ArrowLeft: 37, Escape: 27 }[name];
      await send('Input.dispatchKeyEvent', { type: 'rawKeyDown', key: name, code: name, windowsVirtualKeyCode: code });
      await send('Input.dispatchKeyEvent', { type: 'keyUp', key: name, code: name, windowsVirtualKeyCode: code });
    },
    async touch(x, y, dx, dy) {
      await touches('touchStart', [{ x, y }]);
      for (let i = 1; i <= steps; i++) { await touches('touchMove', [{ x: x + dx * i / steps, y: y + dy * i / steps }]); await sleep(16); }
      await touches('touchEnd', []);
      await sleep(400);
    },
    async tap(x, y) {
      await touches('touchStart', [{ x, y }]);
      await touches('touchEnd', []);
      await sleep(60);
    },
    printMedia: on => send('Emulation.setEmulatedMedia', { media: on ? 'print' : '' }),
  };
  for (const figure of FIGURES) {
    at = figure.metric;
    await run(driver, figure, (ok, line) => results.push([ok, line]), phone);
  }
  let failed = 0;
  for (const [ok, line] of results) {
    if (!ok) failed++;
    console.log(`${ok ? 'ok    ' : 'FAILED'} ${line}`);
  }
  console.log(`${results.length} checks at ${width} by ${height}${phone ? ' as a phone' : ''}` +
    `${text !== 16 ? ` with ${text}px text` : ''}: ${failed ? `${failed} failed` : 'all passed'}`);
  finish(failed ? 1 : 0);
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) await main();

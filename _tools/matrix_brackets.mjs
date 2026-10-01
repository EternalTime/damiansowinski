#!/usr/bin/env node
/* Open every chart of every spacetime on the spacetimes page in headless Chrome, as a reader
   does, and measure how far each matrix's brackets reach past its entries: from the highest
   ink of any entry up to the head of the right bracket, and from the lowest ink of any entry
   down to its foot, in pixels. The smaller of the two is the matrix's clearance, and one
   under --least, a fifth of the formula's em by default, is an error: the page gives a quarter,
   and a glyph's box stands up to a fiftieth of an em past the row it is set in.

     bundle exec jekyll serve
     node _tools/matrix_brackets.mjs http://127.0.0.1:4000

   --metric <id> measures one spacetime and is repeatable; --phone lays the page out as an
   iPhone held upright, 390 by 844; --width and --height set another screen; --text <px> sets
   the browser's default font size, 16 by default; --shot <file> with one --metric saves a
   picture of its first chart's inverse metric with the lines around it. _tools/chrome.mjs starts Chrome.

   An entry's ink is read from MathJax's own boxes: each character's box is its glyph's height
   and depth, and a fraction's rule is its own box, so the measure is of what is drawn and not
   of the row's strut. A stretched bracket or rule inside an entry is its own box; the pieces
   it is assembled from are scaled glyphs whose boxes are not where they are drawn. */
import { writeFileSync } from 'node:fs';
import { launch, sleep } from './chrome.mjs';

const argv = process.argv.slice(2);
function option(name, fallback) {
  const i = argv.indexOf(name);
  return i >= 0 ? argv[i + 1] : fallback;
}
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const only = argv.flatMap((a, i) => (a === '--metric' ? [argv[i + 1]] : []));
const phone = argv.includes('--phone');
const width = Number(option('--width', phone ? 390 : 1440));
const height = Number(option('--height', phone ? 844 : 900));
const text = Number(option('--text', 16));
const least = Number(option('--least', 0.2));
const shot = option('--shot', null);
const LIMIT_S = 120;

const { send, evaluate, finish, onError } = await launch({ width, height, phone, text });
const errors = [];
let opened = 'the list';
onError(message => errors.push(`${opened}: ${message}`));

await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsBusy = 0;
  (function wrap() {
    if (!(window.MathJax && MathJax.typesetPromise && MathJax.startup && MathJax.startup.document)) { setTimeout(wrap, 20); return; }
    var typeset = MathJax.typesetPromise;
    MathJax.typesetPromise = function () {
      window.__mfsBusy++;
      function done() { window.__mfsBusy--; }
      return typeset.apply(this, arguments).then(function (v) { done(); return v; }, function (e) { done(); throw e; });
    };
    window.__mfsWrapped = true;
  })();
  /* Every matrix shown: its label, the em of its formula, and how far the right bracket
     reaches above the highest ink of its entries and below the lowest. */
  window.__mfsBrackets = function (panel) {
    return [].map.call(panel.querySelectorAll('mjx-container mjx-mtable'), function (table) {
      var container = table.closest('mjx-container');
      /* The page stands each table between fences in a box a little taller and deeper than
         itself, so the fences are beside that box where there is one. */
      var box = table.closest('mjx-mpadded') || table;
      var fences = [].filter.call(box.parentNode.children, function (n) { return n.tagName === 'MJX-MO'; });
      var right = fences[fences.length - 1];
      if (!right) return { label: container.textContent.slice(0, 12), fault: 'no right bracket' };
      var b = right.getBoundingClientRect(), top = Infinity, bottom = -Infinity;
      [].forEach.call(table.querySelectorAll('mjx-c, mjx-line, mjx-surd, mjx-stretchy-v, mjx-stretchy-h'), function (c) {
        var r = c.getBoundingClientRect();
        if (!r.width && !r.height) return;
        if (c.tagName === 'MJX-C' && c.closest('mjx-stretchy-v, mjx-stretchy-h')) return;
        top = Math.min(top, r.top); bottom = Math.max(bottom, r.bottom);
      });
      var line = container.closest('.mfs-line'), head = line && line.querySelector('.mfs-line-head');
      var math = container.querySelector('mjx-math'), name = math && math.firstElementChild;
      return { label: (head ? head.textContent : '') || (!name ? 'matrix' : name.tagName === 'MJX-MSUB' ? 'metric'
                 : name.tagName === 'MJX-MSUP' ? 'inverse metric' : 'matrix in prose'),
               em: parseFloat(getComputedStyle(table).fontSize),
               above: top - b.top, below: b.bottom - bottom, height: b.height };
    });
  };
` });

await send('Page.navigate', { url: `${base}/MFS/` });
for (let i = 0; ; i++) {
  if (await evaluate(`!!(window.__mfsWrapped && document.querySelector('.mfs-result'))`)) break;
  if (i > 600) { console.error('The page did not load.'); finish(2); }
  await sleep(100);
}

/* Act, then resolve once the panel holds the chart asked for with every formula typeset, no
   typesetting outstanding, the fonts in and two frames passed. */
function shown(action, chart) {
  return evaluate(`new Promise(function (resolve) {
    var panel = document.getElementById('mfs-content-panel');
    var before = panel.querySelector('.mfs-header');
    var start = performance.now();
    ${action}
    (function poll() {
      var on = panel.querySelector('.mfs-charts .mfs-choice-on'), header = panel.querySelector('.mfs-header');
      var ready = header && header !== before && window.__mfsBusy === 0 &&
        (!on || on.dataset.chart === '${chart}') && !/\\\\[(\\[]/.test(panel.textContent);
      if (!ready) {
        if (performance.now() - start > ${LIMIT_S * 1000}) { resolve(null); return; }
        setTimeout(poll, 25); return;
      }
      document.fonts.ready.then(function () {
        requestAnimationFrame(function () { requestAnimationFrame(function () {
          resolve({ chart: on ? on.textContent : '', matrices: window.__mfsBrackets(panel) });
        }); });
      });
    })();
  })`);
}

/* The inverse metric of the first chart, the last matrix shown, brought to the middle of the
   screen once the panel has come to rest, and its line scrolled to its right bracket where it
   is wider than the screen. */
async function picture() {
  await sleep(1500);
  await evaluate(`(function () {
    var tables = document.querySelectorAll('#mfs-content-panel mjx-container mjx-mtable');
    var table = tables[tables.length - 1], panel = document.getElementById('mfs-content-panel');
    table.scrollIntoView({ block: 'center', inline: 'nearest' });
    for (var n = table.parentNode; n && n !== panel; n = n.parentNode) {
      if (n.scrollWidth > n.clientWidth && getComputedStyle(n).overflowX === 'auto') { n.scrollLeft = n.scrollWidth; break; }
    }
  })()`);
  await sleep(500);
  const png = (await send('Page.captureScreenshot', { format: 'png' })).result.data;
  writeFileSync(shot, Buffer.from(png, 'base64'));
}

const ids = await evaluate(`[].map.call(document.querySelectorAll('.mfs-result'), function (e) { return e.dataset.id; })`);
const rows = [];
function take(id, r, c) {
  if (!r) { errors.push(`${id} / chart ${c}: not ready within ${LIMIT_S}s`); return; }
  for (const m of r.matrices) {
    const at = `${id} / ${r.chart || 'chart'} / ${m.label}`;
    if (m.fault) { errors.push(`${at}: ${m.fault}`); continue; }
    const gap = Math.min(m.above, m.below);
    rows.push({ at, gap, ...m });
    if (gap < least * m.em) {
      errors.push(`${at}: the right bracket clears its entries by ${gap.toFixed(2)}px, under ${(least * m.em).toFixed(2)}px`);
    }
  }
}
for (const id of ids) {
  if (only.length && !only.includes(id)) continue;
  opened = id;
  const first = await shown(`document.querySelector('.mfs-result[data-id="${id}"]').click();`, 0);
  take(id, first, 0);
  if (!first) continue;
  if (shot) await picture();
  const charts = await evaluate(`document.querySelectorAll('#mfs-content-panel .mfs-charts .mfs-choice').length || 1`);
  for (let c = 1; c < charts; c++) {
    take(id, await shown(`document.querySelector('#mfs-content-panel .mfs-charts [data-chart="${c}"]').click();`, c), c);
  }
}

rows.sort((a, b) => a.gap - b.gap);
const layout = `${width}x${height}${phone ? ' phone' : ''}, text ${text}px`;
console.log(`${layout}: ${rows.length} matrices`);
for (const r of rows.slice(0, 5)) {
  console.log(`  ${r.gap.toFixed(2)}px (above ${r.above.toFixed(2)}, below ${r.below.toFixed(2)}, em ${r.em.toFixed(1)})  ${r.at}`);
}
if (rows.length) console.log(`smallest gap: ${rows[0].gap.toFixed(2)}px`);
for (const e of errors) console.error(`ERROR ${e}`);
finish(errors.length || !rows.length ? 1 : 0);

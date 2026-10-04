#!/usr/bin/env node
/* Set every label of every drawing on the spacetimes page as the page sets it, measure the box
   MathJax gives it in headless Chrome, and hold the box the page and the generators estimate
   for it to never being smaller: the labels of the conformal diagrams, of the embedding
   diagrams, of the figures in three dimensions and of the slices marked on a drawing.

     bundle exec jekyll serve
     node _tools/label_sizes.mjs http://127.0.0.1:4000

   The estimate is cdLabelSize() in _layouts/mfs.html, the same function as label_size() in
   _tools/derivations/slices.py, which the page keeps every label's margin of 1.5 em by and the
   generators keep labels clear of each other by, so a label wider than its estimate can stand
   nearer the edge of its drawing than the margin, as "$X = X_s$" on Whittaker's sphere did
   until 3 October 2026. Each label is set in the class the drawing gives it at each size the
   caption takes, 15px and 21px on a phone and a desktop at the usual text and three times
   those at the largest, each once typeset in sight and once typeset unseen, as the page sets
   the views it does not show yet, and measured in ems of that size, the largest kept. MathJax
   sets mathematics to the height of the x of the font around it, 109.9% of the label's em in
   Source Code Pro, and where it cannot measure that, in a view not shown, it takes an x of half
   an em, 113.1%, which is the larger.

   It prints every label wider or taller than its estimate and exits non-zero if there is one;
   --all prints every label with its measure and its estimate, widest overshoot first, and
   --worst <n> the n labels that come nearest their estimate. _tools/chrome.mjs starts Chrome. */
import { readFileSync, readdirSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { launch, sleep } from './chrome.mjs';

const argv = process.argv.slice(2);
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const all = argv.includes('--all');
const worst = argv.includes('--worst') ? +argv[argv.indexOf('--worst') + 1] : 0;
const data = join(dirname(fileURLToPath(import.meta.url)), '..', 'MFS', 'assets', 'data');

// Every label once, by its text and the class the page sets it in, with one place it is drawn.
const labels = new Map();
function add(text, cls, where) {
  const key = cls + '\u0000' + text;
  if (!labels.has(key)) labels.set(key, { text, cls, where });
}
function json(folder) {
  return readdirSync(join(data, folder)).filter(f => f.endsWith('.json'))
    .map(f => [f.slice(0, -5), JSON.parse(readFileSync(join(data, folder, f), 'utf8'))]);
}
for (const [id, file] of json('conformal')) {
  for (const view of file.views || []) {
    for (const L of view.labels || []) add(L.text, L.class, `${id} conformal ${view.id}`);
    for (const mark of view.slices || []) if (mark.place) add(mark.label, 'slice', `${id} conformal ${view.id}`);
  }
}
for (const [id, file] of json('embedding')) {
  for (const view of file.views || []) {
    for (const L of view.figure.labels) {
      if (L.frame) for (const frame of view.movie.frames) add(frame.label, L.class, `${id} embedding ${view.id}`);
      else add(L.text, L.class, `${id} embedding ${view.id}`);
    }
  }
}
for (const [id, file] of json('diagrams')) {
  for (const [system, views] of Object.entries(file.systems || {})) {
    for (const view of views) {
      for (const mark of view.slices || []) if (mark.place) add(mark.label, 'slice', `${id} ${system} ${view.id}`);
    }
  }
  for (const [system, figures] of Object.entries(file.projections || {})) {
    for (const figure of figures) {
      for (const L of figure.labels) add(L.text, L.class, `${id} ${system} ${figure.id}`);
      for (const mark of figure.slices || []) if (mark.place) add(mark.label, 'slice', `${id} ${system} ${figure.id}`);
    }
  }
}
const list = [...labels.values()];

const { send, evaluate, finish } = await launch({});
await send('Page.navigate', { url: `${base}/MFS/` });
for (let i = 0; ; i++) {
  if (await evaluate(`!!(window._mfsLabelSize && window.MathJax && MathJax.typesetPromise && MathJax.startup &&
                         MathJax.startup.document && document.querySelector('.mfs-result'))`)) break;
  if (i > 600) { console.error('The page did not load.'); finish(2); }
  await sleep(100);
}

// A slice's label is set as the page sets it over a spacetime diagram's plot, and every other
// label as a conformal diagram's, in the drawing's own family, at each size the caption takes.
const SIZES = [15, 21, 45, 63, 15, 21, 45, 63], UNSEEN = 4;
const measured = JSON.parse(await evaluate(`new Promise(function (resolve) {
  var list = ${JSON.stringify(list)}, sizes = ${JSON.stringify(SIZES)};
  var grounds = sizes.map(function (size, k) {
    var ground = document.createElement('div');
    ground.style.cssText = 'position:absolute;left:0;top:0;width:20000px;height:10px;visibility:hidden;font-size:' + size + 'px';
    if (k >= ${UNSEEN}) ground.style.display = 'none';
    ground.innerHTML = '<div class="mfs-cd-labels"></div>';
    list.forEach(function (L) {
      var el = document.createElement('span');
      el.className = L.cls === 'slice' ? 'mfs-slice-label' : 'cd-at cd-lab-' + L.cls;
      el.textContent = L.text;
      ground.firstChild.appendChild(el);
    });
    document.body.appendChild(ground);
    return ground;
  });
  MathJax.typesetPromise(grounds).then(function () {
    grounds.forEach(function (ground) { ground.style.display = ''; });
    return document.fonts.ready;
  }).then(function () {
    requestAnimationFrame(function () { requestAnimationFrame(function () {
      var out = list.map(function (L, i) {
        var est = window._mfsLabelSize(L.text), m = { w: 0, h: 0, at: 0, ew: est[0], eh: est[1], error: false };
        grounds.forEach(function (ground, k) {
          var el = ground.firstChild.children[i], r = el.getBoundingClientRect(), em = parseFloat(getComputedStyle(el).fontSize);
          if (r.width / em > m.w) { m.w = r.width / em; m.at = sizes[k] + 'px' + (k >= ${UNSEEN} ? ' unseen' : ''); }
          m.h = Math.max(m.h, r.height / em);
          m.error = m.error || !!el.querySelector('mjx-merror');
          var math = el.querySelector('mjx-container');
          if (math) m.scale = Math.max(m.scale || 0, parseFloat(getComputedStyle(math).fontSize) / parseFloat(getComputedStyle(el).fontSize));
        });
        return m;
      });
      grounds.forEach(function (ground) { ground.remove(); });
      resolve(JSON.stringify(out));
    }); });
  });
})`));

const rows = list.map((L, i) => ({ ...L, ...measured[i], over: Math.max(measured[i].w - measured[i].ew, measured[i].h - measured[i].eh) }));
const fmt = r => `${r.over > 0.005 ? 'wider' : 'clear'} ${r.over.toFixed(3)} em  "${r.text}" (${r.cls}, ${r.where}): ` +
  `scale ${(r.scale || 0).toFixed(3)} set ${r.w.toFixed(3)} by ${r.h.toFixed(3)} em (widest at ${r.at}), estimated ${r.ew.toFixed(2)} by ${r.eh.toFixed(2)}`;
const faults = rows.filter(r => r.over > 0.005 || r.error);
rows.sort((a, b) => b.over - a.over);
if (all) rows.forEach(r => console.log(fmt(r)));
else if (worst) rows.slice(0, worst).forEach(r => console.log(fmt(r)));
for (const r of faults) console.log((r.error ? 'MathJax error in ' : '') + fmt(r));
console.log(`${rows.length} labels measured, ${faults.length} larger than their estimate.`);
finish(faults.length ? 1 : 0);

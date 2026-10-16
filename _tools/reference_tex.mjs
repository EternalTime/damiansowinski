#!/usr/bin/env node
/* Open the spacetimes page, the applets page and the publications page in headless Chrome as a
   reader does and hold every reference against the TeX error the Gott-Li spacetime's 13th
   reference, Elliptic de Sitter space: $dS/\mathbb{Z}_2$, showed until 16 October 2026:

     bundle exec jekyll serve
     node _tools/reference_tex.mjs http://127.0.0.1:4000 references.json

   Each reference shows its mathematics set by MathJax, as many formulas as the bibliography
   writes between dollar signs and none of them an error, and outside them no TeX at all: no
   backslash, no dollar sign, no brace, no TeX quote and no TeX dash. Each reference also shows
   the words of its title as the reader writes them out. The file named second says what to
   look for; _tools/test_reference_tex.py writes it from the bibliography and runs this. The
   spacetimes are each opened from the list, the applets' notes are given every entry of the
   bibliography to show, and the publications are shown once as they are and once from a
   bibliography whose first publication carries a title with mathematics, since none of them
   does. Every console error and page error is an error. It exits non-zero on any of them.
   _tools/chrome.mjs starts Chrome. */
import { readFileSync } from 'node:fs';
import { launch, sleep } from './chrome.mjs';

const base = (process.argv[2] || 'http://127.0.0.1:4000').replace(/\/$/, '');
const spec = JSON.parse(readFileSync(process.argv[3], 'utf8'));

const { send, evaluate, finish, onError } = await launch();
const errors = [], failures = [];
onError(message => errors.push(message));
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsErrors = [];
  window.addEventListener('error', function (e) { window.__mfsErrors.push(String(e.message)); });
  window.addEventListener('unhandledrejection', function (e) { window.__mfsErrors.push(String(e.reason)); });
` });

function check(ok, what) {
  if (!ok) { console.log(`FAIL ${what}`); failures.push(what); }
  return ok;
}
async function until(expression, what, seconds = 60) {
  for (let t = 0; t < seconds * 20; t++) {
    try { if (await evaluate(expression)) return true; } catch {}
    await sleep(50);
  }
  return check(false, `${what} within ${seconds}s`);
}
async function open(path) {
  for (const e of (await evaluate('window.__mfsErrors || []').catch(() => []))) errors.push(`page error: ${e}`);
  await send('Page.navigate', { url: base + path });
}

/* What the page shows of each reference, read in the page: the formulas MathJax set and those
   it set as an error, which an undefined command is, in red, and the text with the formulas
   taken out. `expected` maps an element's id to the reference's expectation. */
const READ = `function (expected) {
  var shown = {};
  Object.keys(expected).forEach(function (id) {
    var el = document.getElementById(id);
    if (!el) return;
    var formulas = el.querySelectorAll('mjx-container');
    var wrong = [].filter.call(formulas, function (f) {
      return f.querySelector('mjx-merror, [data-mjx-error]') ||
        [].some.call(f.querySelectorAll('mjx-mtext, mjx-mi, mjx-mo'), function (n) { return /red/.test(n.getAttribute('style') || ''); });
    });
    var copy = el.cloneNode(true);
    [].forEach.call(copy.querySelectorAll('mjx-container'), function (f) { f.remove(); });
    shown[id] = { formulas: formulas.length, wrong: wrong.length, text: copy.textContent.replace(/\\s+/g, ' ').trim() };
  });
  return shown;
}`;
function typeset(expected) {
  return `(function (expected) {
    return Object.keys(expected).every(function (id) {
      var el = document.getElementById(id);
      return el && el.querySelectorAll('mjx-container').length >= expected[id].formulas;
    });
  })(${JSON.stringify(expected)})`;
}
// TeX left in what a reader sees: a command, a dollar sign, a brace, TeX's quotes and dashes.
const TEX = /\\|\$|[{}]|``|''|`|--/;

function hold(where, expected, shown) {
  let ok = true;
  for (const [id, { key, formulas, words }] of Object.entries(expected)) {
    const s = shown[id];
    if (!check(s, `${where} [${key}] is shown`)) { ok = false; continue; }
    ok = check(s.formulas === formulas, `${where} [${key}] sets ${formulas} formulas: it sets ${s.formulas}: ${s.text}`) && ok;
    ok = check(s.wrong === 0, `${where} [${key}] sets ${s.wrong} formulas as errors: ${s.text}`) && ok;
    const tex = s.text.match(TEX);
    ok = check(!tex, `${where} [${key}] shows ${tex && tex[0]} as TeX: ${s.text}`) && ok;
    for (const word of words) ok = check(s.text.indexOf(word) !== -1, `${where} [${key}] shows "${word}": ${s.text}`) && ok;
  }
  return ok;
}

/* The spacetimes, each chosen from the list as a reader does, the Gott-Li spacetime first. The
   last spacetime's references are marked so that its page is never read for the next one's. */
await open('/MFS/');
await until("window._mfsIndex && window._mfsOpen && document.querySelectorAll('#mfs-search-results .mfs-result').length > 0", 'the list comes up');
let shownReferences = 0;
for (const [id, references] of Object.entries(spec.spacetimes)) {
  await evaluate("(function () { var r = document.querySelector('.mfs-references'); if (r) r.dataset.read = 'yes'; })()");
  await evaluate(`window._mfsOpen(${JSON.stringify(id)})`);
  const expected = {};
  for (const [key, e] of Object.entries(references)) expected[`ref-${key}`] = { key, ...e };
  if (!await until(`(function () {
    var r = document.querySelector('#mfs-content-panel .mfs-references');
    return r && !r.dataset.read && ${typeset(expected)};
  })()`, `${id}'s references are set`)) continue;
  const shown = await evaluate(`(${READ})(${JSON.stringify(expected)})`);
  for (const [refId, { key, number }] of Object.entries(expected)) {
    const s = shown[refId];
    if (s) check(s.text.startsWith(`[${number}] `), `${id} lists [${key}] as [${number}]: ${s.text}`);
  }
  if (hold(id, expected, shown)) shownReferences += Object.keys(expected).length;
}
console.log(`ok   ${shownReferences} references in ${Object.keys(spec.spacetimes).length} spacetimes`);

/* The applets' notes, given every entry of the bibliography to show. */
await open('/applets.html');
await until('window.AppletDocs && window.MathJax && MathJax.typesetPromise', 'the applets come up');
const every = {};
spec.applets.forEach((e, i) => { every[`tex-ref-${e.key}`] = e; });
await evaluate(`(function () {
  var panel = document.createElement('div');
  panel.id = 'tex-notes';
  document.body.appendChild(panel);
  AppletDocs.render(panel, { references: ${JSON.stringify(spec.applets.map(e => e.key))} }, 'tex');
})()`);
if (await until(`document.querySelectorAll('#tex-notes .adoc-reference').length === ${spec.applets.length} && ${typeset(every)}`,
                'the notes show every entry of the bibliography', 180)) {
  if (hold('the applets', every, await evaluate(`(${READ})(${JSON.stringify(every)})`))) {
    console.log(`ok   every one of the ${spec.applets.length} entries of the bibliography in the applets' notes`);
  }
}

/* The publications, in the order of _data/publication_keys.json, then again from a
   bibliography of theirs alone whose first publication carries a title with mathematics. */
async function publications(what, list) {
  await until("document.querySelectorAll('#publications-list .csl-entry').length > 0", `${what} come up`);
  const expected = {};
  list.forEach((e, i) => { expected[`tex-publication-${i}`] = e; });
  await evaluate(`[].forEach.call(document.querySelectorAll('#publications-list .csl-entry'), function (e, i) { e.id = 'tex-publication-' + i; })`);
  if (await until(typeset(expected), `${what} are set`)) {
    if (hold(what, expected, await evaluate(`(${READ})(${JSON.stringify(expected)})`))) {
      console.log(`ok   ${what}`);
    }
  }
}
await open('/publications/');
await publications(`the ${spec.publications.length} publications`, spec.publications);
await send('Page.addScriptToEvaluateOnNewDocument', { source: `(function () {
  var fetched = window.fetch, bib = ${JSON.stringify(spec.withMathematics.bib)};
  window.fetch = function (url) {
    if (String(url).indexOf('references.bib') !== -1) return Promise.resolve(new Response(bib));
    return fetched.apply(this, arguments);
  };
})()` });
await open('/publications/');
await publications('the publications with mathematics in a title', spec.withMathematics.publications);

await sleep(300);
for (const e of (await evaluate('window.__mfsErrors || []'))) errors.push(`page error: ${e}`);
for (const e of errors) console.log(`FAIL ${e}`);
finish(failures.length || errors.length ? 1 : 0);

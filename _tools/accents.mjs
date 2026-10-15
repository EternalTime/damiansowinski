#!/usr/bin/env node
/* Open the spacetimes page and the publications page in headless Chrome as a reader does and
   hold their accents to what the captain asked on 15 October 2026:

     bundle exec jekyll serve
     node _tools/accents.mjs http://127.0.0.1:4000 accents.json

   A search typed a key at a time, with or without its accents and in any capitals, lists the
   same spacetimes and offers the same keywords, so "natario" finds Natário, "godel" Gödel,
   "lemaitre" Lemaître and "nordstrom" Nordström, and "kahler" offers the keyword Kähler. Every
   reference a spacetime lists and every publication shows each word written with a TeX accent
   in the bibliography as that letter with that accent, and no accent as TeX. The file named
   second says what to look for; _tools/test_accents.py writes it from the bibliography and the
   metric files and runs this. Every console error and page error is an error. It exits non-zero
   on any of them. _tools/chrome.mjs starts Chrome. */
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
async function until(expression, what, seconds = 30) {
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
// An accent still written as TeX: a sign accent, or a letter accent before a brace.
const TEX_ACCENT = /\\(['"`^~=.]|[uvHckr]\s*\{)/;

/* The search: each spelling typed a key at a time into an emptied field, as a reader types it,
   and the list's spacetimes and the keywords offered read off the page. */
async function search(text) {
  await evaluate("var f = document.getElementById('mfs-search-input'); f.focus(); f.select()");
  await send('Input.dispatchKeyEvent', { type: 'keyDown', key: 'Backspace', code: 'Backspace', windowsVirtualKeyCode: 8 });
  await send('Input.dispatchKeyEvent', { type: 'keyUp', key: 'Backspace', code: 'Backspace', windowsVirtualKeyCode: 8 });
  for (const key of text) await send('Input.insertText', { text: key });
  return evaluate(`({
    listed: [].map.call(document.querySelectorAll('#mfs-search-results .mfs-result'), function (r) { return r.dataset.id; }),
    offered: [].map.call(document.querySelectorAll('#mfs-keywords .mfs-keyword'), function (k) { return k.textContent; })
  })`);
}
await open('/MFS/');
await until("window._mfsIndex && document.querySelectorAll('#mfs-search-results .mfs-result').length > 0", 'the list comes up');
for (const { spellings, finds } of spec.search) {
  const found = [];
  for (const spelling of spellings) found.push(await search(spelling));
  const [first] = found;
  let ok = check(first.listed.indexOf(finds) !== -1, `"${spellings[0]}" lists ${finds}: it lists ${JSON.stringify(first.listed)}`);
  spellings.forEach((spelling, i) => {
    ok = check(JSON.stringify(found[i]) === JSON.stringify(first),
      `"${spelling}" lists and offers what "${spellings[0]}" does: ${JSON.stringify(found[i])} against ${JSON.stringify(first)}`) && ok;
  });
  if (ok) console.log(`ok   ${spellings.map(s => `"${s}"`).join(', ')} each list ${finds} among ${first.listed.length}`);
}
for (const { typed, offers, carrying } of spec.keywords) {
  const found = await search(typed);
  let ok = check(found.offered.indexOf(offers) !== -1, `"${typed}" offers the keyword ${offers}: it offers ${JSON.stringify(found.offered)}`);
  ok = check(JSON.stringify(found.listed) === JSON.stringify(carrying),
    `"${typed}" lists the spacetimes carrying ${offers}: ${JSON.stringify(found.listed)} against ${JSON.stringify(carrying)}`) && ok;
  if (ok) console.log(`ok   "${typed}" offers ${offers} and lists the ${carrying.length} carrying it`);
}

/* A spacetime's references, each spacetime opened at its own address as a reader following a
   link does. */
let words = 0, wrong = 0;
for (const [id, references] of Object.entries(spec.spacetimes)) {
  await open(`/MFS/?spacetime=${encodeURIComponent(id)}`);
  const keys = Object.keys(references);
  if (!await until(`(function () {
    return ${JSON.stringify(keys)}.every(function (key) {
      var r = document.getElementById('ref-' + key);
      return r && !r.classList.contains('mfs-ref-missing');
    });
  })()`, `${id}'s references come up`)) { wrong++; continue; }
  const shown = await evaluate(`(function () {
    var shown = {};
    ${JSON.stringify(keys)}.forEach(function (key) { shown[key] = document.getElementById('ref-' + key).textContent; });
    return shown;
  })()`);
  for (const [key, expected] of Object.entries(references)) {
    if (!check(!TEX_ACCENT.test(shown[key]), `${id} [${key}] shows an accent as TeX: ${shown[key]}`)) wrong++;
    for (const word of expected) {
      words++;
      if (!check(shown[key].indexOf(word) !== -1, `${id} [${key}] shows ${word}: ${shown[key]}`)) wrong++;
    }
  }
}
if (!wrong) console.log(`ok   ${words} accented words in the references of ${Object.keys(spec.spacetimes).length} spacetimes`);

/* The publications, in the order of _data/publication_keys.json. */
await open('/publications/');
if (await until("document.querySelectorAll('#publications-list .csl-entry').length > 0", 'the publications come up')) {
  const shown = await evaluate("[].map.call(document.querySelectorAll('#publications-list .csl-entry'), function (e) { return e.textContent; })");
  check(shown.length === spec.publications.length, `the publications page shows ${spec.publications.length} publications: it shows ${shown.length}`);
  let ok = true;
  spec.publications.forEach((expected, i) => {
    ok = check(!TEX_ACCENT.test(shown[i] || ''), `publication ${i + 1} shows an accent as TeX: ${shown[i]}`) && ok;
    for (const word of expected) ok = check((shown[i] || '').indexOf(word) !== -1, `publication ${i + 1} shows ${word}: ${shown[i]}`) && ok;
  });
  if (ok) console.log(`ok   every accent in the ${shown.length} publications`);
}

await sleep(300);
for (const e of (await evaluate('window.__mfsErrors || []'))) errors.push(`page error: ${e}`);
for (const e of errors) console.log(`FAIL ${e}`);
finish(failures.length || errors.length ? 1 : 0);

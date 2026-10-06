#!/usr/bin/env node
/* Open a spacetime at its own address, the spacetimes page with ?spacetime= and its id, in
   headless Chrome as a reader following a link does, and hold that its contents come up and
   the page raises no error on the way:

     bundle exec jekyll serve
     node _tools/deep_link.mjs http://127.0.0.1:4000 kerr

   On 6 October 2026 such a link threw "Cannot read properties of undefined (reading
   'length')" before the search's keywords were set up, and the page never came up.
   _tools/test_deep_link.py runs it in the suite. It exits non-zero on any page error or console
   error, or when the spacetime's name is not in its panel within 30 seconds.
   _tools/chrome.mjs starts Chrome. */
import { launch, sleep } from './chrome.mjs';

const base = (process.argv[2] || 'http://127.0.0.1:4000').replace(/\/$/, '');
const id = process.argv[3] || 'kerr';

const { send, evaluate, finish, onError } = await launch();
const errors = [];
onError(message => errors.push(message));
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsErrors = [];
  window.addEventListener('error', function (e) { window.__mfsErrors.push(String(e.message)); });
  window.addEventListener('unhandledrejection', function (e) { window.__mfsErrors.push(String(e.reason)); });
` });
await send('Page.navigate', { url: `${base}/MFS/?spacetime=${encodeURIComponent(id)}` });

const shown = `(function () {
  var entry = (window._mfsIndex || []).find(function (m) { return m.id === ${JSON.stringify(id)}; });
  var panel = document.getElementById('mfs-content-panel');
  return !!(entry && panel && panel.textContent.indexOf(entry.name) !== -1);
})()`;
let up = false;
for (let t = 0; t < 600 && !up; t++) {
  try { up = await evaluate(shown); } catch {}
  if (!up) await sleep(50);
}
await sleep(500);
for (const e of (await evaluate('window.__mfsErrors || []'))) errors.push(`page error: ${e}`);

if (!up) console.log(`FAIL ${id} did not come up at its own address within 30s`);
for (const e of errors) console.log(`FAIL ${e}`);
if (up && !errors.length) console.log(`ok   ${id} came up at its own address without an error`);
finish(up && !errors.length ? 0 : 1);

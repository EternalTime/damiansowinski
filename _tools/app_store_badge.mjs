#!/usr/bin/env node
/* Open the spacetimes page in headless Chrome as a reader does and hold Apple's App Store badge
   to what the captain asked on 15 October 2026: seen without hunting, on a desktop in the title's
   band before the exit sign and centred on its line, on a phone in a row of its own under the
   title and the sign, Apple's artwork at its own proportions and 40 pixels tall at the usual
   text size, with no shadow and no filter, linking to the application's page on the App Store:

     bundle exec jekyll serve
     node _tools/app_store_badge.mjs http://127.0.0.1:4000
     node _tools/app_store_badge.mjs http://127.0.0.1:4000 --phone

   Once the panels have come in it is drawn whole, inside the window, clear of the title, the
   sign and the panels, and a press at its middle lands on it. Every console error and page
   error is an error. It exits non-zero on any of them. _tools/test_app_store_badge.py runs it
   at three sizes; _tools/chrome.mjs starts Chrome. */
import { launch, sleep } from './chrome.mjs';

const argv = process.argv.slice(2);
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const phone = argv.includes('--phone');
const size = (argv.find(a => /^--size=\d+x\d+$/.test(a)) || (phone ? '--size=390x844' : '--size=1440x900')).slice(7).split('x').map(Number);
const APP = 'https://apps.apple.com/us/app/my-favorite-spacetimes/id6810529195';

const { send, evaluate, finish, onError } = await launch({ width: size[0], height: size[1], phone });
const errors = [], failures = [];
onError(message => errors.push(message));
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsErrors = [];
  window.addEventListener('error', function (e) { window.__mfsErrors.push(String(e.message)); });
  window.addEventListener('unhandledrejection', function (e) { window.__mfsErrors.push(String(e.reason)); });
` });
function check(ok, what) {
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${size.join(' by ')}: ${what}`);
  if (!ok) failures.push(what);
}

await send('Page.navigate', { url: `${base}/MFS/` });
let shown = false;
for (let t = 0; t < 400 && !shown; t++) {
  try { shown = await evaluate("getComputedStyle(document.getElementById('mfs-appstore')).opacity === '1'"); } catch {}
  if (!shown) await sleep(50);
}
check(shown, 'the badge comes in with the panels');
await sleep(800);

const b = await evaluate(`(function () {
  function box(e) { var r = e.getBoundingClientRect(); return { left: r.left, top: r.top, right: r.right, bottom: r.bottom }; }
  var row = document.getElementById('mfs-appstore'), link = row.querySelector('a'), img = row.querySelector('img');
  var b = box(img), under = document.elementFromPoint((b.left + b.right) / 2, (b.top + b.bottom) / 2);
  var styles = [row, link, img].map(function (e) { var s = getComputedStyle(e); return [s.boxShadow, s.filter, s.textShadow]; });
  return {
    img: b, link: box(link), href: link.href, target: link.target, rel: link.rel, alt: img.alt,
    natural: [img.naturalWidth, img.naturalHeight], src: img.getAttribute('src'),
    pressed: under === img, styles: styles,
    title: box(document.getElementById('mfs-title')), exit: box(document.getElementById('mfs-exit')),
    panels: ['mfs-search-panel', 'mfs-coffee-panel'].map(function (id) { return box(document.getElementById(id)); }),
    view: [window.innerWidth, window.innerHeight]
  };
})()`);
const apart = (p, q) => p.right <= q.left || q.right <= p.left || p.bottom <= q.top || q.bottom <= p.top;
const w = b.img.right - b.img.left, h = b.img.bottom - b.img.top;

check(b.href === APP && b.target === '_blank' && /noopener/.test(b.rel), `it links to ${APP} in a new tab: ${b.href}`);
check(b.alt === 'Download on the App Store' && /app_store_button_white\.svg$/.test(b.src), `it is Apple's white badge: ${b.src}`);
check(b.natural[0] > 0, 'its artwork loaded');
check(Math.abs(h - 40) < 0.5 && Math.abs(w / h - 119.66407 / 40) < 0.01, `it is 40 pixels tall at Apple's proportions: ${w.toFixed(1)} by ${h.toFixed(1)}`);
check(b.styles.every(s => s[0] === 'none' && s[1] === 'none' && s[2] === 'none'), `it has no shadow and no filter: ${JSON.stringify(b.styles)}`);
check(Math.abs((b.link.right - b.link.left) - w) < 0.5 && Math.abs((b.link.bottom - b.link.top) - h) < 0.5, 'a press lands on it only where it is drawn');
check(b.img.left >= 0 && b.img.top >= 0 && b.img.right <= b.view[0] && b.img.bottom <= b.view[1], 'it is inside the window');
check(b.pressed, 'a press at its middle lands on it');
check([b.title, b.exit, ...b.panels].every(r => apart(b.img, r)), 'it is clear of the title, the exit sign and the panels');
if (phone) {
  check(b.img.top >= Math.max(b.title.bottom, b.exit.bottom) && b.img.bottom <= b.panels[0].top,
        'it stands in its own row between the title and the list');
  check(Math.abs(b.img.left - b.title.left) < 1, `its row starts where the title's words do: ${b.img.left} against ${b.title.left}`);
} else {
  check(b.img.right <= b.exit.left && b.img.left >= b.title.right, 'it stands between the title and the exit sign');
  check(Math.abs((b.img.top + b.img.bottom) / 2 - (b.exit.top + b.exit.bottom) / 2) < 1, "it is centred on the exit sign's line");
  check(b.img.bottom <= b.panels[0].top, 'it stands above the panels');
}

for (const e of (await evaluate('window.__mfsErrors || []'))) errors.push(`page error: ${e}`);
for (const e of errors) console.log(`FAIL ${e}`);
finish(failures.length || errors.length ? 1 : 0);

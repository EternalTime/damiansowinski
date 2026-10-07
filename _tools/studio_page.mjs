#!/usr/bin/env node
/* Open the spacetimes page in headless Chrome as a reader does and hold "Studio" and the page
   it opens, "More from Owl's Nest Creations", to what _tools/README.md under "More from Owl's
   Nest Creations" says they do:

     bundle exec jekyll serve
     node _tools/studio_page.mjs http://127.0.0.1:4000
     node _tools/studio_page.mjs http://127.0.0.1:4000 --phone

   "Studio" stands before the exit sign in the sign's own face and size, as far from it as the
   sign stands from the window's edge, on its line, over nothing else and under nothing else; on
   a phone it stands under the sign at the sign's end, and the title keeps its two lines. A
   press on it, with the list showing and with a spacetime open, puts the page in the spacetime's
   panel: the heading, each application the studio site's list of its games shows, its name and
   its line over its icon, drawn, a link to its App Store page in a tab of its own, and under them
   the captain's two lines to the studio site and this one, and nothing else, every word set as
   the coffee panel's paragraph is. The
   list is back, with no spacetime marked open, and the spacetime that was open opens again from
   its name. With --phone the page is laid out as an iPhone held upright, with --side as one on
   its side, --width sets a desktop window's width and --text the browser's text size in pixels; on a phone the press brings the page into view below the list. Every console error and
   page error is an error, and it exits non-zero on any of them. _tools/chrome.mjs starts Chrome. */
import { launch, sleep } from './chrome.mjs';

const argv = process.argv.slice(2);
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const phone = argv.includes('--phone') || argv.includes('--side');
const side = argv.includes('--side');
const text = argv.includes('--text') ? Number(argv[argv.indexOf('--text') + 1]) : 16;
const wide = argv.includes('--width') ? Number(argv[argv.indexOf('--width') + 1]) : 1440;
const OWN_STORE_ID = '6810529195';

const { send, evaluate, finish, onError } = await launch({
  width: side ? 844 : phone ? 390 : wide, height: side ? 390 : phone ? 844 : 900, phone, text
});
const errors = [], failures = [];
onError(message => errors.push(message));
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsErrors = [];
  window.addEventListener('error', function (e) { window.__mfsErrors.push(String(e.message)); });
  window.addEventListener('unhandledrejection', function (e) { window.__mfsErrors.push(String(e.reason)); });
` });

function check(ok, what) {
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${what}`);
  if (!ok) failures.push(what);
}
async function until(expression, what, seconds = 20) {
  for (let t = 0; t < seconds * 20; t++) {
    if (await evaluate(expression)) return true;
    await sleep(50);
  }
  check(false, `${what} within ${seconds}s`);
  return false;
}
const box = id => evaluate(`(function () { var r = document.getElementById(${JSON.stringify(id)}).getBoundingClientRect();
  return { left: r.left, right: r.right, top: r.top, bottom: r.bottom, width: r.width, height: r.height }; })()`);
// A press where a reader's finger or pointer lands, the middle of what is pressed, and only if
// that is what stands there.
async function press(selector, what) {
  // A page brought into view is scrolled there smoothly, and a press lands where things stand
  // once it has come to rest.
  let was = null;
  for (let t = 0; t < 100; t++) {
    const now = await evaluate('window.scrollY');
    if (now === was) break;
    was = now;
    await sleep(150);
  }
  const at = await evaluate(`(function () {
    var el = document.querySelector(${JSON.stringify(selector)});
    if (!el) return null;
    el.scrollIntoView({ block: 'center' });
    var r = el.getBoundingClientRect(), x = r.left + r.width / 2, y = r.top + r.height / 2;
    var hit = document.elementFromPoint(x, y);
    return { x: x, y: y, mine: !!hit && (hit === el || el.contains(hit)) };
  })()`);
  check(!!at && at.mine, `${what} is what stands under a press at its middle`);
  if (!at) return;
  if (phone) {
    await send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x: at.x, y: at.y }] });
    await send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  } else {
    await send('Input.dispatchMouseEvent', { type: 'mousePressed', x: at.x, y: at.y, button: 'left', clickCount: 1 });
    await send('Input.dispatchMouseEvent', { type: 'mouseReleased', x: at.x, y: at.y, button: 'left', clickCount: 1 });
  }
}
const panelIn = `(function () { var t = document.getElementById('mfs-content-panel').style.transform;
  return t === 'translateX(0px)' || t === 'translateX(0)'; })()`;
const studioShown = `(${panelIn} || getComputedStyle(document.documentElement).getPropertyValue('--mfs-phone').trim() === '1') && !!document.querySelector('#mfs-content-panel .mfs-studio')`;

// The studio's list as it stands now, read here apart from the page: released games with an
// App Store address carrying an identifier, and never this application.
const apps = (await (await fetch('https://owlsnestcreations.com/data/games.json')).json()).games
  .map(game => ({ game, id: (/\/id([0-9]+)$/.exec(game.appStoreUrl || '') || [])[1] }))
  .filter(({ game, id }) => game.released === true && id && id !== OWN_STORE_ID)
  .map(({ game, id }) => ({ store_id: id, name: game.name, line: game.tagline }));
const LINKS = [["See what's hatching in the nest...", 'https://owlsnestcreations.com'],
               ["Check out Damian's research page...", 'https://damiansowinski.com']];

await send('Page.navigate', { url: `${base}/MFS/` });
await until("document.readyState === 'complete' && document.querySelectorAll('#mfs-search-results .mfs-result').length > 0", 'the list is drawn');
// The banner's words flicker in and the list's panels slide in.
await until(`document.getElementById('mfs-studio-sign').textContent === 'Studio' && getComputedStyle(document.getElementById('mfs-search-panel')).pointerEvents === 'all'`,
            'Studio and the list have come in');
await sleep(1500);

// Where Studio stands.
const more = await box('mfs-studio-sign'), exit = await box('mfs-exit'), title = await box('mfs-title');
const face = await evaluate(`(function () {
  var m = getComputedStyle(document.getElementById('mfs-studio-sign')), e = getComputedStyle(document.getElementById('mfs-exit'));
  var letters = [].map.call(document.querySelectorAll('#mfs-studio-sign span:not(.mfs-word)'), function (s) { return getComputedStyle(s).color; });
  return { same: m.fontFamily === e.fontFamily && m.fontSize === e.fontSize && m.letterSpacing === e.letterSpacing && m.color === e.color,
           more: [m.fontFamily, m.fontSize, m.color].join(' '), exit: [e.fontFamily, e.fontSize, e.color].join(' '),
           settled: letters.every(function (c) { return c === m.color; }) };
})()`);
check(face.same, `Studio is set as the exit sign is (${face.more} against ${face.exit})`);
check(face.settled, 'every letter of Studio has come to rest in the sign\'s colour');
// The phone's column is the layout in force on a phone and on a desktop at a text size too large
// for panels side by side, which the page says in --mfs-phone.
const column = await evaluate("getComputedStyle(document.documentElement).getPropertyValue('--mfs-phone').trim() === '1'");
if (column) {
  check(more.top >= exit.bottom && Math.abs(more.right - exit.right) <= 1,
        `Studio stands under the sign, at its end (${more.top.toFixed(1)} under ${exit.bottom.toFixed(1)}, ending at ${more.right.toFixed(1)} against ${exit.right.toFixed(1)})`);
  const lines = await evaluate(`(function () { var t = document.getElementById('mfs-title');
    return t.getBoundingClientRect().height / parseFloat(getComputedStyle(t).lineHeight); })()`);
  if (!side && text === 16) check(Math.round(lines) === 2, `the title keeps its two lines beside the signs (${lines.toFixed(2)})`);
} else {
  check(Math.abs(more.top - exit.top) <= 1 && Math.abs(more.bottom - exit.bottom) <= 1,
        `Studio stands on the sign's line (${more.top.toFixed(1)} to ${more.bottom.toFixed(1)} against ${exit.top.toFixed(1)} to ${exit.bottom.toFixed(1)})`);
  check(more.right < exit.left, `Studio stands before the sign (${more.right.toFixed(1)} before ${exit.left.toFixed(1)})`);
  const edge = (await evaluate('document.documentElement.clientWidth')) - exit.right;
  check(Math.abs(exit.left - more.right - edge) <= 1,
        `Studio stands as far from the sign as the sign from the window's edge (${(exit.left - more.right).toFixed(1)} against ${edge.toFixed(1)})`);
}
check(more.left >= title.right || more.top >= title.bottom, 'Studio lies over no part of the title');

async function holdsTheStudioPage(when) {
  await until(studioShown, `${when}: the page of the studio's other applications shows`);
  // An icon is the studio site's picture, fetched as the page is drawn, so its arrival is
  // waited for rather than raced.
  await until(`[].every.call(document.querySelectorAll('#mfs-content-panel .mfs-studio img'), function (i) { return i.complete; })`,
              `${when}: every icon has arrived or failed`);
  await sleep(500);
  const page = await evaluate(`(function () {
    var panel = document.getElementById('mfs-content-panel');
    return {
      parts: [].map.call(panel.children, function (c) { return c.classList[0]; }),
      heading: panel.querySelector('.mfs-header .mfs-title').textContent,
      links: [].map.call(panel.querySelectorAll('.mfs-studio > .mfs-studio-link'), function (a) {
        return [a.textContent, a.getAttribute('href'), a.target];
      }),
      order: [].map.call(panel.querySelectorAll('.mfs-studio > *'), function (a) { return a.classList.contains('mfs-studio-link') ? 'line' : 'app'; }).join(' '),
      apps: [].map.call(panel.querySelectorAll('.mfs-studio > :not(.mfs-studio-link)'), function (app) {
        var a = app.querySelector('a.mfs-studio-store'), img = a && a.querySelector('img'), words = app.querySelector('.mfs-studio-text');
        var w = words.getBoundingClientRect(), i = (img || a).getBoundingClientRect();
        return { tag: a ? a.tagName : null, href: a && a.getAttribute('href'), target: a && a.target, rel: a && a.rel,
                 name: app.querySelector('.mfs-studio-name').textContent, line: app.querySelector('.mfs-studio-line').textContent,
                 drawn: !!img && img.complete && img.naturalWidth > 0, side: img ? img.getBoundingClientRect().width : 0,
                 under: i.top >= w.bottom, middle: Math.abs((i.left + i.right) / 2 - (w.left + w.right) / 2) <= 1 };
      }),
      // Every word of the page against the coffee panel's paragraph, the page's about panel:
      // its face, weight, colour, spacing and leading, at the size of the panel's written areas.
      faces: (function () {
        var about = getComputedStyle(document.getElementById('mfs-coffee-text'));
        var prose = parseFloat(getComputedStyle(panel.querySelector('.mfs-studio')).getPropertyValue('--mfs-prose')) *
                    parseFloat(getComputedStyle(document.documentElement).fontSize);
        return [].map.call(panel.querySelectorAll('.mfs-studio-text'), function (t) {
          var c = getComputedStyle(t), lead = function (s) { return parseFloat(s.lineHeight) / parseFloat(s.fontSize); };
          var same = c.fontFamily === about.fontFamily && c.fontWeight === about.fontWeight && c.fontStyle === about.fontStyle &&
                     c.color === about.color && Math.abs(lead(c) - lead(about)) < 0.01 &&
                     Math.abs(parseFloat(c.letterSpacing) / parseFloat(c.fontSize) - parseFloat(about.letterSpacing) / parseFloat(about.fontSize)) < 0.001 &&
                     Math.abs(parseFloat(c.fontSize) - prose) < 0.01 && c.textDecorationLine === 'none' &&
                     c.textShadow === 'none' && c.backgroundImage === 'none';
          return same ? null : t.textContent + ': ' + [c.fontFamily, c.fontWeight, c.color, c.lineHeight, c.letterSpacing, c.fontSize].join(' ') +
            ' against ' + [about.fontFamily, about.fontWeight, about.color, about.lineHeight, about.letterSpacing, prose + 'px'].join(' ');
        }).filter(Boolean);
      })(),
      active: document.querySelectorAll('#mfs-search-results .mfs-result-active').length,
      address: location.search,
      listShows: !document.getElementById('mfs-list-pane').classList.contains('mfs-pane-off')
    };
  })()`);
  check(JSON.stringify(page.parts) === JSON.stringify(['mfs-header', 'mfs-studio']), `${when}: the page holds its heading and its links and nothing else (${page.parts.join(', ')})`);
  check(page.heading === "More from Owl's Nest Creations", `${when}: the heading reads More from Owl's Nest Creations`);
  check(page.apps.length === apps.length, `${when}: one link for each of the studio's ${apps.length} other applications`);
  page.apps.forEach((shown, i) => {
    const app = apps[i] || {};
    check(shown.tag === 'A' && shown.href === `https://apps.apple.com/app/id${app.store_id}`,
          `${when}: ${shown.name} goes to ${shown.href}`);
    check(shown.target === '_blank' && /noopener/.test(shown.rel), `${when}: ${shown.name} opens in a tab of its own`);
    check(shown.name === app.name && shown.line === app.line, `${when}: ${shown.name} says its name and its line, "${shown.line}"`);
    check(shown.drawn && Math.abs(shown.side - 60) < 0.5, `${when}: ${shown.name}'s icon is drawn, ${shown.side} px across`);
    check(shown.under && shown.middle, `${when}: ${shown.name}'s icon stands in the middle under its words, as the coffee's button does`);
  });
  check(page.faces.length === 0, `${when}: every word of the page is set as the coffee panel's paragraph is${page.faces.length ? ' (' + page.faces.join('; ') + ')' : ''}`);
  check(JSON.stringify(page.links) === JSON.stringify(LINKS.map(([w, a]) => [w, a, '_blank'])),
        `${when}: the captain's two lines stand in his words, each opening its site in a tab of its own (${JSON.stringify(page.links)})`);
  check(page.order === apps.map(() => 'app').concat(['line', 'line']).join(' '), `${when}: the two lines stand under the applications (${page.order})`);
  check(page.active === 0, `${when}: no spacetime is marked open in the list`);
  check(page.address === '', `${when}: the address names no spacetime (${page.address})`);
  check(page.listShows, `${when}: the list shows beside it`);
}

await press('#mfs-studio-sign', 'Studio');
await holdsTheStudioPage('from the list');

// From an open spacetime: the page takes the spacetime's place, and the spacetime comes back.
const first = await evaluate("document.querySelector('#mfs-search-results .mfs-result').dataset.id");
const firstName = (await (await fetch(`${base}/MFS/assets/data/metrics/${first}.json`)).json()).name;
await press(`#mfs-search-results .mfs-result[data-id="${first}"]`, `${first}'s name in the list`);
await until(`(${panelIn} || getComputedStyle(document.documentElement).getPropertyValue('--mfs-phone').trim() === '1') && !document.querySelector('#mfs-content-panel .mfs-studio') && !!document.querySelector('#mfs-content-panel .mfs-header .mfs-title')`, `${first} opens`, 60);
await sleep(800);
await press('#mfs-studio-sign', 'Studio over an open spacetime');
await holdsTheStudioPage(`over ${first}`);
await press(`#mfs-search-results .mfs-result[data-id="${first}"]`, `${first}'s name in the list again`);
await until(`!document.querySelector('#mfs-content-panel .mfs-studio') && !!document.querySelector('#mfs-content-panel .mfs-header .mfs-title')`, `${first} opens again`, 60);
const drawnAgain = await evaluate("document.querySelector('#mfs-content-panel .mfs-header .mfs-title').textContent");
check(drawnAgain === firstName, `${first} is drawn again after the studio's page, as ${drawnAgain}`);

const thrown = await evaluate('window.__mfsErrors');
errors.push(...thrown.map(e => `page error: ${e}`));
errors.forEach(e => check(false, e));
console.log(failures.length ? `${failures.length} failed` : 'every check passed');
finish(failures.length ? 1 : 0);

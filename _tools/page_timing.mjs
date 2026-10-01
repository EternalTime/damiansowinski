#!/usr/bin/env node
/* Open every spacetime on the spacetimes page in headless Chrome, as a reader does, and time
   how long each takes to be ready: from the click on its name in the list, or the choice of a
   chart in the selector, until its mathematics is typeset and the page answers again.

   It measures in the page itself, so the time includes whatever keeps the tab busy after
   MathJax is done, which is where the Natário warp drive's minute went until 25 September
   2026, and it reports the longest single task that blocked the page, which is the thing a
   reader feels as a frozen tab. It needs Chrome and a served copy of the site, and nothing
   else; Node 22 or later carries the WebSocket it talks to Chrome over:

     bundle exec jekyll serve
     node _tools/page_timing.mjs http://127.0.0.1:4000

   --metric <id> times one spacetime and is repeatable; --phone lays the page out as an
   iPhone held upright, 390 by 844; --width and --height set another screen; --text <px>
   sets the browser's default font size, 16 by default, as a reader who enlarges text does;
   --cpu <n> slows the processor n times over, as Chrome's own tools do to stand in for a
   slower phone; --still asks for reduced motion, as a reader's system setting does, under
   which the list and the contents must swap without a slide; and --budget <seconds>, 10 by default, is the time past which a chart
   counts as failing. _tools/chrome.mjs starts Chrome, and CHROME names the browser when it
   is not where macOS keeps it.

   It also watches every number and label of every drawing shown, and the plot of each
   spacetime diagram, from the first frame the reader sees them until two frames after the
   last font has loaded, as each spacetime opens, as each chart is chosen and as each view of
   a diagram is chosen, and counts any that moves or is shown or hidden after that first frame
   as an error: until 29 September 2026 the first spacetime to use one of MathJax's fonts,
   Bianchi's fraktur or Malament-Hogarth's small capitals, fitted its numbers in the face
   drawn in its place and moved them by up to 4 pixels when the font arrived.

   In each of those states it holds every number, name of an axis and label of every drawing
   shown, a spacetime diagram, a conformal diagram, a figure in three dimensions and an
   embedding diagram, to the margin of 1.5 em of the drawing's words that the page gives them,
   and counts any that stands nearer the edge of the drawing's ground as an error. It holds
   them to the captain's order of 30 September 2026 as well: every word at least the size of
   the caption beside its drawing and in the page's white, no two labels over each other, and
   a spacetime diagram's tick marks at least half an em of its numbers long and a tenth of one
   thick. --print does the same for the print copy of each chart, as the page lays it out for
   paper once its print button is pressed, the browser's print dialog left closed.

   Nothing on the page glows, as the captain asked on 29 September 2026: "No glow anywhere".
   In each of those states it counts as an error every element whose computed style carries
   a text or box shadow or a filter that blurs or casts a shadow, every rule of the page's
   stylesheets that gives one, whatever its selector, every shadow drawn on a canvas and every
   one a script sets, even for a frame, such as the title's letters as they flicker in. The
   chosen button of each row and the spacetime shown in the list are pink, and every kind of
   button is pink while it is pressed and hovered over, which Chrome's own tools force on it.

   The spacetime's name, the word above it and every button keep one size in px whatever the
   reader's text size, as the captain asked on 29 September 2026, and in each state any of
   them set larger than that size counts as an error.

   Wherever the list and the spacetime stand side by side, in each state the list's panel
   starts on the spacetime panel's top and the coffee panel ends on its bottom, with 6px
   between the two on the left, and the space between the Exit sign and the spacetime panel
   is half what it was before 29 September 2026, as the captain asked that day; any of them a
   pixel out counts as an error.

   The Conventions section follows the chart, as the captain asked on 29 September 2026: in
   each state the paragraph under the signature reads exactly the chosen chart's convention
   followed by the one its spacetime shares across charts, as the metric file gives them,
   each formula read back from MathJax as the TeX it was set from, and choosing another chart
   changes it. Any other text there counts as an error.

   The list's panel shows the contents of the spacetime that is open, as the captain asked on
   1 October 2026. Each spacetime is opened as a reader does it, from the list searched and
   scrolled, and in each state the contents must read exactly the headings of the page as
   drawn, in its order, with the list slid out of the panel and out of the keyboard's reach.
   An entry must bring its section to the top of what shows of the spacetime, under the name
   where the name stays in sight, or as far up as the page scrolls. "‹ All spacetimes" must
   bring the list back with the search as it was typed and the list where it was scrolled to,
   and the open spacetime's name in the list must bring its contents back without drawing it
   again. Last, the page is opened at the address of a spacetime, ?spacetime= and its id,
   and must come up on that spacetime with its contents showing and nothing having slid.

   It prints one line per chart, slowest last, then every page error and console error the
   page raised and every label that moved, and exits non-zero if any chart missed the budget
   or never became ready, or if the page raised any error at all. */
import { launch, sleep } from './chrome.mjs';

const argv = process.argv.slice(2);
function option(name, fallback) {
  const i = argv.indexOf(name);
  return i >= 0 ? argv[i + 1] : fallback;
}
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const only = argv.flatMap((a, i) => (a === '--metric' ? [argv[i + 1]] : []));
const phone = argv.includes('--phone');
const width = +option('--width', phone ? 390 : 1440);
const height = +option('--height', phone ? 844 : 900);
const text = +option('--text', 16);
const cpu = +option('--cpu', 1);
const budget = +option('--budget', 10);
const print = argv.includes('--print');
const still = argv.includes('--still');
const LIMIT_S = Math.max(300, budget * 10);

const { send, evaluate, finish, onError } = await launch({ width, height, phone, text, cpu });
/* Every console.error and every error the browser logs itself, such as a file that failed to
   load, with the spacetime open when it came. The page's own errors are gathered in the page,
   below, since Chrome reports "ResizeObserver loop completed with undelivered notifications"
   to the page's error handlers alone. */
if (still) await send('Emulation.setEmulatedMedia', { features: [{ name: 'prefers-reduced-motion', value: 'reduce' }] });
const errors = [];
let opened = 'the list';
onError(message => errors.push(`${opened}: ${message}`));

/* Before the page's own scripts run: count the typesetting MathJax has in hand, keep every
   task that blocked the page for more than 50ms, which is what the browser calls a long task,
   and keep every error the page raises. */
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsTiming = { busy: 0, long: [], errors: [], slides: 0 };
  document.addEventListener('transitionrun', function (e) {
    if (e.target.classList && e.target.classList.contains('mfs-pane') && e.propertyName === 'transform') window.__mfsTiming.slides++;
  }, true);
  window.addEventListener('error', function (e) {
    window.__mfsTiming.errors.push(e.message + (e.filename ? ' (' + e.filename + ':' + e.lineno + ')' : ''));
  });
  window.addEventListener('unhandledrejection', function (e) {
    window.__mfsTiming.errors.push('unhandled rejection: ' + (e.reason && e.reason.stack || e.reason));
  });
  try {
    new PerformanceObserver(function (list) {
      list.getEntries().forEach(function (e) { window.__mfsTiming.long.push([e.startTime, e.duration]); });
    }).observe({ entryTypes: ['longtask'] });
  } catch (e) {}
  (function wrap() {
    if (!(window.MathJax && MathJax.typesetPromise && MathJax.startup && MathJax.startup.document)) { setTimeout(wrap, 20); return; }
    var typeset = MathJax.typesetPromise;
    MathJax.typesetPromise = function () {
      window.__mfsTiming.busy++;
      function done() { window.__mfsTiming.busy--; }
      return typeset.apply(this, arguments).then(function (v) { done(); return v; }, function (e) { done(); throw e; });
    };
    window.__mfsTiming.ready = true;
  })();
  /* Where every number and label of every drawing shown stands, and each spacetime diagram's
     plot, from the top left of its figure, since the panel slides in, and whether it is shown. */
  window.__mfsLabels = function (panel) {
    var out = {};
    [].forEach.call(panel.querySelectorAll('.mfs-nr-figure, .mfs-cd-figure'), function (fig, i) {
      if (!fig.getClientRects().length) return;
      var b = fig.getBoundingClientRect();
      [].forEach.call(fig.querySelectorAll('.nr-plot, .nr-tick, .mfs-slice-label, .cd-at'), function (el, j) {
        if (!el.getClientRects().length) return;
        var r = el.getBoundingClientRect();
        out[i + ' ' + j] = { name: el.className.split(' ')[0] + ' "' + el.textContent.trim().slice(0, 24) + '"',
          at: [r.left - b.left, r.top - b.top, r.right - b.left, r.bottom - b.top],
          shown: getComputedStyle(el).visibility !== 'hidden' };
      });
    });
    return out;
  };
  /* From the first frame at which shown() holds, every frame is held against that first one
     until stop() is called, and each label that moved by more than a twentieth of a pixel, the
     noise of the panel's slide, or was shown or hidden, is kept once with how far it went. */
  window.__mfsWatch = function (panel, shown) {
    var first = null, start = 0, moved = {}, on = true;
    function compare() {
      var now = window.__mfsLabels(panel);
      Object.keys(now).forEach(function (k) {
        var a = first[k], b = now[k];
        if (!a || moved[k]) return;
        var d = Math.max.apply(null, a.at.map(function (v, i) { return Math.abs(v - b.at[i]); }));
        if (d > 0.05 || a.shown !== b.shown) {
          moved[k] = a.name + (a.shown !== b.shown ? (b.shown ? ' shown' : ' hidden') : ' moved ' + d.toFixed(2) + 'px') +
            ' at ' + Math.round(performance.now() - start) + 'ms after first shown';
        }
      });
    }
    (function frame() {
      if (!on) return;
      if (first) compare();
      else if (shown()) { first = window.__mfsLabels(panel); start = performance.now(); }
      requestAnimationFrame(frame);
    })();
    return function stop() {
      on = false;
      if (first) compare();
      return Object.keys(moved).map(function (k) { return moved[k]; });
    };
  };
  /* The words of a drawing's ground, a spacetime diagram's .mfs-nr-figure or the
     .mfs-cd-drawing of a conformal diagram, a figure in three dimensions or an embedding
     diagram, that a reader sees: its numbers, the names of its axes and its labels. */
  window.__mfsWords = function (ground) {
    var words = ground.matches('.mfs-nr-figure') ? '.nr-tick, .nr-name, .nr-ref, .mfs-slice-label' : '.cd-at';
    return [].filter.call(ground.querySelectorAll(words), function (el) {
      return el.getClientRects().length && getComputedStyle(el).visibility !== 'hidden';
    });
  };
  /* What a drawing is, for a message: its kind and the name its SVG gives it. */
  window.__mfsDrawingName = function (ground) {
    var figure = ground.closest('.mfs-cd-figure'), svg = ground.querySelector('svg');
    var kind = !figure ? 'spacetime diagram' : figure.classList.contains('mfs-em-figure') ? 'embedding diagram'
      : figure.classList.contains('mfs-pj-figure') ? 'figure' : 'conformal diagram';
    var name = svg && svg.getAttribute('aria-label') || '';
    return kind + ' "' + name.replace(/^[^,]*, /, '').slice(0, 48) + '"';
  };
  /* Every drawing shown in root whose words are not all readable, as the captain asked on
     30 September 2026: every word at least the size of the caption beside it and in the
     page's white, --text-bright, no two over each other, and a spacetime diagram's tick marks
     at least half an em of its numbers long and a tenth of one thick. One line for each
     drawing, with the numbers. On paper, where every word prints black, the colour is not held. */
  window.__mfsReadable = function (root, paper) {
    var probe = document.body.appendChild(document.createElement('i'));
    probe.style.color = 'var(--text-bright)';
    var white = getComputedStyle(probe).color, out = [];
    probe.remove();
    [].forEach.call(root.querySelectorAll('.mfs-nr-figure, .mfs-cd-drawing'), function (ground) {
      if (!ground.getClientRects().length) return;
      var view = ground.closest('.mfs-nr-view'), caption = view && view.querySelector('.mfs-nr-caption p');
      var least = caption ? parseFloat(getComputedStyle(caption).fontSize) : 0;
      var words = window.__mfsWords(ground), small = null, colours = {}, faults = [];
      words.forEach(function (el) {
        var size = parseFloat(getComputedStyle(el).fontSize);
        if (!small || size < small.size) small = { size: size, text: el.textContent.trim().slice(0, 16) };
        if (!paper) [el].concat([].slice.call(el.querySelectorAll('*'))).some(function (e) {
          var c = getComputedStyle(e).color;
          if (c !== white) { colours[c] = (colours[c] || 0) + 1; return true; }
        });
      });
      if (small && small.size < least - 0.01) {
        faults.push('words from ' + small.size.toFixed(1) + 'px ("' + small.text + '"), below the caption, ' + least.toFixed(1) + 'px');
      }
      var grey = Object.keys(colours);
      if (grey.length) {
        faults.push(grey.reduce(function (n, c) { return n + colours[c]; }, 0) + ' of ' + words.length +
          ' words not white ' + white + ' but ' + grey.join(', '));
      }
      if (ground.matches('.mfs-nr-figure')) {
        var tick = ground.querySelector('.nr-tick'), em = tick ? parseFloat(getComputedStyle(tick).fontSize) : 0;
        var marks = [].map.call(ground.querySelectorAll('.nr-tickmark'), function (m) {
          var r = m.getBoundingClientRect(), x = m.closest('.nr-xticks');
          return { long: x ? r.height : r.width, thick: x ? r.width : r.height };
        });
        if (!marks.length) {
          var svg = ground.querySelector('.mfs-nr-svg'), scale = svg.getBoundingClientRect().width / svg.viewBox.baseVal.width;
          marks = [].map.call(svg.querySelectorAll('line.nr-axis'), function (m) {
            var r = m.getBoundingClientRect();
            return { long: Math.max(r.width, r.height), thick: parseFloat(getComputedStyle(m).strokeWidth) * scale };
          });
        }
        var longest = marks.reduce(function (a, m) { return Math.min(a, m.long); }, Infinity);
        var thickest = marks.reduce(function (a, m) { return Math.min(a, m.thick); }, Infinity);
        // Less a 64th of a pixel, the unit a browser lays out in: a tenth of 42px is drawn 4.1875px.
        if (marks.length && (longest < 0.5 * em - 0.016 || thickest < 0.1 * em - 0.016)) {
          faults.push('tick marks ' + longest.toFixed(1) + 'px long and ' + thickest.toFixed(2) + 'px thick, below ' +
            (0.5 * em).toFixed(1) + 'px and ' + (0.1 * em).toFixed(2) + 'px');
        }
      }
      // Two labels over each other by more than half a pixel each way; they may touch.
      var boxes = words.map(function (el) { return el.getBoundingClientRect(); });
      boxes.forEach(function (a, i) {
        boxes.slice(0, i).forEach(function (b, j) {
          if (a.left < b.right - 0.5 && b.left < a.right - 0.5 && a.top < b.bottom - 0.5 && b.top < a.bottom - 0.5) {
            faults.push('"' + words[j].textContent.trim().slice(0, 16) + '" and "' + words[i].textContent.trim().slice(0, 16) + '" overlap');
          }
        });
      });
      if (faults.length) out.push(window.__mfsDrawingName(ground) + ': ' + faults.join('; '));
    });
    return out;
  };
  /* Every number, name of an axis and label of each drawing shown in root that stands nearer
     the edge of the drawing's ground than 1.5 em of its words, the margin .mfs-nr-figure and
     .mfs-cd-drawing give them, with how near it stands, in those ems. */
  window.__mfsMargins = function (root) {
    var out = [];
    [].forEach.call(root.querySelectorAll('.mfs-nr-figure, .mfs-cd-drawing'), function (fig) {
      if (!fig.getClientRects().length) return;
      var b = fig.getBoundingClientRect(), em = parseFloat(getComputedStyle(fig).fontSize);
      window.__mfsWords(fig).forEach(function (el) {
        var r = el.getBoundingClientRect();
        var near = Math.min(r.left - b.left, r.top - b.top, b.right - r.right, b.bottom - r.bottom);
        // Less half a pixel, to which a browser may round where a turned or shifted box stands.
        if (near < 1.5 * em - 0.5) {
          out.push(window.__mfsDrawingName(fig) + ': ' + el.className.split(' ')[0] + ' "' + el.textContent.trim().slice(0, 24) + '" stands ' +
            (near / em).toFixed(2) + ' em from the edge of its diagram, inside the margin of 1.5 em');
        }
      });
    });
    return out;
  };
  /* Every glow on the page: an element whose computed style casts a shadow or blurs, a rule
     of a stylesheet that gives one, and a canvas that draws one. What a script sets is caught
     at the next frame, and a canvas as it is set, and both are kept in __mfsGlowSeen. */
  (function () {
    var seen = window.__mfsGlowSeen = [];
    function name(el) {
      return el.tagName.toLowerCase() + (el.id ? '#' + el.id : '') +
        (typeof el.className === 'string' && el.className ? '.' + el.className.trim().split(/\\s+/).join('.') : '') +
        (el.textContent ? ' "' + el.textContent.trim().slice(0, 24) + '"' : '');
    }
    function glow(style) {
      var out = [];
      if (style.textShadow && style.textShadow !== 'none') out.push('text-shadow ' + style.textShadow);
      if (style.boxShadow && style.boxShadow !== 'none') out.push('box-shadow ' + style.boxShadow);
      if (/blur|drop-shadow/.test(style.filter || '')) out.push('filter ' + style.filter);
      return out;
    }
    window.__mfsGlowOf = function (el) {
      var out = glow(getComputedStyle(el));
      if (/^(filter|fegaussianblur|fedropshadow)$/i.test(el.tagName)) out.push('an SVG ' + el.tagName);
      return out.map(function (g) { return name(el) + ' has ' + g; });
    };
    window.__mfsGlow = function (root) {
      var out = [];
      [root].concat([].slice.call(root.getElementsByTagName('*'))).forEach(function (el) {
        out = out.concat(window.__mfsGlowOf(el));
      });
      [].forEach.call(document.styleSheets, function (sheet) {
        var where = sheet.href || 'the page', list;
        try { list = sheet.cssRules; } catch (e) { return; }
        (function rules(list) {
          [].forEach.call(list, function (rule) {
            if (rule.cssRules) rules(rule.cssRules);
            if (rule.style) glow(rule.style).forEach(function (g) {
              out.push('the rule ' + (rule.selectorText || rule.cssText.slice(0, 40)) + ' in ' + where + ' gives ' + g);
            });
          });
        })(list);
      });
      return out.concat(seen.splice(0));
    };
    var ctx = CanvasRenderingContext2D.prototype;
    ['shadowBlur', 'shadowOffsetX', 'shadowOffsetY', 'globalCompositeOperation'].forEach(function (k) {
      var d = Object.getOwnPropertyDescriptor(ctx, k);
      Object.defineProperty(ctx, k, { configurable: true, enumerable: d.enumerable, get: d.get, set: function (v) {
        if (k === 'globalCompositeOperation' ? /lighter|screen/.test(v) : +v) {
          var what = 'canvas' + (this.canvas.id ? '#' + this.canvas.id : '') + ' draws with ' + k + ' ' + v;
          if (seen.indexOf(what) < 0) seen.push(what);
        }
        d.set.call(this, v);
      } });
    });
    var changed = new Set(), frame = 0;
    new MutationObserver(function (records) {
      records.forEach(function (r) { changed.add(r.target); });
      frame = frame || requestAnimationFrame(function () {
        frame = 0;
        changed.forEach(function (el) {
          if (el.isConnected) window.__mfsGlowOf(el).forEach(function (g) { if (seen.indexOf(g) < 0) seen.push(g); });
        });
        changed.clear();
      });
    }).observe(document, { attributes: true, attributeFilter: ['style', 'class'], subtree: true });
    /* The chosen button of every row, and the spacetime shown in the list, that is not pink. */
    window.__mfsPink = function () {
      var probe = document.body.appendChild(document.createElement('i'));
      probe.style.color = 'var(--pink-light)';
      var pink = getComputedStyle(probe).color, out = [];
      probe.remove();
      [].forEach.call(document.querySelectorAll('#mfs-content-panel .mfs-choice-on, .mfs-result-active'), function (el) {
        // A button turns pink over its colour's transition, which is measured at its end.
        el.getAnimations().forEach(function (a) { a.finish(); });
        var cs = getComputedStyle(el), border = el.classList.contains('mfs-choice');
        if (cs.color !== pink || (border && cs.borderTopColor !== pink)) {
          out.push(name(el) + ' is chosen and ' + cs.color + (border ? ' bordered ' + cs.borderTopColor : '') + ', not pink ' + pink);
        }
      });
      return out;
    };
  })();
  /* The name, the word above it or a button set larger than its one size in px, the
     desktop's and then the phone's; fitWords() may set one smaller, never larger. */
  window.__mfsSizes = function () {
    var phone = getComputedStyle(document.documentElement).getPropertyValue('--mfs-phone').trim() === '1', out = [];
    [['#mfs-content-panel .mfs-header .mfs-title', 40, 26], ['#mfs-content-panel .mfs-header .mfs-label', 21, 14],
     ['#mfs-print-btn', 15, 12], ['#mfs-content-panel .mfs-choice', 15, 12], ['.mfs-result', 18, 18],
     ['.mfs-toc-entry', 18, 18], ['.mfs-toc-back', 18, 18]].forEach(function (s) {
      [].forEach.call(document.querySelectorAll(s[0]), function (el) {
        var size = parseFloat(getComputedStyle(el).fontSize), fixed = phone ? s[2] : s[1];
        if (el.getClientRects().length && size > fixed + 0.01) {
          out.push(s[0] + ' "' + el.textContent.trim().slice(0, 24) + '" is set at ' + size + 'px, above its ' + fixed + 'px');
        }
      });
    });
    return out;
  };
  /* Wherever the panels stand side by side: the top of the list's panel on the top of the
     spacetime's panel and the bottom of the coffee panel on its bottom, within a pixel, the
     list's panel 6px above the coffee panel, and the spacetime panel's top halfway between
     the foot of the Exit sign and where it stood before 29 September 2026, 150px down or 60px
     below the title's foot, or 30px below the title's foot, whichever is lower. */
  window.__mfsPanelEdges = function () {
    if (getComputedStyle(document.documentElement).getPropertyValue('--mfs-phone').trim() === '1') return [];
    function box(id) { return document.getElementById(id).getBoundingClientRect(); }
    var list = box('mfs-search-panel'), coffee = box('mfs-coffee-panel'), right = box('mfs-content-panel');
    var sign = box('mfs-exit').bottom, foot = box('mfs-title').bottom, out = [];
    function near(what, a, b) { if (Math.abs(a - b) > 1) out.push(what + ' is at ' + a.toFixed(1) + ', not ' + b.toFixed(1)); }
    near("the top of the list's panel", list.top, right.top);
    near('the bottom of the coffee panel', coffee.bottom, right.bottom);
    near("the gap between the list's panel and the coffee panel", coffee.top - list.bottom, 6);
    var old = Math.max(150, Math.ceil(foot + 60));
    near("the top of the spacetime's panel", right.top, Math.max(sign + (old - sign) / 2, foot + 30));
    return out;
  };
  /* Which of the list panel's two views shows, once any slide has ended: the one asked for
     fills the panel inside its padding and answers the keyboard, and the other stands wholly
     outside it, hidden and inert. */
  window.__mfsShowing = function (contents) {
    var list = document.getElementById('mfs-list-pane'), toc = document.getElementById('mfs-toc-pane'), out = [];
    [list, toc].forEach(function (pane) { pane.getAnimations().forEach(function (a) { a.finish(); }); });
    var room = document.getElementById('mfs-panes').getBoundingClientRect();
    var shown = contents ? toc : list, off = contents ? list : toc;
    var what = contents ? 'the contents' : 'the list', other = contents ? 'the list' : 'the contents';
    var a = shown.getBoundingClientRect(), b = off.getBoundingClientRect();
    if (Math.abs(a.left - room.left) > 0.5 || Math.abs(a.right - room.right) > 0.5) out.push(what + ' do not fill the panel, from ' + a.left.toFixed(1) + ' to ' + a.right.toFixed(1));
    if (shown.inert || getComputedStyle(shown).visibility !== 'visible') out.push(what + ' are out of reach');
    if (b.right > room.left + 0.5 && b.left < room.right - 0.5) out.push(other + ' still stand in the panel, from ' + b.left.toFixed(1) + ' to ' + b.right.toFixed(1));
    if (!off.inert || getComputedStyle(off).visibility !== 'hidden') out.push(other + ' can still be reached');
    return out;
  };
  /* The contents against the headings of the page as drawn, and the contents showing. */
  window.__mfsContents = function () {
    function texts(selector) {
      return [].map.call(document.querySelectorAll(selector), function (el) { return el.textContent.trim(); });
    }
    var heads = texts('#mfs-content-panel .mfs-section-label'), entries = texts('#mfs-toc .mfs-toc-entry'), out = [];
    if (!heads.length || heads.join('|') !== entries.join('|')) {
      out.push('the contents read ' + JSON.stringify(entries) + ', not the headings ' + JSON.stringify(heads));
    }
    return out.concat(window.__mfsShowing(true));
  };
  window.__mfsSettled = function () {
    return document.fonts.ready.then(function () {
      return new Promise(function (r) { requestAnimationFrame(function () { requestAnimationFrame(r); }); });
    });
  };` });
await send('Page.navigate', { url: `${base}/MFS/` });
for (let i = 0; i < 600; i++) {
  if (await evaluate(`!!(window.__mfsTiming && window.__mfsTiming.ready && document.querySelector('.mfs-result'))`)) break;
  await sleep(100);
}
const listed = await evaluate(`[].map.call(document.querySelectorAll('.mfs-result'), function (e) { return [e.dataset.id, e.textContent]; })`);
const ids = (listed || []).map(l => l[0]), names = Object.fromEntries(listed || []);
if (!ids || !ids.length) { console.error(`No spacetimes listed at ${base}/MFS/.`); finish(2); }
for (const id of only) if (!ids.includes(id)) { console.error(`No spacetime ${id} in the list.`); finish(2); }

/* Act, then resolve once the panel has been drawn again, holds the chart asked for with every
   formula typeset and no typesetting outstanding, and two frames have passed, so the
   observers that run once MathJax is done have run too. Every wait is a timer inside the
   page, so while a task holds the page nothing here moves, and the time taken includes that
   task. The labels are watched from the first frame at which the panel stands in place with
   its mathematics set and its figures in sight, until the page has settled and the fonts are
   in; that wait is not timed. */
function timed(action, chart) {
  return evaluate(`new Promise(function (resolve) {
    var t = window.__mfsTiming, panel = document.getElementById('mfs-content-panel');
    var before = panel.querySelector('.mfs-header');
    var stop = window.__mfsWatch(panel, function () {
      var header = panel.querySelector('.mfs-header');
      return header && header !== before && t.busy === 0 && !panel.classList.contains('mfs-fonts-wait') &&
        /^translateX\\(0(px)?\\)$/.test(panel.style.transform);
    });
    var start = performance.now();
    ${action}
    (function poll() {
      var on = panel.querySelector('.mfs-charts .mfs-choice-on'), header = panel.querySelector('.mfs-header');
      var ready = header && header !== before && t.busy === 0 &&
        (!on || on.dataset.chart === '${chart}') && !/\\\\[(\\[]/.test(panel.textContent);
      if (!ready) {
        if (performance.now() - start > ${LIMIT_S * 1000}) { stop(); resolve(null); return; }
        setTimeout(poll, 25); return;
      }
      requestAnimationFrame(function () { requestAnimationFrame(function () { setTimeout(function () {
        var end = performance.now();
        var longest = t.long.filter(function (l) { return l[0] + l[1] > start; })
                            .reduce(function (m, l) { return Math.max(m, l[1]); }, 0);
        var result = { seconds: (end - start) / 1000, longest: longest / 1000,
                       chart: on ? on.textContent : '',
                       lines: panel.querySelectorAll('.mfs-line').length,
                       elements: panel.getElementsByTagName('*').length };
        window.__mfsSettled().then(function () {
          result.moved = stop();
          result.margins = window.__mfsMargins(panel).concat(window.__mfsReadable(panel), window.__mfsGlow(document.documentElement), window.__mfsPink(),
            window.__mfsSizes(), window.__mfsPanelEdges(), window.__mfsContents());
          resolve(result);
        });
      }, 0); }); });
    })();
  })`);
}

/* Choose each view of each diagram but the one first shown, and watch its labels from the
   frame after the choice until the page has settled. */
async function views() {
  const moved = await evaluate(`(async function () {
    var panel = document.getElementById('mfs-content-panel'), moved = [];
    var buttons = [].slice.call(panel.querySelectorAll('.mfs-diagram > .mfs-section-head .mfs-choice'));
    for (var i = 0; i < buttons.length; i++) {
      if (buttons[i].classList.contains('mfs-choice-on')) continue;
      var button = buttons[i], clicked = false;
      var stop = window.__mfsWatch(panel, function () { return clicked; });
      button.click();
      clicked = true;
      await window.__mfsSettled();
      await new Promise(function (r) { requestAnimationFrame(r); });
      stop().forEach(function (m) { moved.push(button.textContent.trim() + ': ' + m); });
      window.__mfsMargins(panel).concat(window.__mfsReadable(panel), window.__mfsGlow(panel), window.__mfsPink(), window.__mfsSizes())
        .forEach(function (m) { moved.push(button.textContent.trim() + ': ' + m); });
    }
    return moved;
  })()`);
  for (const m of moved) errors.push(`${opened}: view ${m}`);
}

/* Press the print button with the browser's print dialog held closed, wait until the print
   copy is set, and hold its spacetime diagrams to their margin as paper lays them out. */
async function printCopy() {
  const ready = await evaluate(`new Promise(function (resolve) {
    var start = performance.now(), printed = false;
    window.print = function () { printed = true; resolve(true); };
    document.getElementById('mfs-print-btn').click();
    (function wait() {
      if (printed) return;
      if (performance.now() - start > ${LIMIT_S * 1000}) { resolve(false); return; }
      setTimeout(wait, 50);
    })();
  })`);
  if (!ready) { errors.push(`${opened}: the print copy was never set`); return; }
  await send('Emulation.setEmulatedMedia', { media: 'print' });
  for (const m of await evaluate(`window.__mfsMargins(document.getElementById('mfs-print-root'))
      .concat(window.__mfsReadable(document.getElementById('mfs-print-root'), true),
              window.__mfsGlow(document.getElementById('mfs-print-root')))`)) {
    errors.push(`${opened}: print: ${m}`);
  }
  await send('Emulation.setEmulatedMedia', { media: '' });
}

/* Hold each kind of button hovered over and pressed, as Chrome's own tools do, the first time
   a spacetime shows one, and count one that is not pink then, or that glows, as an error. */
const BUTTONS = new Set(['#mfs-content-panel .mfs-choice:not(.mfs-choice-on)', '#mfs-content-panel .mfs-choice-on',
  '#mfs-print-btn', '#mfs-content-panel .mfs-turn-reset', '.mfs-result:not(.mfs-result-active)',
  '.mfs-toc-entry', '.mfs-toc-back']);
async function pressed() {
  if (!BUTTONS.size) return;
  await send('DOM.enable');
  await send('CSS.enable');
  const root = (await send('DOM.getDocument', { depth: 0 })).result.root.nodeId;
  for (const selector of BUTTONS) {
    const node = (await send('DOM.querySelector', { nodeId: root, selector })).result?.nodeId;
    if (!node) continue;
    BUTTONS.delete(selector);
    await send('CSS.forcePseudoState', { nodeId: node, forcedPseudoClasses: ['hover', 'active'] });
    const found = await evaluate(`(function () {
      var el = document.querySelector(${JSON.stringify(selector)});
      el.getAnimations().forEach(function (a) { a.finish(); });
      var probe = document.body.appendChild(document.createElement('i'));
      probe.style.color = 'var(--pink-light)';
      var pink = getComputedStyle(probe).color, cs = getComputedStyle(el);
      probe.remove();
      var out = window.__mfsGlowOf(el);
      var border = el.matches('.mfs-result, .mfs-toc-entry, .mfs-toc-back') ? pink : cs.borderTopColor;
      if (cs.color !== pink || border !== pink) {
        out.push(${JSON.stringify(selector)} + ' pressed is ' + cs.color + ' bordered ' + cs.borderTopColor + ', not pink ' + pink);
      }
      return out;
    })()`);
    await send('CSS.forcePseudoState', { nodeId: node, forcedPseudoClasses: [] });
    for (const m of found) errors.push(`${opened}: ${m}`);
  }
  await send('CSS.disable');
  await send('DOM.disable');
}

/* Read the Conventions paragraph as it stands, each formula as the TeX MathJax set it from and
   the soft hyphens that divide a long word taken out, and hold it to the chart's convention
   followed by the spacetime's, as the metric file gives them. Returns what it read, so the
   next chart can be held to reading differently. */
async function convention(id, chart, before) {
  const found = await evaluate(`(async function () {
    var data = await (await fetch(${JSON.stringify(`${base}/MFS/assets/data/metrics/${id}.json`)})).json();
    var coord = (data.coordinates || [])[${chart}] || {};
    var expected = [coord.convention, data.convention].filter(function (t) { return t && t.trim(); })
      .map(function (t) { return t.trim(); }).join(' ');
    var box = document.querySelector('#mfs-content-panel .mfs-convention');
    var paragraphs = box ? [].filter.call(box.querySelectorAll('p'), function (p) {
      return !p.querySelector('.mfs-convention-k');
    }) : [];
    var roots = new Map();
    paragraphs.forEach(function (p) {
      MathJax.startup.document.getMathItemsWithin(p).forEach(function (m) { roots.set(m.typesetRoot, m.math); });
    });
    function read(node) {
      if (roots.has(node)) return '$' + roots.get(node) + '$';
      if (node.nodeType === 3) return node.data;
      return [].map.call(node.childNodes, read).join('');
    }
    var shown = paragraphs.map(read).join(' ').replace(/\u00AD/g, '').replace(/\\s+/g, ' ').trim();
    return { expected: expected, shown: shown, paragraphs: paragraphs.length };
  })()`);
  const at = `${opened} / chart ${chart}: conventions`;
  if (found.paragraphs !== (found.expected ? 1 : 0)) errors.push(`${at} stand in ${found.paragraphs} paragraphs`);
  if (found.shown !== found.expected) errors.push(`${at} read ${JSON.stringify(found.shown)}, not ${JSON.stringify(found.expected)}`);
  if (before !== undefined && found.shown === before) errors.push(`${at} did not change with the chart`);
  return found.shown;
}

/* Open a spacetime as a reader does: back to the list where the contents show, a letter of
   its name typed into the search, the list scrolled to its end, and the name pressed. What the
   search holds and where the list stands once pressed is kept for contents() to hold the list
   to when it comes back. */
function opening(id) {
  return `(function () {
    var input = document.getElementById('mfs-search-input'), results = document.getElementById('mfs-search-results');
    if (document.getElementById('mfs-list-pane').inert) document.getElementById('mfs-toc-back').click();
    input.value = ${JSON.stringify(names[id].trim().slice(0, 1).toLowerCase())};
    input.dispatchEvent(new Event('input'));
    results.scrollTop = results.scrollHeight;
    document.querySelector('.mfs-result[data-id="${id}"]').click();
    window.__mfsList = { value: input.value, scroll: results.scrollTop, slides: window.__mfsTiming.slides };
  })();`;
}

/* Press the last, the middle and the first entry of the contents and hold each section to its
   place, then the link back and hold the list to how it was left, then the open spacetime's
   name, which must bring the contents back over the same page. The search is left empty. */
async function contents() {
  const found = await evaluate(`(async function () {
    var out = [], kept = window.__mfsList, panel = document.getElementById('mfs-content-panel');
    var input = document.getElementById('mfs-search-input'), results = document.getElementById('mfs-search-results');
    var phone = getComputedStyle(document.documentElement).getPropertyValue('--mfs-phone').trim() === '1';
    var still = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    if (!still && window.__mfsTiming.slides === kept.slides) out.push('the contents came without a slide');
    if (still && window.__mfsTiming.slides) out.push('the views slid, ' + window.__mfsTiming.slides + ' times, under reduced motion');
    function at() { return phone ? window.scrollY : panel.scrollTop; }
    // Until what scrolls has stood still for a fifth of a second.
    function settled() {
      return new Promise(function (resolve) {
        var last = at(), since = performance.now(), start = since;
        (function poll() {
          var now = performance.now();
          if (at() !== last) { last = at(); since = now; }
          if (now - since > 200 || now - start > 5000) { resolve(); return; }
          requestAnimationFrame(poll);
        })();
      });
    }
    var entries = [].slice.call(document.querySelectorAll('#mfs-toc .mfs-toc-entry'));
    var labels = [].slice.call(panel.querySelectorAll('.mfs-section-label'));
    var order = [entries.length - 1, Math.floor(entries.length / 2), 0];
    for (var k = 0; k < order.length; k++) {
      var i = order[k];
      if (!entries[i] || !labels[i]) continue;
      entries[i].click();
      await settled();
      var head = panel.querySelector('.mfs-header'), loose = head.classList.contains('mfs-header-loose');
      var edge = phone ? 0 : panel.getBoundingClientRect().top + panel.clientTop + parseFloat(getComputedStyle(panel).paddingTop);
      var top = loose ? edge : head.getBoundingClientRect().bottom;
      var section = (labels[i].closest('.mfs-section') || labels[i]).getBoundingClientRect().top;
      var end = phone ? window.scrollY >= document.documentElement.scrollHeight - window.innerHeight - 1
                      : panel.scrollTop >= panel.scrollHeight - panel.clientHeight - 1;
      if (Math.abs(section - top) > 1 && !(end && section > top)) {
        out.push('the entry "' + entries[i].textContent + '" leaves its section at ' + section.toFixed(1) + ', not ' + top.toFixed(1));
      }
    }
    document.getElementById('mfs-toc-back').click();
    out = out.concat(window.__mfsShowing(false));
    if (input.value !== kept.value) out.push('the list came back with the search ' + JSON.stringify(input.value) + ', not ' + JSON.stringify(kept.value));
    if (Math.abs(results.scrollTop - kept.scroll) > 0.5) out.push('the list came back scrolled to ' + results.scrollTop + ', not ' + kept.scroll);
    if (!document.querySelector('.mfs-result-active')) out.push('the list came back without the open spacetime marked');
    input.value = '';
    input.dispatchEvent(new Event('input'));
    var page = panel.querySelector('.mfs-header');
    document.querySelector('.mfs-result-active').click();
    if (panel.querySelector('.mfs-header') !== page) out.push('the open spacetime was drawn again from the list');
    out = out.concat(window.__mfsContents());
    if (phone) { window.scrollTo(0, 0); await settled(); }
    return out;
  })()`);
  for (const m of found) errors.push(`${opened}: ${m}`);
}

/* Open the page at a spacetime's own address: it must come up on that spacetime, its contents
   showing from the start, with no slide between the two views. */
async function direct(id) {
  opened = `${id} by its address`;
  await send('Page.navigate', { url: `${base}/MFS/?spacetime=${id}` });
  let found = null;
  for (let i = 0; i < LIMIT_S * 10 && !found; i++) {
    await sleep(100);
    found = await evaluate(`(function () {
      var t = window.__mfsTiming, panel = document.getElementById('mfs-content-panel');
      if (!(t && t.ready && panel && panel.querySelector('.mfs-header') && t.busy === 0 &&
            /^translateX\\(0(px)?\\)$/.test(panel.style.transform))) return null;
      var out = window.__mfsContents();
      if (t.slides) out.push('the contents slid in, ' + t.slides + ' slides');
      var active = document.querySelector('.mfs-result-active');
      if (!active || active.dataset.id !== ${JSON.stringify(id)}) out.push('the list marks ' + (active ? active.dataset.id : 'nothing'));
      if (new URLSearchParams(location.search).get('spacetime') !== ${JSON.stringify(id)}) out.push('the address reads ' + location.search);
      return out;
    })()`);
  }
  if (!found) errors.push(`${opened}: never ready`);
  for (const m of found || []) errors.push(`${opened}: ${m}`);
  await pageErrors();
}

async function pageErrors() {
  for (const e of await evaluate(`window.__mfsTiming.errors.splice(0)`)) errors.push(`${opened}: page error: ${e}`);
}
await pageErrors();

const results = [];
for (const id of ids) {
  if (only.length && !only.includes(id)) continue;
  opened = id;
  const first = await timed(opening(id), 0);
  results.push({ id, ...(first || { chart: '?' }), failed: !first });
  if (!first) continue;
  for (const m of first.moved.concat(first.margins)) errors.push(`${opened} / ${first.chart}: ${m}`);
  let conventions = await convention(id, 0);
  await contents();
  await pressed();
  await views();
  if (print) await printCopy();
  const charts = await evaluate(`document.querySelectorAll('#mfs-content-panel .mfs-charts .mfs-choice').length || 1`);
  for (let c = 1; c < charts; c++) {
    const r = await timed(`document.querySelector('#mfs-content-panel .mfs-charts [data-chart="${c}"]').click();`, c);
    results.push({ id, ...(r || { chart: String(c) }), failed: !r });
    if (!r) continue;
    for (const m of r.moved.concat(r.margins)) errors.push(`${opened} / ${r.chart}: ${m}`);
    conventions = await convention(id, c, conventions);
    await views();
    if (print) await printCopy();
  }
  await pageErrors();
}

const last = results.filter(r => !r.failed).pop();
if (last) await direct(last.id);

results.sort((a, b) => (a.failed - b.failed) || (a.seconds - b.seconds));
let bad = 0;
for (const r of results) {
  const over = r.failed || r.seconds > budget;
  if (over) bad++;
  const name = `${r.id} / ${r.chart}`.padEnd(44);
  console.log(r.failed
    ? `${name}  never ready within ${LIMIT_S}s`
    : `${name} ${r.seconds.toFixed(2).padStart(6)}s  longest task ${r.longest.toFixed(2).padStart(5)}s  ` +
      `${String(r.lines).padStart(4)} lines ${String(r.elements).padStart(7)} elements${over ? '  OVER BUDGET' : ''}`);
}
for (const e of errors) console.log(e);
for (const selector of BUTTONS) console.log(`No spacetime showed ${selector} to press.`);
console.log(`${results.length} charts at ${width} by ${height}${phone ? ' as a phone' : ''}` +
  `${text !== 16 ? ` with ${text}px text` : ''}${still ? ' with reduced motion' : ''}${print ? ' and their print copies' : ''}${cpu > 1 ? `, processor slowed ${cpu} times` : ''}: ` +
  (bad ? `${bad} over the ${budget}s budget` : `all within the ${budget}s budget`) +
  `, ${errors.length ? `${errors.length} errors` : 'no errors'}`);
finish(bad || errors.length ? 1 : 0);

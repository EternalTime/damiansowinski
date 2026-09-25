#!/usr/bin/env node
/* Typeset every TeX string the spacetimes page sets, and report each one MathJax cannot.

   sympy reads the mathematics and never the typesetting, so a published value can be right
   and still print as an error box. This sets, through the same TeX input jax the page loads:

     metrics      every published value and the same value with a minus sign in front,
                  which the page prints when it merges a component with its opposite, the
                  line elements, domains, geodesics and scalars, the signature, and the
                  mathematics inside every history, convention and parameter description;
     diagrams     every label, tick, caption and declared input of the null ray diagrams;
     conformal    every label, legend, caption, restriction, setting and input of the
                  conformal diagrams.

   It needs mathjax-full, which does not belong in the repository:

     npm install --prefix /tmp/mfs-node mathjax-full
     node _tools/derivations/tex_check.cjs /tmp/mfs-node

   and exits non-zero, naming the file, the place and the string, if anything fails. */
'use strict';
const fs = require('fs');
const path = require('path');
const {createRequire} = require('module');

const prefix = process.argv[2] || '/tmp/mfs-node';
const need = createRequire(path.join(path.resolve(prefix), 'node_modules', 'mathjax-full', 'package.json'));
const {mathjax} = need('mathjax-full/js/mathjax.js');
const {TeX} = need('mathjax-full/js/input/tex.js');
const {CHTML} = need('mathjax-full/js/output/chtml.js');
const {liteAdaptor} = need('mathjax-full/js/adaptors/liteAdaptor.js');
const {RegisterHTMLHandler} = need('mathjax-full/js/handlers/html.js');
need('mathjax-full/js/input/tex/AllPackages.js');

// The page loads tex-mml-chtml.js, whose TeX input carries these packages and autoloads the
// rest; noundefined is what turns an unknown command into the red text caught below.
const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const tex = new TeX({packages: ['base', 'ams', 'newcommand', 'noundefined', 'require', 'autoload', 'configmacros']});
const doc = mathjax.document('', {InputJax: tex, OutputJax: new CHTML()});

const ROOT = path.resolve(__dirname, '..', '..');
const DATA = path.join(ROOT, 'MFS', 'assets', 'data');

// Mathematics set as it stands, and prose whose mathematics sits between dollar signs.
const MATH = new Set(['line_element', 'value', 'ricci_scalar', 'kretschmann', 'signature', 'domains', 'geodesics']);
const PROSE = new Set(['history', 'convention', 'description', 'caption', 'input', 'settings', 'restriction',
                       'label', 'xlabel', 'ylabel', 'legend', 'text', 'families']);

function segments(text) {
  const out = [];
  const re = /\$\$([\s\S]+?)\$\$|\$([^$]+)\$/g;
  let m;
  while ((m = re.exec(text))) out.push(m[1] !== undefined ? m[1] : m[2]);
  return out;
}

const items = [];

function walk(value, where, key, file) {
  if (typeof value === 'string') {
    if (MATH.has(key)) {
      const body = value.replace(/^\$|\$$/g, '');
      items.push([file, where, body]);
      if (key === 'value') items.push([file, where + ' negated', '-' + body]);
    } else if (PROSE.has(key)) {
      segments(value.replace(/<[^>]*>/g, '')).forEach(function(s) { items.push([file, where, s]); });
    }
    return;
  }
  if (Array.isArray(value)) {
    value.forEach(function(v, i) { walk(v, where + '[' + i + ']', key, file); });
    return;
  }
  if (value && typeof value === 'object') {
    for (const k of Object.keys(value)) walk(value[k], where + '.' + k, k, file);
  }
}

for (const folder of ['metrics', 'diagrams', 'conformal']) {
  const dir = path.join(DATA, folder);
  if (!fs.existsSync(dir)) continue;
  for (const name of fs.readdirSync(dir).sort()) {
    if (!name.endsWith('.json') || / \d+\.json$/.test(name)) continue;
    walk(JSON.parse(fs.readFileSync(path.join(dir, name), 'utf8')), '', '', folder + '/' + name);
  }
}

let bad = 0;
for (const [file, where, source] of items) {
  const html = adaptor.outerHTML(doc.convert(source, {display: false}));
  if (/<mjx-merror|data-mjx-error|mathcolor="red"|color:\s*red/i.test(html)) {
    bad++;
    const reason = (html.match(/data-mjx-error="([^"]*)"/) || [])[1] || 'an undefined command';
    console.error(`${file} ${where}: ${reason}\n    ${source}`);
  }
}
console.log(`${items.length} TeX strings set, ${bad} failed`);
process.exit(bad ? 1 : 0);

/* A BibTeX field from assets/data/references.bib as a reader reads it, the one reader the
   spacetimes page (_layouts/mfs.html), the publications page (publications.markdown) and the
   applets' notes (assets/js/applet-shell.js) all include, and which
   _tools/build_agent_data.py's unlatex follows for the agent layer.

   Mathematics between dollar signs is left exactly as written, for MathJax to set on the page.
   Until 16 October 2026 the braces were taken off inside it too, so $dS/\mathbb{Z}_2$ reached
   MathJax as \mathbbZ_2 and was printed as an error in red, and $10^{-13}$ came out 10^-13.

   Outside it, a TeX accent is put on whatever letter it stands over, as \'e, \'{e}, {\'e}, \c{c}
   or \v c are written: the accent becomes Unicode's combining mark after the letter and the two
   are composed into one character. Until 15 October 2026 each accent had a table of the vowels
   it went on, so Stojkovi\'c, Bizo\'n, Chru\'sciel and voja\^go were printed as typed, and
   before then \'e was read as an umlaut, so Poincaré came out Poincarë. \i and \j are the
   dotless letters TeX puts an accent on, so \'{\i} is í. \ss is ß and \& is &, TeX's dashes
   -- and --- are both an en dash, since the site carries no em dash, `` and '' are its double
   quotes and ` its single opening quote, which are written as the plain quotes the rest of a
   title carries, and the braces that keep a word's capitals are taken off. */
var UNLATEX_MARKS = {
  "'": '́', '"': '̈', '`': '̀', '^': '̂', '~': '̃', '=': '̄',
  '.': '̇', u: '̆', v: '̌', H: '̋', c: '̧', k: '̨', r: '̊'
};
function unlatex(str) {
  return str.split(/(\$[^$]*\$)/).map(function(part, i) {
    if (i % 2) return part;
    return part
      .replace(/\\(?:(['"`^~=.])|([uvHckr])(?=\s*\{|\s+[A-Za-z\\]))\s*(?:\{\s*(\\[ij](?![A-Za-z])|[A-Za-z])\s*\}|(\\[ij](?![A-Za-z])|[A-Za-z]))/g,
        function(_, sign, word, braced, bare) {
          var letter = (braced || bare).replace(/^\\/, '');
          return (letter + UNLATEX_MARKS[sign || word]).normalize('NFC');
        })
      .replace(/\\ss\b/g, 'ß')
      .replace(/\\&/g, '&')
      .replace(/---?/g, '\u2013')
      .replace(/``|''/g, '"')
      .replace(/`/g, "'")
      .replace(/\{([^}]*)\}/g, '$1');
  }).join('').trim();
}

---
layout: default
title: Publications
permalink: /publications/
---

<a href="{{ '/' | relative_url }}" style="float: right;">← Back to Home</a>
<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Publications (with no pay-walls!)
</h1>

---
<br>
<div id="publications-list"></div>

<script>
(function() {
  var KEYS = [
    'sowinski2026beh',
    'sowinski2025exo',
    'sowinski2025causal',
    'sowinski2025configurational',
    'sowinski2024information',
    'sowinski2023semantic',
    'sowinski2022consensus',
    'sowinski2021poroelasticity',
    'sowinski2018instantons',
    'sowinski2017information',
    'sowinski2016complexity',
    'pinero2026information',
    'quillen2025notions',
    'mcgarry2025vivo',
    'mcgarry2022mapping',
    'mcgarry2021heterogenous',
    'hannum2022correlated',
    'jyoti2022quantifying',
    'bowen2022visual',
    'gleiser2018oscillons',
    'gleiser2018mapmaking',
    'gleiser2015information',
    'gleiser2013information'
  ];

  function unlatex(str) {
    return str
      .replace(/\\["']([aeiouAEIOUy])/g, function(_, c) {
        return { a:'ä',e:'ë',i:'ï',o:'ö',u:'ü',A:'Ä',E:'Ë',I:'Ï',O:'Ö',U:'Ü',y:'ÿ' }[c] || c;
      })
      .replace(/\\'([aeiouAEIOUy])/g, function(_, c) {
        return { a:'á',e:'é',i:'í',o:'ó',u:'ú',A:'Á',E:'É',I:'Í',O:'Ó',U:'Ú',y:'ý' }[c] || c;
      })
      .replace(/\\`([aeiouAEIOU])/g, function(_, c) {
        return { a:'à',e:'è',i:'ì',o:'ò',u:'ù',A:'À',E:'È',I:'Ì',O:'Ò',U:'Ù' }[c] || c;
      })
      .replace(/\\~([nNaAoO])/g, function(_, c) {
        return { n:'ñ',N:'Ñ',a:'ã',A:'Ã',o:'õ',O:'Õ' }[c] || c;
      })
      .replace(/\\c\{?([cCsS])\}?/g, function(_, c) {
        return { c:'ç',C:'Ç',s:'ş',S:'Ş' }[c] || c;
      })
      .replace(/\\ss\b/g, 'ß')
      .replace(/\{([^}]*)\}/g, '$1')
      .trim();
  }

  function parseBibtex(entry) {
    var typeKey = entry.match(/^@(\w+)\s*\{\s*([\w:]+)\s*,/);
    var type = typeKey ? typeKey[1].toLowerCase() : 'misc';
    var fields = {};
    var fieldRe = /(\w+)\s*=\s*(?:\{([^}]*(?:\{[^}]*\}[^}]*)*)\}|"([^"]*)")/g;
    var m;
    while ((m = fieldRe.exec(entry)) !== null) {
      fields[m[1].toLowerCase()] = unlatex(m[2] !== undefined ? m[2] : m[3]);
    }
    return { type: type, fields: fields };
  }

  function formatEntry(entry) {
    var p = parseBibtex(entry);
    var f = p.fields;
    var author  = f.author  || '';
    var title   = f.title   || '';
    var year    = f.year    || '';
    var journal = f.journal || f.booktitle || f.publisher || '';
    var volume  = f.volume  || '';
    var number  = f.number  || '';
    var pages   = f.pages   || '';
    var url     = f.url     || '';
    var html = '';
    if (author)  html += '<span class="pub-author">'  + author  + '</span> ';
    if (year)    html += '<span class="pub-year">('   + year    + ')</span>. ';
    if (title)   html += '<span class="pub-title">'   + title   + '</span>. ';
    if (journal) html += '<span class="pub-journal">' + journal + '</span>';
    if (volume)  html += ', <span class="pub-volume">' + volume + '</span>';
    if (number)  html += '(' + number + ')';
    if (pages)   html += ', ' + pages;
    if (url)     html += ' <a href="' + url + '" target="_blank" class="pub-pdf">[PDF]</a>';
    return html;
  }

  fetch('{{ "/assets/data/references.bib" | relative_url }}')
    .then(function(r) { return r.text(); })
    .then(function(text) {
      var bibCache = {};
      text.split(/(?=@\w+\s*\{)/).forEach(function(entry) {
        entry = entry.trim();
        if (!entry) return;
        var m = entry.match(/^@\w+\s*\{\s*([\w:]+)\s*,/);
        if (m) bibCache[m[1]] = entry;
      });
      var output = KEYS.map(function(key) {
        var entry = bibCache[key];
        if (!entry) return '<div class="csl-entry pub-missing">' + key + '</div>';
        return '<div class="csl-entry">' + formatEntry(entry) + '</div>';
      }).join('');
      document.getElementById('publications-list').innerHTML = output;
    })
    .catch(function(e) {
      console.error('Failed to load references.bib', e);
      document.getElementById('publications-list').innerHTML = '<p>Failed to load publications.</p>';
    });
})();
</script>

<style>
  .csl-entry {
    margin-bottom: 1.4em;
    line-height: 1.6;
  }
  .pub-author  { font-weight: normal; }
  .pub-year    { opacity: 0.65; }
  .pub-title   { font-style: italic; }
  .pub-journal { opacity: 0.8; }
  .pub-volume  { font-weight: normal; }
  .pub-pdf     { font-size: 0.75em; font-weight: bold; margin-left: 0.3em; }
  .pub-missing { opacity: 0.5; font-style: italic; }
</style>

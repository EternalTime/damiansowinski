#!/usr/bin/env python3
"""The rules a pattern can hold of the prose a reader of the spacetimes page meets.

`text_problems` holds texts to them, and `test_build_mfs_data.py` holds every text on the site to
them. _tools/README.md carries the rules these stand for.
"""

import re

# Words that in this prose can only mean the collection's own machinery or the page
# describing itself, which a reader who came for a spacetime is never told about.
# _tools/README.md carries the rule this stands for.
MACHINERY = (
    r"(?i)\bentr(?:y|ies)\b",
    r"(?i)\b(?:this|the) collection\b",
    r"(?i)\b(?:displayed|published|printed|shown|listed) (?:here|above|below)\b",
    r"(?i)\bthe history describes\b",
    r"`",
)
MACHINERY_OUTSIDE_A_HISTORY = (
    r"(?i)\bpublish\w*",
    r"(?i)\bprint(?:ed|s)?\b",
    r"(?i)(?<!building )\bblocks?\b",
    r"(?i)\bvariants?\b",
)

# "is the distance the metric gives": a noun handed back by what gives it, which reads
# as artificial where "is the metric distance" says the same thing. The sentences and
# the templates say it the compact way. _tools/README.md carries the rule this stands for.
_WORD = r"(?:(?!(?:of|and|or|against|to|in|on|at|by)\b)[\w'-]+ )"
GIVES = rf"(?i)\bthe {_WORD}{{1,2}}the {_WORD}{{0,2}}(?:gives?|yields?|provides?|returns?)\b"
# Nor is the metric said to give anything: an equation is stated as it stands, "where $f = 1$,
# $ds^2 = -c^2dt^2$", never "where $f = 1$ and the metric gives $ds^2 = -c^2dt^2$".
METRIC_GIVES = r"(?i)\bthe metric gives?\b"

# The captain's voice, ~/VOICE.md, as far as a pattern can hold it, over every paragraph of every
# history and convention and every caption, note, restriction band and sentence stated in place
# of a drawing. Mathematics and quotations are left out first: a quotation is its speaker's own
# words. _tools/README.md carries the rules these stand for.
_NOUN = r"(?:(?!(?:of|and|or|against|to|in|on|at|by)\b)[\w'-]+ )"
VOICE = (
    ("an em dash", r"—"),
    ("a dash standing alone as punctuation", r"(?:^|\s)-(?:\s|$)"),
    ("a noun handed back to what gives it", rf"(?i)\bthe {_NOUN}{{1,2}}the {_NOUN}{{0,2}}"
                                            r"(?:gives?|yields?|provides?|returns?|affords?|supplies|supply"
                                            r"|produces?|delivers?|offers?|draws?|declares?)\b"),
    ("a paper, a result or a field cast as the actor",
     r"(?i)\b(?:papers?|articles?|literature|results?|stud(?:y|ies)|textbooks?|numerical relativity)\b "
     r"(?:reports?|shows?|landed|lands|argues?|finds?|reads?|tests?|demonstrates?|tells?|writes?|is written)\b"
     r"|\breads? as\b|\bwere read\b"),
    ("a drawing cast as the actor", r"(?i)\b(?:diagrams?|drawing|figure|caption|section|page)\b "
                                    r"(?:declares?|draws?|shows?|says|tells)\b"),
    ("a staged reveal", r"(?i)(?:^|[.;:,] |\b(?:and|but|so) )what\b[^.;:$]*?\b(?:is|are|was)\b"
                        r"|\bwhat comes out\b"),
    ("a sentence explaining what is not drawn",
     r"(?i)\b(?:is|are) (?:not|never) (?:drawn|shown|plotted)\b|\bnothing (?:is |are )?(?:drawn|shown|plotted)\b"
     r"|\b(?:drawn|shown) instead\b|\bno \w+(?: \w+)? (?:is|are) (?:drawn|shown|given|plotted)\b"),
    ("a thing said to be what it is not", r"(?i)\b(?:is|are) not\b[^.;:,]{1,80}\bbut\b"),
    ("an absolute for emphasis", r"(?i)\bevery single\b|\bat all\b|\bnever fails\b|\beach and every\b|\bwhatsoever\b"
                                 r"|\bactual(?:ly)?\b"),
    ("a machine's tell", r"(?i)\bit is important to note\b|\bit(?:'s| is) worth not(?:ing|e)\b|\bin essence\b"
                         r"|\binteresting question\b|\bdive in\b|\bfascinating\b"
                         r"|(?:^|\. )(?:Crucially|Importantly|Notably|Interestingly),|\bdelve\b|\butili[sz]e\b"
                         r"|\bleverage\b|\bbookkeeping\b|\bfunctionals?\b|\bin today's world\b|\bin modern times\b"
                         r"|\badditionally\b"),
    ("stacked transitions", r"(?i)\b(?:moreover|furthermore),? (?:moreover|furthermore)\b"),
    ("a sentence about the page", r"(?i)\bas (?:noted|shown|mentioned|described|stated|discussed) "
                                  r"(?:above|below|earlier)\b|\bthis (?:section|caption|note|page|entry)\b"),
    ("a slide into the future tense", r"(?i)\bwe will\b|\bwe'll\b"),
    ("a list without its Oxford comma", r"\$M\$, \$M\$ (?:and|or) \$M\$"),
    ("a figure of speech standing in for the claim",
     r"(?i)\bstands? on its own\b|\bhas its say\b|\bkeeps? its secrets\b|\ba question with answers\b"),
)
# A caption, note or restriction band is a compact noun phrase with its values in parentheses, as in
# the captain's own "A spherically symmetric distribution of dust collapsing from rest ($R_0 = 2\,r_s$),
# each point in the diagram a 2-sphere." It never opens "This is", never says "This is the whole of",
# and nothing in it "stands for" anything. _tools/README.md carries the rule these stand for.
CAPTION_VOICE = (
    ("\"This is the whole of\"", r"(?i)\bthis is the whole of\b"),
    ("\"stands for a sphere\"", r"(?i)\bstands? for (?:a|one|one such) (?:\w+ )?sphere\b"),
    ("\"stands for\"", r"(?i)\bstands? for\b"),
    ("a sentence opening \"This is\" in place of a noun phrase", r"^This is\b"),
)
# A hyphen joins two names, a name and a word, or a designation; beyond those it is part of a
# spelling only in these terms, and an ordinary compound is rewritten without it. A term goes on the
# list only where its standard spelling carries the hyphen, as the no-boundary proposal and e-foldings
# do; a sphere or a geometry of n dimensions is written with its numeral, 3-sphere, 4-sphere, 3-geometry,
# as a designation, and needs no entry. The rule is for compounds such as "future-directed",
# "non-singular" and "time-like", which plain English writes as one word or two.
HYPHENATED_TERMS = {"anti-de", "anti-trapped", "plane-fronted", "pp-wave", "pp-waves", "scalar-tensor",
                    "three-brane", "three-branes", "no-boundary", "e-folding", "e-foldings", "s-wave",
                    "super-entropic", "quasi-topological"}


def without_mathematics_or_quotations(text):
    return re.sub(r'"[^"]*"', '"Q"', re.sub(r"\$\$.+?\$\$|\$[^$]+\$", "$M$", text, flags=re.DOTALL))


def refused_hyphens(text):
    """Every hyphenated word of a text, mathematics, quotations and citations taken out, that joins
    no name or designation and is no established term, a term's possessive read as the term."""
    text = re.sub(r"\[[^\]]*\]", " ", without_mathematics_or_quotations(text))
    return [word for word in re.findall(r"[\w'’]+(?:-[\w'’]+)+", text)
            if not any(part[0].isupper() or part[0].isdigit() for part in word.split("-"))
            and re.sub(r"['’]s$", "", word.lower()) not in HYPHENATED_TERMS]


# A thing defined by a contrast with what it is not: "X rather than Y", "instead of", "not a X but Y",
# "is X, not Y" or "less X, more Y".
CONTRAST = (r"(?i)\brather than\b|\binstead\b"
            r"|\bnot (?:a|an|the|just|only|merely|simply)\b[^.;:]{0,80}?\bbut\b"
            r"|\b(?:is|are|was|were)\b[^.;:,]{1,60},\s*not (?!long\b|yet\b)\w"
            r"|\bless\b[^.;:,]{1,40},\s*more\b")

# Physics said to require something in the words of money.
COST = r"(?i)\bcost(?:s|ly|ing)?\b"

# The rules a text may be held to, by name, each finding what in a text breaks it.
RULES = {
    "machinery": lambda text: [("a word for the collection's machinery", m)
                               for p in MACHINERY for m in [re.search(p, text)] if m],
    "machinery outside a history": lambda text: [("a word for the collection's machinery", m)
                                                 for p in MACHINERY_OUTSIDE_A_HISTORY for m in [re.search(p, text)] if m],
    "gives": lambda text: [("a noun handed back to what gives it", m)
                           for p in (GIVES, METRIC_GIVES) for m in [re.search(p, text)] if m],
    "voice": lambda text: [(rule, m) for rule, p in VOICE
                           for m in [re.search(p, without_mathematics_or_quotations(text))] if m],
    "caption": lambda text: [(rule, m) for rule, p in CAPTION_VOICE
                             for m in [re.search(p, without_mathematics_or_quotations(text))] if m],
    "hyphens": lambda text: [("a compound spelled with a hyphen", word) for word in refused_hyphens(text)],
    "contrast": lambda text: [("a thing defined by what it is not", m)
                              for m in re.finditer(CONTRAST, without_mathematics_or_quotations(text))],
    "cost": lambda text: [("physics said in the words of money", m) for m in [re.search(COST, text)] if m],
}


def text_problems(texts):
    """What is wrong with each text, given as its place, its words and the names of the RULES it is
    held to, each problem the place, the rule broken and the words that break it, in the order
    given; an empty list when every text keeps to its rules."""
    problems = []
    for where, text, rules in texts:
        for name in rules:
            for rule, found in RULES[name](text):
                words = found if isinstance(found, str) else found.group(0)
                problems.append(f"{where} has {rule}: {words!r}")
    return problems

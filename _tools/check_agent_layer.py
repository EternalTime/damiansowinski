#!/usr/bin/env python3
"""Check the agent layer of a built site: llms.txt, the JSON data, robots.txt and JSON-LD.

    bundle exec jekyll build
    python3 _tools/check_agent_layer.py _site
    python3 _tools/check_agent_layer.py _site --external   # also fetch every outside link

Every address on damiansowinski.com that llms.txt or llms-full.txt names must be a file in
the build, every outside address in llms.txt must answer when --external is given (those
in llms-full.txt come from the pages and are only reported), both JSON files must parse and carry what agents query, robots.txt must point at
a sitemap that parses, and each JSON-LD block must use only types and properties that the
schema.org vocabulary defines, each property on a type it belongs to. The vocabulary is
downloaded once and kept in ~/.cache/damiansowinski/.
"""

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import unquote, urlsplit

SITE_URL = "https://damiansowinski.com"
VOCABULARY_URL = "https://schema.org/version/latest/schemaorg-current-https.jsonld"
VOCABULARY_CACHE = Path.home() / ".cache" / "damiansowinski" / "schemaorg-current-https.jsonld"
MARKDOWN_LINK = re.compile(r"\]\(([^)\s]+)\)")
JSON_LD = re.compile(r'<script type="application/ld\+json">\s*(.*?)\s*</script>', re.DOTALL)
# Pages whose JSON-LD the agent layer adds, and the types each must describe.
JSON_LD_PAGES = {"index.html": "Person", "publications/index.html": "CollectionPage", "MFS/index.html": "Dataset"}
# Sites that answer every script with a refusal, whatever the page.
BOT_WALLS = {"www.linkedin.com": {999}, "x.com": {400, 403}, "doi.org": {403}}


class Problems(list):
    def check(self, condition, message):
        if not condition:
            self.append(message)


def local_file(site, url):
    path = unquote(urlsplit(url).path)
    target = site / path.lstrip("/")
    return target / "index.html" if path.endswith("/") else target


def check_links(site, name, problems, external):
    text = (site / name).read_text(encoding="utf-8")
    links = MARKDOWN_LINK.findall(text) + [u.rstrip(".,;:") for u in re.findall(r"(?<=\s)(https://[^\s)]+)", text)]
    outside = set()
    for url in dict.fromkeys(links):
        if url.startswith("mailto:"):
            continue
        if url.startswith(SITE_URL + "/") or url == SITE_URL:
            problems.check(local_file(site, url).is_file(), f"{name}: {url} is not in the build")
        else:
            problems.check(url.startswith("https://"), f"{name}: {url} is not an https address")
            outside.add(url)
    if external:
        with ThreadPoolExecutor(max_workers=16) as pool:
            statuses = dict(zip(sorted(outside), pool.map(fetch_status, sorted(outside))))
        for url, status in statuses.items():
            if status < 400 or status in BOT_WALLS.get(urlsplit(url).hostname, set()):
                continue
            message = f"{name}: {url} answered {status}"
            # llms.txt is written for agents and must hold; llms-full.txt carries the pages'
            # own outside links, which are the pages' business and are only reported.
            if name == "llms.txt":
                problems.append(message)
            else:
                print(f"warning: {message}", file=sys.stderr)
    return len(links)


def fetch_status(url):
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (link check)"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return response.status
    except urllib.error.HTTPError as exc:
        return exc.code
    except (urllib.error.URLError, TimeoutError):
        return 599


def check_llms(site, problems, external):
    text = (site / "llms.txt").read_text(encoding="utf-8")
    lines = text.splitlines()
    problems.check(lines[0].startswith("# ") and not lines[0].startswith("## "),
                   "llms.txt: the first line must be the site's name as an H1")
    problems.check(any(line.startswith("> ") for line in lines[:4]),
                   "llms.txt: a blockquote summary must follow the H1")
    sections = [line[3:] for line in lines if line.startswith("## ")]
    for section in ("Research", "Publications", "My Favorite Spacetimes", "Python libraries",
                    "Applets", "Notes", "Contact"):
        problems.check(section in sections, f"llms.txt: no '## {section}' section")
    for library in ("pyCA", "pyCE", "pyCoop", "pyEDW", "pyGD", "LEAFS"):
        problems.check(f"[{library}: " in text, f"llms.txt: the {library} documentation is not linked")
    for data in ("/data/spacetimes.json", "/data/publications.json", "/llms-full.txt"):
        problems.check(f"({SITE_URL}{data})" in text, f"llms.txt: {data} is not linked")
    problems.check("{{" not in text and "{%" not in text, "llms.txt: Liquid was left unrendered")
    count = check_links(site, "llms.txt", problems, external)
    count += check_links(site, "llms-full.txt", problems, external)
    return count


def check_data(site, problems):
    spacetimes = json.loads((site / "data" / "spacetimes.json").read_text(encoding="utf-8"))
    problems.check(spacetimes["count"] == len(spacetimes["spacetimes"]) > 0, "spacetimes.json: bad count")
    for s in spacetimes["spacetimes"]:
        for field in ("id", "name", "description", "tags", "signature", "page_url", "data_url", "charts"):
            problems.check(s.get(field), f"spacetimes.json: {s.get('id')} has no {field}")
        problems.check(local_file(site, s["data_url"]).is_file(), f"spacetimes.json: {s['data_url']} is missing")
        for chart in s["charts"]:
            for field in ("coordinates", "line_element", "metric_components"):
                problems.check(chart.get(field), f"spacetimes.json: {s['id']}/{chart.get('id')} has no {field}")
    publications = json.loads((site / "data" / "publications.json").read_text(encoding="utf-8"))
    problems.check(publications["count"] == len(publications["publications"]) > 0, "publications.json: bad count")
    for p in publications["publications"]:
        for field in ("title", "authors", "year", "venue", "url"):
            problems.check(p.get(field), f"publications.json: {p.get('key')} has no {field}")
        if p.get("pdf_url", "").startswith(SITE_URL):
            problems.check(local_file(site, p["pdf_url"]).is_file(), f"publications.json: {p['pdf_url']} is missing")
    return len(spacetimes["spacetimes"]), len(publications["publications"])


def check_robots(site, problems):
    robots = (site / "robots.txt").read_text(encoding="utf-8")
    problems.check(re.search(r"(?m)^User-agent: \*\nAllow: /$", robots), "robots.txt: '*' is not allowed /")
    problems.check(not re.search(r"(?mi)^Disallow: *\S", robots), "robots.txt: something is disallowed")
    sitemap = re.search(r"(?m)^Sitemap: (\S+)$", robots)
    problems.check(sitemap, "robots.txt: no Sitemap line")
    if sitemap:
        path = local_file(site, sitemap.group(1))
        problems.check(path.is_file(), f"robots.txt: {sitemap.group(1)} is not in the build")
        if path.is_file():
            ET.parse(path)


def load_vocabulary():
    if not VOCABULARY_CACHE.exists():
        VOCABULARY_CACHE.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(VOCABULARY_URL, timeout=60) as response:
            VOCABULARY_CACHE.write_bytes(response.read())
    graph = json.loads(VOCABULARY_CACHE.read_text(encoding="utf-8"))["@graph"]

    def names(value):
        values = value if isinstance(value, list) else [value] if value else []
        return {v["@id"].removeprefix("schema:") for v in values}

    classes, properties = {}, {}
    for node in graph:
        kind = node["@type"] if isinstance(node["@type"], list) else [node["@type"]]
        name = node["@id"].removeprefix("schema:")
        if "rdfs:Class" in kind:
            classes[name] = names(node.get("rdfs:subClassOf"))
        if "rdf:Property" in kind:
            properties[name] = names(node.get("schema:domainIncludes"))
    return classes, properties


def ancestors(classes, name):
    seen, stack = set(), [name]
    while stack:
        current = stack.pop()
        if current not in seen:
            seen.add(current)
            stack.extend(classes.get(current, ()))
    return seen


def check_node(node, where, vocabulary, problems):
    classes, properties = vocabulary
    if isinstance(node, list):
        for item in node:
            check_node(item, where, vocabulary, problems)
        return
    if not isinstance(node, dict):
        return
    kind = node.get("@type")
    if kind is None:
        return
    problems.check(kind in classes, f"{where}: {kind} is not a schema.org type")
    lineage = ancestors(classes, kind)
    for key, value in node.items():
        if key.startswith("@"):
            continue
        if key not in properties:
            problems.append(f"{where}: {key} is not a schema.org property")
            continue
        problems.check(properties[key] & lineage, f"{where}: {key} does not belong on a {kind}")
        check_node(value, f"{where} > {key}", vocabulary, problems)


def check_json_ld(site, problems):
    vocabulary = load_vocabulary()
    checked = 0
    for page, expected in JSON_LD_PAGES.items():
        blocks = [json.loads(b) for b in JSON_LD.findall((site / page).read_text(encoding="utf-8"))]
        types = [b.get("@type") for b in blocks]
        problems.check(expected in types, f"{page}: no {expected} JSON-LD")
        for block in blocks:
            problems.check(block.get("@context") == "https://schema.org", f"{page}: JSON-LD without the schema.org context")
            check_node(block, f"{page} {block.get('@type')}", vocabulary, problems)
            checked += 1
    return checked


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("site", type=Path, help="the built site, such as _site")
    parser.add_argument("--external", action="store_true", help="fetch every link outside the site")
    args = parser.parse_args(argv)

    problems = Problems()
    links = check_llms(args.site, problems, args.external)
    spacetimes, publications = check_data(args.site, problems)
    check_robots(args.site, problems)
    blocks = check_json_ld(args.site, problems)
    for problem in problems:
        print(f"problem: {problem}", file=sys.stderr)
    if problems:
        return 1
    print(f"{links} links, {spacetimes} spacetimes, {publications} publications and "
          f"{blocks} JSON-LD blocks check out")
    return 0


if __name__ == "__main__":
    sys.exit(main())

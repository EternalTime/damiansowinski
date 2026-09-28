# Damian Sowinski

## The agent layer

The site carries a layer for AI agents that no visitor sees: nothing in it is linked from a page or drawn on one.

- `/llms.txt` is the index for language models in the [llms.txt](https://llmstxt.org) format, built by Jekyll from `llms.txt` at the top of the repository.
- `/llms-full.txt` is the text of every main page as Markdown in one file.
  One file was chosen over a Markdown twin beside each page because GitHub Pages turns any `.md` file into HTML, and because one generated file is simpler to keep in step.
- `/data/spacetimes.json` is the My Favorite Spacetimes catalogue and `/data/publications.json` the publications, both published by Jekyll from `_data/generated/`.
- `/robots.txt` welcomes every crawler and AI agent and points at `/sitemap.xml`, which `jekyll-sitemap` builds.
- The home page carries schema.org JSON-LD for a Person, the publications page an ItemList of ScholarlyArticles and the spacetimes page a Dataset, from `_includes/agent-json-ld.html`.

Everything in it comes from the sources the pages themselves read, and nothing is copied by hand.
The publications page reads its order from `_data/publication_keys.json` and the code page its libraries from `_data/libraries.json`; `_data/person.json` holds who Damian is.
`llms.txt` and the JSON-LD are rendered from those and from `_data/generated/` on every build.

`_data/generated/` and `llms-full.txt` are written by one command, run from the top of the repository after a page, a paper, a library or a spacetime changes:

    python3 _tools/build_agent_data.py

`python3 -m unittest discover -s _tools` fails while either is out of date, as `python3 _tools/build_agent_data.py --check` does.
To check a build the way an agent reads it, including that every address in `llms.txt` answers and that the JSON-LD uses only what schema.org defines:

    bundle exec jekyll build
    python3 _tools/check_agent_layer.py _site --external

A spacetime has no address of its own, since the catalogue chooses one in the browser, so each entry gives the catalogue as its page and its metric file as its data.

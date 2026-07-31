#!/usr/bin/env python3
"""Phase-5: assemble index_of_entities.html from the template + the catalog.

The catalog is inlined as a <script type="application/json"> block, so the
result is one self-contained file that opens from file:// with no network
access, no build step and no dependencies — the same constraint the scene
window in posim/src/scene/scene.html works under.

Stdlib only.
"""

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(ROOT, "tools", "app_template.html")
CATALOG = os.path.join(ROOT, "index_data", "catalog.json")
OUT = os.path.join(ROOT, "index_of_entities.html")


def main():
    tpl = open(TEMPLATE, encoding="utf-8").read()
    data = json.load(open(CATALOG, encoding="utf-8"))

    # Minified, and with the two sequences that can break out of a <script>
    # block neutralised. `</script` is the real one; `<!--` is the historical
    # one, and both are cheap to rule out.
    blob = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    blob = blob.replace("</", "<\\/").replace("<!--", "<\\!--")

    if "/*__CATALOG__*/" not in tpl:
        sys.exit("template is missing the /*__CATALOG__*/ placeholder")
    html = "<!doctype html>\n<html lang=\"en\">\n" + tpl.replace("/*__CATALOG__*/", blob) \
        + "\n</html>\n"

    open(OUT, "w", encoding="utf-8").write(html)
    entries = [e for e in data if not e.get("_meta")]
    print("%s  —  %.0f KB, %d entries, %d examples"
          % (OUT, len(html) / 1024, len(entries),
             sum(len(e["examples"]) for e in entries)))


if __name__ == "__main__":
    main()

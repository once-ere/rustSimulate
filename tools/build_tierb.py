#!/usr/bin/env python3
"""Phase-7: Tier-B entries — the first-party Rust API.

Every public item of physical_object, special_functions, quantum and posim
becomes an entry carrying its real signature, its doc comment, and its
file:line. Examples are TEMPLATED ONLY where the signature shape can be
turned into Rust that is certain to compile; everything else is an honest
stub. Generating a plausible-looking snippet that does not build would be
worse than admitting there isn't one — the whole index rests on the claim
that its examples run.

Coverage is reported, not hidden: the script prints how many entries got
examples and how many did not, and the status page counts the stubs.

Stdlib only.
"""

import json
import re
import sys
from collections import Counter, defaultdict

TIER_B = ["physical_object", "special_functions", "quantum", "posim"]
LEVELS = ["trivial", "intermediate", "advanced", "expert"]

# --- how to build a receiver of each type, and what it costs to set up ----
RECEIVER = {
    "physical_object": ("physical_object::physical_object::physical_object",
        "let subject = physical_object::new_point(0, 2.0, Vec3::new(1.0, 0.0, 0.0), "
        "Vec3::new(0.0, 3.0, 0.0));"),
    "PhysicalObjectSystem": ("physical_object::system::PhysicalObjectSystem",
        "let subject = PhysicalObjectSystem::new(\n"
        "    vec![physical_object::new_point(0, 1.0, Vec3::new(1.0, 0.0, 0.0), "
        "Vec3::new(0.0, 1.0, 0.0))],\n    1.0,\n);"),
    "Vec3": ("physical_object::linalg::Vec3", "let subject = Vec3::new(1.0, 2.0, 3.0);"),
    "Mat3": ("physical_object::linalg::Mat3", "let subject = Mat3::identity();"),
    "Quat": ("physical_object::linalg::Quat", "let subject = Quat::identity();"),
    "Complex64": ("special_functions::complex::Complex64",
        "let subject = Complex64::new(3.0, -4.0);"),
    "Grid": ("quantum::qm1d::Grid", "let subject = Grid::new(-8.0, 8.0, 60)?;"),
    "Grid2": ("quantum::qm2d::Grid2",
        "let subject = Grid2::new(-4.0, 4.0, 20, -4.0, 4.0, 20)?;"),
    "Grid3": ("quantum::qm3d::Grid3",
        "let subject = Grid3::new(-3.0, 3.0, 10, -3.0, 3.0, 10, -3.0, 3.0, 10)?;"),
}

# --- a value of each argument type, chosen to be inside every domain ------
ARG = {
    "f64": "0.5", "i32": "1", "usize": "1", "u8": "1", "bool": "true",
    "Vec3": "Vec3::new(1.0, 0.0, 0.0)", "Mat3": "Mat3::identity()",
    "Quat": "Quat::identity()", "C": "Complex64::new(1.0, 0.5)",
    "Complex64": "Complex64::new(1.0, 0.5)",
    "&[f64]": "&[1.0, 2.0, 3.0]", "&Vec3": "&Vec3::new(1.0, 0.0, 0.0)",
    "&Mat3": "&Mat3::identity()",
}
# order arguments for the special-function families, whose first argument is
# an ORDER and whose second is the argument: (0, 2.0) is inside every domain
# these functions have, which (1, 0.5) is not for the singular ones.
SPECIAL_ARGS = {
    ("i32", "f64"): ("0", "2.0"),
    ("i32", "C"): ("0", "Complex64::new(2.0, 0.5)"),
    ("f64", "C"): ("0.5", "Complex64::new(2.0, 0.5)"),
    ("f64", "f64"): ("0.5", "2.0"),
    ("C", "C"): ("Complex64::new(0.5, 0.0)", "Complex64::new(2.0, 0.5)"),
    ("usize", "f64"): ("4", "2.0"),
}

SIG = re.compile(r"^pub (?:const )?fn (\w+)\s*\((.*?)\)\s*(?:->\s*(.+?))?\s*$", re.S)


def parse(sig):
    m = SIG.match(" ".join(sig.split()))
    if not m:
        return None
    name, args, ret = m.group(1), m.group(2).strip(), (m.group(3) or "()").strip()
    takes_self = args.startswith("&self") or args.startswith("&mut self") \
        or args.startswith("self")
    mut_self = args.startswith("&mut self")
    if takes_self:
        args = args.split(",", 1)[1] if "," in args else ""
    types = []
    for part in re.split(r",(?![^<(\[]*[>)\]])", args):
        part = part.strip().rstrip(",")
        if not part:
            continue
        if ":" not in part:
            return None
        types.append(part.split(":", 1)[1].strip())
    return {"name": name, "self": takes_self, "mut": mut_self,
            "args": types, "ret": ret}


def call_args(p, crate):
    if crate == "special_functions" and len(p["args"]) == 2:
        key = tuple(p["args"])
        if key in SPECIAL_ARGS:
            return list(SPECIAL_ARGS[key])
    out = []
    for t in p["args"]:
        if t not in ARG:
            return None
        out.append(ARG[t])
    return out


def snippet(item, p):
    """Rust that compiles, or None when the shape is not one we can be sure of."""
    crate, name = item["crate"], item["bare"]
    owner = item["name"].split("::")[0] if "::" in item["name"] else None
    args = call_args(p, crate)
    if args is None:
        return None
    fallible = p["ret"].startswith("Result<")
    q = "?" if fallible else ""
    if p["self"]:
        if owner not in RECEIVER:
            return None
        _, ctor = RECEIVER[owner]
        if p["mut"]:
            ctor = ctor.replace("let subject", "let mut subject")
        body = ctor + f"\nlet value = subject.{name}({', '.join(args)}){q};"
    else:
        if owner:                       # an associated fn that is not a method
            return None
        body = f"let value = {name}({', '.join(args)}){q};"
    if p["ret"] in ("()", "Result<(), String>"):
        body = body.replace("let value = ", "")
    else:
        # A binding, not a println: `Airy` and `Uniform` are plain data types
        # with no Debug, and a snippet that only compiles for the types that
        # happen to derive it is a snippet that lies about the rest.
        body = body.replace("let value = ", "let _value = ")
    return body


def use_lines(item, p):
    """Import lines for one snippet.

    Every path is written `::crate::...`. That is not decoration: the union
    struct is named `physical_object` in lower case by specification, so the
    moment a scope imports it the bare name `physical_object` resolves to the
    STRUCT and every later `physical_object::…` path fails with "is a struct,
    not a module". CLAUDE.md hard rule #8 says to prefix the others with `::`,
    and that is exactly what this does.
    """
    owner = item["name"].split("::")[0] if "::" in item["name"] else None
    uses = set()
    if owner and owner in RECEIVER:
        uses.add("::" + RECEIVER[owner][0])
    if item["crate"] == "physical_object" or owner in ("Vec3", "Mat3", "Quat"):
        uses.add("::physical_object::linalg::{Vec3, Mat3, Quat}")
    if not p["self"]:
        # items declared in lib.rs sit at the CRATE ROOT — there is no
        # `::lib::` module to path through
        mod = "" if item["module"] in ("lib", "main") else item["module"] + "::"
        uses.add(f"::{item['crate']}::{mod}{item['bare']}")
    for t in p["args"] + [p["ret"]]:
        if "Complex64" in t or t == "C":
            uses.add("::special_functions::complex::Complex64")
        if "Vec3" in t or "Mat3" in t or "Quat" in t:
            uses.add("::physical_object::linalg::{Vec3, Mat3, Quat}")
    if owner == "physical_object" or (item["crate"] == "physical_object"
                                      and "physical_object" in item["signature"]):
        uses.add("::physical_object::physical_object::physical_object")
    # the PhysicalObjectSystem constructor builds a physical_object, so it
    # needs the struct in scope too
    if owner in RECEIVER and "physical_object::new_point" in RECEIVER[owner][1]:
        uses.add("::physical_object::physical_object::physical_object")
    # the braced linalg import already brings in Vec3/Mat3/Quat; importing one
    # of them again by name is E0252, not a harmless duplicate
    if "::physical_object::linalg::{Vec3, Mat3, Quat}" in uses:
        for t in ("Vec3", "Mat3", "Quat"):
            uses.discard("::physical_object::linalg::" + t)
    return sorted(uses)


def main():
    items = [json.loads(l) for l in open("index_data/rust_items.jsonl")]
    pub = [d for d in items if d["crate"] in TIER_B and d["visibility"] == "pub"
           and d["kind"] != "mod"]

    # Tier-A builtins with the same bare name are the notebook-reachable twin
    cat = json.load(open("index_data/catalog.json"))
    a_by_name = {e["name"]: e["id"] for e in cat
                 if not e.get("_meta") and e["kind"] == "builtin"}

    out, made, skipped = [], 0, Counter()
    for d in pub:
        kind = "function" if d["kind"] == "fn" else "type"
        owner = d["name"].split("::")[0] if "::" in d["name"] else None
        eid = f"rs.{d['crate']}.{d['module']}.{d['name'].replace('::', '.')}"
        examples = []
        if d["kind"] == "fn" and d["crate"] == "posim":
            # posim is a binary crate with no lib target, so nothing can
            # `use posim::…`. Its items are catalogued and cross-linked to the
            # notebook commands they implement, but they get no Rust snippet.
            skipped["posim is a bin crate (not importable)"] += 1
        elif d["kind"] == "fn":
            p = parse(d["signature"])
            if p is None:
                skipped["unparsed signature"] += 1
            else:
                code = snippet(d, p)
                if code:
                    uses = use_lines(d, p)
                    examples.append({
                        "level": "trivial", "medium": "rust",
                        "code": "\n".join("use " + u + ";" for u in uses) + "\n\n" + code,
                        "expected": None, "verified": None, "runner": "cargo build"})
                    made += 1
                else:
                    skipped["shape not templatable"] += 1
        else:
            skipped[d["kind"]] += 1

        see = []
        if d["bare"] in a_by_name:
            see.append(a_by_name[d["bare"]])
        doc = d["doc"] or ""
        out.append({
            "id": eid, "name": d["name"], "kind": kind, "tier": "B",
            "aliases": [], "indexKeys": [(d["bare"][0] if d["bare"] else "Σ").upper()],
            "summary": (doc.split(".")[0][:150] + "." if doc
                        else f"{d['kind']} `{d['bare']}` in `{d['crate']}::{d['module']}`."),
            "definition": (doc or "No doc comment in the source.")
                + f"  Declared in the `{d['module']}` module of the `{d['crate']}` crate.",
            "syntax": [d["signature"]],
            "parameters": [], "returns": None, "errors": [],
            "locations": [{"file": d["file"], "line": d["line"],
                           "role": "definition"}],
            "examples": examples, "seeAlso": see, "invariants": [],
            "status": "complete" if examples else "stub",
        })

    json.dump(out, open("index_data/entries_tierb.json", "w"), indent=1)
    print(f"{len(out)} Tier-B entries -> index_data/entries_tierb.json")
    print(f"  with a templated example: {made}")
    print(f"  stubs: {len(out) - made}")
    for k, n in skipped.most_common():
        print(f"    {n:4d}  {k}")


if __name__ == "__main__":
    main()

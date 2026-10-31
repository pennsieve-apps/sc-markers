#!/usr/bin/env python3
"""
Development tool (not run by the processor, needs network): fetches every work the app cites
from its source and rewrites the `references:` of CITATION.cff from what came back, so no
reference is written from memory.

  papers    Crossref (title, full author list in order, journal, year, volume, issue, pages),
            then PubMed E-utilities: esearch by DOI, esummary; title, first author and year
            must match, else no PMID is written and a warning is printed.
  datasets  DataCite (title, creators, year, publisher, licence) and the Pennsieve Discover API
            (contributors with ORCID, licence, version).

  python scripts/verify_citations.py           # print the check
  python scripts/verify_citations.py --cff     # also rewrite CITATION.cff's references
"""
from __future__ import annotations

import html
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

PAPERS = {
    "scanpy": "10.1186/s13059-017-1382-0",
    "anndata": "10.21105/joss.04371",
    "wilcoxon": "10.2307/3001968",
    "mann_whitney": "10.1214/aoms/1177730491",
    "bh": "10.1111/j.2517-6161.1995.tb02031.x",
    "pullin": "10.1186/s13059-024-03183-0",
    "squair": "10.1038/s41467-021-25960-2",
    "zimmerman": "10.1038/s41467-021-21038-1",
    "zhang": "10.1016/j.cels.2019.07.012",
    "heumos": "10.1038/s41576-023-00586-w",
    "matplotlib": "10.1109/MCSE.2007.55",
    "drg_c2_paper": "10.1093/brain/awaf370",
}
# Demo data: DOI -> Discover dataset id.
DATASETS = {"10.26275/cxuf-ftav": 475}


def get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "pennsieve-apps citation check (mailto:support@pennsieve.io)"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as fh:
                return json.load(fh)
        except Exception:  # noqa: BLE001
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"cannot fetch {url}")


def plain(text: str) -> str:
    text = html.unescape(re.sub(r"<[^>]+>", "", text or ""))
    return re.sub(r"\s+", " ", text).strip()


def norm(text: str) -> str:
    text = unicodedata.normalize("NFKD", plain(text)).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]", "", text)


def crossref(doi: str) -> dict:
    m = get("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]
    authors = []
    for a in m.get("author", []):
        if a.get("family"):
            authors.append({"family": plain(a["family"]), "given": plain(a.get("given", ""))})
        elif a.get("name"):
            authors.append({"name": plain(a["name"])})
    year = None
    for k in ("published-print", "issued", "published-online"):
        dp = (m.get(k) or {}).get("date-parts")
        if dp and dp[0] and dp[0][0]:
            year = dp[0][0]
            break
    return {"doi": doi, "type": "article", "title": plain((m.get("title") or [""])[0]), "authors": authors,
            "journal": plain((m.get("container-title") or [""])[0]), "year": year, "volume": m.get("volume", ""),
            "issue": m.get("issue", ""), "pages": m.get("page") or m.get("article-number") or ""}


def pubmed(ref: dict) -> str:
    q = urllib.parse.quote(f"{ref['doi']}[doi]")
    for attempt in range(5):   # E-utilities answers an error body without idlist when rate-limited
        res = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}&retmode=json")
        time.sleep(0.6 * (attempt + 1))
        if "idlist" in res.get("esearchresult", {}):
            break
    else:
        raise RuntimeError(f"PubMed esearch failed for {ref['doi']}")
    ids = res["esearchresult"]["idlist"]
    if len(ids) != 1:
        return ""
    s = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={ids[0]}&retmode=json")["result"][ids[0]]
    time.sleep(0.4)
    title_ok = norm(s.get("title", ""))[:60] == norm(ref["title"])[:60]
    years = {str(ref["year"]), str(ref["year"] - 1)}
    year_ok = any(y in s.get("pubdate", "") + s.get("epubdate", "") for y in years)
    # First personal author: Crossref may list a consortium first; PubMed writes "Family Initials".
    first = next((a for a in ref["authors"] if a.get("family")), ref["authors"][0])
    first_name = first.get("family") or first.get("name") or ""
    pm_names = [a.get("name", "") for a in (s.get("authors") or [])]
    pm_families = [n.rsplit(" ", 1)[0] if " " in n else n for n in pm_names[:3]]
    pm_first = pm_names[0] if pm_names else ""
    author_ok = any(norm(first_name) == norm(f) for f in pm_families)
    n_ok = len(s.get("authors") or []) == len(ref["authors"])
    if not (title_ok and year_ok and author_ok):
        print(f"  PubMed {ids[0]} does not match {ref['doi']}: title {title_ok}, year {year_ok}, first author "
              f"{author_ok} ({pm_first!r} vs {first_name!r})", file=sys.stderr)
        return ""
    if not n_ok:
        print(f"  note: PubMed lists {len(s.get('authors') or [])} authors for {ref['doi']}, Crossref "
              f"{len(ref['authors'])} (consortium members are counted differently)", file=sys.stderr)
    return ids[0]


def dataset(doi: str, discover_id: int) -> dict:
    a = get("https://api.datacite.org/dois/" + urllib.parse.quote(doi))["data"]["attributes"]
    d = get(f"https://api.pennsieve.io/discover/datasets/{discover_id}")
    if d.get("doi", "").lower() != doi.lower():
        raise RuntimeError(f"Discover {discover_id} has DOI {d.get('doi')}, not {doi}")
    authors = [{"family": c["lastName"], "given": " ".join(x for x in (c.get("firstName"), c.get("middleInitial")) if x),
                "orcid": c.get("orcid") or ""} for c in d.get("contributors", [])]
    creators = [c.get("name") for c in a.get("creators", [])]
    if len(creators) != len(authors):
        print(f"  {doi}: DataCite has {len(creators)} creators, Discover {len(authors)} contributors", file=sys.stderr)
    lic = [r.get("rightsIdentifier", "") for r in a.get("rightsList", [])]
    pub = a.get("publisher")
    return {"doi": doi, "type": "data", "title": plain(a["titles"][0]["title"]), "authors": authors,
            "publisher": pub if isinstance(pub, str) else (pub or {}).get("name", ""), "year": a.get("publicationYear"),
            "license": "CC-BY-4.0" if "cc-by-4.0" in lic else ";".join(lic), "version": str(d.get("version", "")),
            "url": f"https://doi.org/{doi}", "discover": discover_id}


def q(s) -> str:
    return json.dumps(str(s), ensure_ascii=False)


def cff_refs(refs: list[dict]) -> str:
    lines = ["references:"]
    for r in refs:
        lines += [f"  - type: {r['type']}", f"    title: {q(r['title'])}", "    authors:"]
        for a in r["authors"]:
            if a.get("name"):
                lines.append(f"      - name: {q(a['name'])}")
                continue
            lines.append(f"      - family-names: {q(a['family'])}")
            if a.get("given"):
                lines.append(f"        given-names: {q(a['given'])}")
            if a.get("orcid"):
                lines.append(f"        orcid: {q('https://orcid.org/' + a['orcid'])}")
        if r["type"] == "data":
            lines += [f"    publisher:\n      name: {q(r['publisher'])}", f"    year: {r['year']}",
                      f"    version: {q(r['version'])}", f"    license: {r['license']}", f"    url: {q(r['url'])}",
                      f"    doi: {r['doi']}", f"    notes: {q('Demo data; Pennsieve Discover dataset ' + str(r['discover']))}"]
            continue
        lines += [f"    journal: {q(r['journal'])}", f"    year: {r['year']}"]
        for k in ("volume", "issue"):
            if r.get(k):
                lines.append(f"    {k}: {q(r[k])}")
        if r.get("pages"):
            p = str(r["pages"]).split("-", 1)
            lines.append(f"    start: {q(p[0])}")
            if len(p) > 1:
                lines.append(f"    end: {q(p[1])}")
        lines.append(f"    doi: {r['doi']}")
        if r.get("pmid"):
            lines.append(f"    notes: {q('PMID ' + r['pmid'])}")
    return "\n".join(lines) + "\n"


def verify() -> list[dict]:
    refs = []
    for key, doi in PAPERS.items():
        r = crossref(doi)
        r["pmid"] = pubmed(r)
        refs.append(r)
        print(f"{key:14} {r['year']} {len(r['authors']):3} authors  PMID {r['pmid'] or '-':9} {r['journal'][:28]:28} "
              f"{r['volume']}({r['issue']}):{r['pages']}  {r['title'][:60]}")
    for doi, did in DATASETS.items():
        r = dataset(doi, did)
        refs.append(r)
        print(f"{'data':14} {r['year']} {len(r['authors']):3} authors  Discover {did} v{r['version']} {r['license']}  "
              f"{r['title'][:60]}")
    return refs


if __name__ == "__main__":
    refs = verify()
    if "--cff" in sys.argv:
        cff = ROOT / "CITATION.cff"
        head = cff.read_text().split("references:")[0]
        cff.write_text(head + cff_refs(refs))
        print(f"wrote {cff}")

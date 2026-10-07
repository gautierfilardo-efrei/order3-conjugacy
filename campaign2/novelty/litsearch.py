import json, os, re, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET

HDR = {"User-Agent": "python-urllib"}

def _get(url, timeout=40):
    req = urllib.request.Request(url, headers=HDR)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

def crossref(q, rows=8):
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(
        {"query.bibliographic": q, "rows": rows, "select": "DOI,title,author,issued,container-title,type"})
    try:
        d = json.loads(_get(url))
    except Exception as e:
        return [{"error": str(e)[:200]}]
    out = []
    for it in d.get("message", {}).get("items", []):
        out.append({
            "src": "crossref", "doi": it.get("DOI"),
            "title": (it.get("title") or [""])[0],
            "authors": [f"{a.get('family','')}, {a.get('given','')}" for a in it.get("author", [])][:6],
            "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "venue": (it.get("container-title") or [""])[0], "type": it.get("type")})
    return out

def arxiv(q, rows=8):
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode(
        {"search_query": q, "max_results": rows})
    try:
        xmltxt = _get(url)
    except Exception as e:
        return [{"error": str(e)[:200]}]
    ns = {"a": "http://www.w3.org/2005/Atom", "ar": "http://arxiv.org/schemas/atom"}
    out = []
    for e in ET.fromstring(xmltxt).findall("a:entry", ns):
        doi = e.find("ar:doi", ns)
        out.append({
            "src": "arxiv", "id": e.find("a:id", ns).text.strip(),
            "title": re.sub(r"\s+", " ", e.find("a:title", ns).text.strip()),
            "authors": [x.find("a:name", ns).text for x in e.findall("a:author", ns)][:6],
            "year": e.find("a:published", ns).text[:4],
            "doi": doi.text if doi is not None else None,
            "summary": re.sub(r"\s+", " ", e.find("a:summary", ns).text.strip())[:600]})
    return out

def openalex(q, rows=8, key=None):
    key = key or os.environ.get("OPENALEX_API_KEY")
    if not key:
        return [{"error": "no OPENALEX_API_KEY"}]
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(
        {"search": q, "per-page": rows, "api_key": key,
         "select": "id,doi,title,display_name,publication_year,authorships,primary_location,cited_by_count,type"})
    try:
        d = json.loads(_get(url))
    except Exception as e:
        return [{"error": str(e)[:200]}]
    out = []
    for it in d.get("results", []):
        loc = it.get("primary_location") or {}
        out.append({
            "src": "openalex", "oa_id": it.get("id"),
            "doi": (it.get("doi") or "").replace("https://doi.org/", "") or None,
            "title": it.get("display_name"),
            "authors": [a.get("author", {}).get("display_name") for a in it.get("authorships", [])][:6],
            "year": it.get("publication_year"),
            "venue": (loc.get("source") or {}).get("display_name"),
            "cited": it.get("cited_by_count"), "type": it.get("type")})
    return out

def resolve_doi(doi):
    """Verify a DOI via Crossref metadata; returns dict or None."""
    try:
        d = json.loads(_get("https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")))["message"]
    except Exception as e:
        return None
    return {"doi": d.get("DOI"), "title": (d.get("title") or [""])[0],
            "authors": [f"{a.get('family','')}, {a.get('given','')}" for a in d.get("author", [])],
            "year": (d.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "venue": (d.get("container-title") or [""])[0], "volume": d.get("volume"),
            "issue": d.get("issue"), "pages": d.get("page"), "publisher": d.get("publisher"),
            "type": d.get("type")}

def show(hits, n=8):
    for h in hits[:n]:
        if "error" in h:
            print("  ERR", h["error"]); continue
        print(f"  [{h['src']}] {h.get('year')} | {str(h.get('title'))[:95]} | {'; '.join(str(x) for x in h.get('authors', [])[:3])} | {h.get('venue') or h.get('id','')} | doi={h.get('doi')}")

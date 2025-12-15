import json
import random
import re
import time
from urllib.parse import urljoin, urlparse, parse_qs

import requests
from bs4 import BeautifulSoup


BASE = "https://kmt.vander-lingen.nl"
LIST_TEMPLATE = "https://kmt.vander-lingen.nl/data/reaction/doi/{doi}/start/{start}"
DEFAULT_DOI = "10.1021/jacsau.4c01276"

# Common solvent mapping for nicer names
PREFERRED_SOLVENTS = {
    "ClCCl": "dichloromethane",
    "CO[H]": "methanol",
    "C1CCCO1": "tetrahydrofuran",
    "CC(=O)C": "acetone",
    "CC#N": "acetonitrile",
    "CCOCC": "diethyl ether",
    "C1(=CC=CC=C1)C": "toluene",
    "O": "water",
}

name_cache = {}

# Known names for key compounds in this DOI
KNOWN_COMPOUND_NAMES = {
    # Farnesyl acetate and derivatives
    r"C(C)(=O)OC\\C=C(/C)\\CC\\C=C(/C)\\CCC=C(C)C": "(2E,6E)-Farnesyl Acetate",
    r"C(C)(=O)OC\\C=C(\\CC\\C=C(\\CC\\C=C(\\C=O)/C)/C)/C": "(2E,6E,10E)-12-oxo-3,7,11-trimethyldodeca-2,6,10-trien-1-yl acetate",
    r"C(C)(=O)OC\\C=C(\\CC\\C=C(\\CC\\C=C(\\CO)/C)/C)/C": "(2E,6E,10E)-12-Hydroxy-3,7,11-trimethyldodeca-2,6,10-trien-1-yl acetate",
    r"C(C)(=O)OC\\C=C(\\CC\\C=C(\\CCC1OC1(C)CO)/C)/C": "(2E,6E)-9-(3-(Hydroxymethyl)-3-methyloxiran-2-yl)-3,7-dimethylnona-2,6-dien-1-yl acetate",
    r"C(C1=CC=CC=C1)(=O)OC[C@@]1(O[C@H]1CC\\C(=C\\CC\\C(=C\\COC(C)=O)\\C)\\C)C": "((2S,3S)-3-((3E,7E)-9-Acetoxy-3,7-dimethylnona-3,7-dien-1-yl)-2-methyloxiran-2-yl)methyl benzoate",
}

def _norm_smiles(s):
    return s.replace("\\", "").replace("/", "")


def make_session():
    s = requests.Session()
    s.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
            "(KHTML, like Gecko) Chrome/120.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.9",
            "Connection": "keep-alive",
        }
    )
    return s


def extract_reactions_from_list(html):
    rxns = []
    # Inline JS array pattern
    for m in re.finditer(r"reactions\.push\(\s*['\"]([\s\S]*?)['\"]\s*\)", html):
        s = m.group(1)
        if s:
            rxns.append(s)
    # HTML attribute fallback
    for m in re.finditer(r"data-reaction-smiles\s*=\s*['\"]([^'\"]+)['\"]", html):
        val = m.group(1)
        if val:
            rxns.append(val)
    return rxns


def find_next_page(html):
    soup = BeautifulSoup(html, "lxml")
    for a in soup.select("a"):
        t = (a.get_text() or "").strip().lower()
        if t in {"next", ">", "»"}:
            href = a.get("href")
            if href:
                return urljoin(BASE, href)
    return None


def extract_text_blocks(el):
    texts = []
    for s in el.stripped_strings:
        texts.append(s)
    return texts


def parse_details_page(html):
    soup = BeautifulSoup(html, "lxml")
    text = "\n".join(soup.stripped_strings)
    reactant_smiles = []
    solvents = []
    product_smiles = []
    product_names = []

    for table in soup.select("table"):
        rows = table.select("tr")
        for r in rows:
            cells = r.find_all(["th", "td"])[:2]
            if len(cells) < 2:
                continue
            key = (cells[0].get_text(" ", strip=True) or "").lower()
            val = cells[1]
            vals = extract_text_blocks(val)
            if "reactant" in key and "smiles" in key:
                for v in vals:
                    if v and v != "SMILES":
                        reactant_smiles.append(v)
            elif "reactant" in key and "solvent" in key:
                solvents.extend([v for v in vals if v])
            elif "product" in key and "smiles" in key:
                for v in vals:
                    if v and v != "SMILES":
                        product_smiles.append(v)
            elif "product" in key and ("name" in key or "product" == key):
                product_names.extend([v for v in vals if v])

    if not reactant_smiles or not solvents or not product_smiles:
        for dl in soup.select("dl"):
            items = list(dl.children)
            for i in range(0, len(items) - 1, 2):
                dt = items[i]
                dd = items[i + 1]
                if getattr(dt, "name", None) != "dt" or getattr(dd, "name", None) != "dd":
                    continue
                key = (dt.get_text(" ", strip=True) or "").lower()
                vals = extract_text_blocks(dd)
                if "reactant" in key and "smiles" in key:
                    reactant_smiles.extend([v for v in vals if v and v != "SMILES"])
                elif "reactant" in key and "solvent" in key:
                    solvents.extend([v for v in vals if v])
                elif "product" in key and "smiles" in key:
                    product_smiles.extend([v for v in vals if v and v != "SMILES"])
                elif "product" in key and ("name" in key or "product" == key):
                    product_names.extend([v for v in vals if v])

    if not reactant_smiles:
        for m in re.finditer(r"SMILES\s*[:=]\s*([^\s]+)", text, flags=re.I):
            reactant_smiles.append(m.group(1))

    data = {
        "reactant_smiles": sorted(set([s for s in reactant_smiles if s])),
        "solvents": sorted(set([s for s in solvents if s])),
        "product_smiles": sorted(set([s for s in product_smiles if s])),
        "product_name": product_names[0] if product_names else None,
    }
    return data


def parse_reaction_string(s):
    parts = s.split(">")
    while len(parts) < 3:
        parts.append("")
    reactants = [p.strip() for p in parts[0].split(".") if p.strip()]
    solvents = [p.strip() for p in parts[1].split(".") if p.strip()]
    products = [p.strip() for p in parts[2].split(".") if p.strip()]
    return {
        "smiles": s,
        "reactant_smiles": reactants,
        "solvents": solvents,
        "product_smiles": products,
    }


def extract_doi_from_arg(arg):
    try:
        if arg.startswith("http"):
            p = urlparse(arg)
            m = re.search(r"/doi/(.+?)(?:/start|$)", p.path)
            if m:
                return m.group(1)
            return None
        return arg
    except Exception:
        return None


def extract_dois_from_archive(html):
    dois = []
    soup = BeautifulSoup(html, "lxml")
    for a in soup.select("a[href]"):
        href = a.get("href") or ""
        if "/data/reaction/doi/" in href:
            m = re.search(r"/doi/(.+?)(?:/start|$)", href)
            if m:
                dois.append(m.group(1))
    return sorted(set(dois))


def resolve_name_with_pubchem(smiles):
    key = f"name:{smiles}"
    if key in name_cache:
        return name_cache[key]
    try:
        u = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
            + requests.utils.quote(smiles, safe="")
            + "/property/IUPACName/JSON"
        )
        r = requests.get(u, timeout=30)
        if r.status_code == 200:
            j = r.json()
            props = j.get("PropertyTable", {}).get("Properties", [])
            if props:
                nm = props[0].get("IUPACName")
                if nm:
                    name_cache[key] = nm
                    return nm
        # fallback to synonyms
        u2 = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
            + requests.utils.quote(smiles, safe="")
            + "/synonyms/JSON"
        )
        r2 = requests.get(u2, timeout=30)
        if r2.status_code == 200:
            j2 = r2.json()
            info = j2.get("InformationList", {}).get("Information", [])
            if info:
                syns = info[0].get("Synonym", [])
                for s in syns:
                    if not s:
                        continue
                    ls = s.lower()
                    if ls.startswith("cid"):
                        continue
                    if re.fullmatch(r"\d{2,7}-\d{2}-\d", s):
                        continue
                    name_cache[key] = s
                    return s
    except Exception:
        pass
    name_cache[key] = None
    return None


def resolve_name(smiles):
    # direct known mapping first
    ns = _norm_smiles(smiles)
    for k, nm in KNOWN_COMPOUND_NAMES.items():
        if smiles == k or ns == _norm_smiles(k):
            return nm
    nm = resolve_name_with_pubchem(smiles)
    if nm:
        return nm
    # try CACTUS iupac_name
    key = f"cactus:{smiles}"
    if key in name_cache:
        return name_cache[key]
    try:
        u = (
            "https://cactus.nci.nih.gov/chemical/structure/"
            + requests.utils.quote(smiles, safe="")
            + "/iupac_name"
        )
        r = requests.get(u, timeout=30)
        if r.status_code == 200:
            txt = r.text.strip()
            if txt and "Not Found" not in txt:
                name_cache[key] = txt
                return txt
    except Exception:
        pass
    name_cache[key] = None
    return None


def pick_primary_solvent(solvent_smiles):
    for s in solvent_smiles:
        if s in PREFERRED_SOLVENTS:
            return s, PREFERRED_SOLVENTS[s]
    if solvent_smiles:
        s0 = solvent_smiles[0]
        nm = resolve_name_with_pubchem(s0)
        return s0, nm
    return None, None


def scrape_all(max_pages=15, doi=None):
    s = make_session()
    start_url = LIST_TEMPLATE.format(doi=(doi or DEFAULT_DOI), start=0)
    url = start_url
    results = []
    seen_pages = set()
    pages = 0
    while url and url not in seen_pages and pages < max_pages:
        seen_pages.add(url)
        r = s.get(url, timeout=30)
        if r.status_code != 200:
            break
        html = r.text
        rxn_strings = extract_reactions_from_list(html)
        for rs in rxn_strings:
            item = parse_reaction_string(rs)
            item["page_url"] = url
            rsmi = item["reactant_smiles"]
            psmi = item["product_smiles"]
            s_smiles, s_name = pick_primary_solvent(item["solvents"])
            results.append(
                {
                    "solvent": s_name,
                    "reactant_smiles": rsmi,
                    "solvent_smiles": [s_smiles] if s_smiles else [],
                    "product_smiles": psmi,
                }
            )
        next_url = find_next_page(html)
        url = next_url
        pages += 1
        time.sleep(random.uniform(0.6, 1.5))
    return results


def main():
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--debug-list":
        s = make_session()
        doi = DEFAULT_DOI
        if len(sys.argv) > 3 and sys.argv[2] == "--doi":
            doi = sys.argv[3]
        r = s.get(LIST_TEMPLATE.format(doi=doi, start=0), timeout=30)
        soup = BeautifulSoup(r.text, "lxml")
        html = r.text
        m = re.search(r"id=\"title-0\"", html)
        if m:
            i = m.start()
            print("HTML AROUND title-0:")
            print(html[max(0, i-300): i+300])
        m2 = re.search(r"id=\"display-0\"", html)
        if m2:
            i = m2.start()
            print("HTML AROUND display-0:")
            print(html[max(0, i-400): i+400])
        for m in re.finditer(r"data-reaction-smiles=\"(.*?)\"", html):
            val = m.group(1)
            if val:
                print("REACTION_SMILES", val[:120])
        details_links = [a for a in soup.select("a") if a.get_text(strip=True).lower() == "details"]
        print("DETAILS COUNT", len(details_links))
        if details_links:
            a = details_links[0]
            print("DETAILS ATTRS", json.dumps({k: a.get(k) for k in a.attrs}, ensure_ascii=False))
            parent = a.find_parent()
            for _ in range(5):
                if parent is None:
                    break
                txt = parent.get_text(" ", strip=True)
                print("PARENT", parent.name, len(txt))
                parent = parent.find_parent()
            print("CARD TEXT SAMPLE")
            parent = details_links[0].find_parent("div")
            if parent:
                print(parent.get_text(" \n", strip=True)[:1000])
        for a in soup.select("a")[:50]:
            t = a.get_text(strip=True)
            if t.lower() == "details":
                attrs = {k: a.get(k) for k in a.attrs}
                print("DETAILS", json.dumps(attrs, ensure_ascii=False))
            else:
                print(t, "->", a.get("href"))
        print("SCRIPTS")
        scripts = soup.select("script")
        for sc in scripts:
            src = sc.get("src")
            if src:
                u = urljoin(BASE, src)
                print("SCRIPT", u)
                try:
                    sr = s.get(u, timeout=30)
                    body = sr.text
                    if "title-" in body or "data-reaction" in body or "window.open" in body:
                        print("JS HINT", u)
                        idx = body.find("title-")
                        if idx != -1:
                            print(body[max(0, idx-200): idx+200])
                    for m in re.finditer(r"fetch\(\s*['\"](.*?)['\"]", body):
                        print("FETCH", m.group(1))
                    for m in re.finditer(r"axios\.(get|post)\(\s*['\"](.*?)['\"]", body):
                        print("AXIOS", m.group(2))
                    for m in re.finditer(r"\$\.ajax\(\s*\{[\s\S]*?url\s*:\s*['\"](.*?)['\"]", body):
                        print("AJAX", m.group(1))
                    for m in re.finditer(r"\$\.get\(\s*['\"](.*?)['\"]", body):
                        print("JQGET", m.group(1))
                    if "data-reaction-smiles" in body:
                        print("HAS data-reaction-smiles IN JS")
                    if "dataModal" in body:
                        print("HAS dataModal IN JS")
                except Exception:
                    pass
            else:
                body = sc.get_text() or ""
                if body.strip():
                    print("INLINE SCRIPT")
                    print(body[:1000])
        return
    import sys
    max_pages = 15
    args = sys.argv[1:]
    targets = []
    combined_out = None
    i = 0
    archive_limit = None
    while i < len(args):
        a = args[i]
        if a == "--max-pages" and i + 1 < len(args):
            try:
                max_pages = int(args[i + 1])
            except Exception:
                pass
            i += 2
        elif a == "--doi" and i + 1 < len(args):
            targets.append(args[i + 1])
            i += 2
        elif a == "--combined-out" and i + 1 < len(args):
            combined_out = args[i + 1]
            i += 2
        elif a == "--archive" and i + 1 < len(args):
            arch_url = args[i + 1]
            s = make_session()
            try:
                r = s.get(arch_url, timeout=30)
                if r.status_code == 200:
                    ds = extract_dois_from_archive(r.text)
                    if archive_limit is not None:
                        ds = ds[:archive_limit]
                    targets.extend(ds)
            except Exception:
                pass
            i += 2
        elif a == "--archive-limit" and i + 1 < len(args):
            try:
                archive_limit = int(args[i + 1])
            except Exception:
                archive_limit = None
            i += 2
        else:
            if "archive" in a:
                s = make_session()
                try:
                    r = s.get(a, timeout=30)
                    if r.status_code == 200:
                        ds = extract_dois_from_archive(r.text)
                        if archive_limit is not None:
                            ds = ds[:archive_limit]
                        targets.extend(ds)
                except Exception:
                    pass
            else:
                targets.append(a)
            i += 1
    if not targets:
        targets = [DEFAULT_DOI]
    if combined_out is None:
        combined_out = "kmt_reactions_combined.json"
    all_results = []
    for t in targets:
        d = extract_doi_from_arg(t)
        if not d:
            continue
        data = scrape_all(max_pages=max_pages, doi=d)
        all_results.extend(data)
        print(d)
        print(len(data))
    if all_results:
        with open(combined_out, "w", encoding="utf-8") as f:
            json.dump(all_results, f, ensure_ascii=False, indent=2)
        print(combined_out)
        print(len(all_results))


if __name__ == "__main__":
    main()


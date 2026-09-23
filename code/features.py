#!/usr/bin/env python3
"""Feature extraction for RDCB. Sources: Ensembl VEP REST (VCF-format batch),
gnomAD v4 GraphQL (AF + gene constraint), local ClinVar subset (same-codon
evidence from OTHER variants only - G5). Disk-cached, resume-safe."""
import json, os, re, time, csv, sys, urllib.request, urllib.error

BASE = "/home/sandbox/rare-disease-variant-pipeline"
CACHE = os.path.join(BASE, "data/cache")
os.makedirs(os.path.join(CACHE, "vep"), exist_ok=True)
os.makedirs(os.path.join(CACHE, "gnomad"), exist_ok=True)
os.makedirs(os.path.join(CACHE, "constraint"), exist_ok=True)

def http_post(url, payload, tries=5):
    data = json.dumps(payload).encode()
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data,
                headers={"Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=240) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 400:
                raise
            time.sleep(2 ** i + 1)
        except Exception:
            time.sleep(2 ** i + 1)
    raise RuntimeError("failed after retries: " + url)

def vid(c, p, r, a):
    return f"{c}-{p}-{r}-{a}"

# ---------------- VEP ----------------
def to_region(c, p, r, a):
    """VCF-anchored -> VEP region format. Strips common prefix."""
    i = 0
    while i < len(r) and i < len(a) and r[i] == a[i]:
        i += 1
    r2, a2 = r[i:], a[i:]
    start = p + i
    if not r2:
        return f"{c} {start + 1} {start} -/{a2} +"
    end = start + len(r2) - 1
    return f"{c} {start} {end} {r2}/{a2 or '-'} +"


def vep_batch(variants):
    """variants: list of (chrom,pos,ref,alt). Returns dict vid -> feature dict."""
    out, todo = {}, []
    for c, p, r, a in variants:
        v = vid(c, p, r, a)
        fp = os.path.join(CACHE, "vep", v + ".json")
        if os.path.exists(fp):
            out[v] = json.load(open(fp))
        else:
            todo.append((c, p, r, a))
    for i in range(0, len(todo), 50):
        chunk = todo[i:i + 50]
        lines = [to_region(c, int(p), r, a) for c, p, r, a in chunk]
        res = http_post("https://rest.ensembl.org/vep/human/region", {"variants": lines})
        region2vid = {to_region(c, int(p), r, a): vid(c, p, r, a) for c, p, r, a in chunk}
        by_input = {}
        for entry in res:
            k = region2vid.get(entry.get("input", "").strip())
            if not k:
                continue
            tcs = entry.get("transcript_consequences", [])
            tc = (next((t for t in tcs if t.get("mane_select")), None)
                  or next((t for t in tcs if t.get("canonical") == 1), None)
                  or (tcs[0] if tcs else {}))
            feat = {
                "most_severe": entry.get("most_severe_consequence", ""),
                "consequence": tc.get("consequence_terms", []),
                "gene": tc.get("gene_symbol", ""),
                "transcript": tc.get("transcript_id", ""),
                "sift": tc.get("sift_prediction", ""),
                "polyphen": tc.get("polyphen_prediction", ""),
                "protein_start": tc.get("protein_start"),
                "amino_acids": tc.get("amino_acids", ""),
                "codons": tc.get("codons", ""),
                "hgvsc": tc.get("hgvsc", ""),
                "hgvsp": tc.get("hgvsp", ""),
                "impact": tc.get("impact", ""),
            }
            by_input[k] = feat
        for c, p, r, a in chunk:
            v = vid(c, p, r, a)
            feat = by_input.get(v, {"most_severe": "", "consequence": [], "gene": "",
                                    "sift": "", "polyphen": "", "protein_start": None,
                                    "amino_acids": "", "codons": "", "hgvsc": "",
                                    "hgvsp": "", "impact": "", "transcript": ""})
            json.dump(feat, open(os.path.join(CACHE, "vep", v + ".json"), "w"))
            out[v] = feat
        time.sleep(0.34)
    return out

# ---------------- gnomAD AF ----------------
def gnomad_batch(variants):
    """variants: list of (chrom,pos,ref,alt). Returns dict vid -> {af, ac, an}."""
    out, todo = {}, []
    for c, p, r, a in variants:
        v = vid(c, p, r, a)
        fp = os.path.join(CACHE, "gnomad", v + ".json")
        if os.path.exists(fp):
            out[v] = json.load(open(fp))
        else:
            todo.append(v)
    for i in range(0, len(todo), 40):
        chunk = todo[i:i + 40]
        q = "query { " + " ".join(
            f'v{j}: variant(variantId: "{v}", dataset: gnomad_r4) '
            "{ variant_id genome { af ac an } exome { af ac an } }"
            for j, v in enumerate(chunk)) + " }"
        res = http_post("https://gnomad.broadinstitute.org/api", {"query": q})
        data = res.get("data", {})
        for j, v in enumerate(chunk):
            d = data.get(f"v{j}")
            if d is None:
                feat = {"af": 0.0, "ac": 0, "an": 0, "found": False}
            else:
                g = d.get("genome") or {}
                e = d.get("exome") or {}
                af = max(g.get("af") or 0.0, e.get("af") or 0.0)
                ac = (g.get("ac") or 0) + (e.get("ac") or 0)
                an = max(g.get("an") or 0, e.get("an") or 0)
                feat = {"af": af, "ac": ac, "an": an, "found": True}
            json.dump(feat, open(os.path.join(CACHE, "gnomad", v + ".json"), "w"))
            out[v] = feat
        time.sleep(0.5)
    return out

# ---------------- gene constraint ----------------
def constraint_batch(genes):
    out, todo = {}, []
    for g in genes:
        fp = os.path.join(CACHE, "constraint", g + ".json")
        if os.path.exists(fp):
            out[g] = json.load(open(fp))
        else:
            todo.append(g)
    for i in range(0, len(todo), 20):
        chunk = todo[i:i + 20]
        q = "query { " + " ".join(
            f'g{j}: gene(gene_symbol: "{g}", reference_genome: GRCh38) '
            "{ gnomad_constraint { oe_lof oe_lof_upper mis_z lof_z } }"
            for j, g in enumerate(chunk)) + " }"
        res = http_post("https://gnomad.broadinstitute.org/api", {"query": q})
        data = res.get("data", {})
        for j, g in enumerate(chunk):
            d = (data.get(f"g{j}") or {}).get("gnomad_constraint") or {}
            feat = {"oe_lof_upper": d.get("oe_lof_upper"), "mis_z": d.get("mis_z"),
                    "lof_z": d.get("lof_z"), "oe_lof": d.get("oe_lof")}
            json.dump(feat, open(os.path.join(CACHE, "constraint", g + ".json"), "w"))
            out[g] = feat
        time.sleep(0.5)
    return out

# ---------------- ClinVar same-codon (OTHER variants only) ----------------
AA = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
      "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
      "Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*"}
PROT_RE = re.compile(r"p\.(?:\()?([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|=)?")

def build_codon_index():
    """gene -> prot_pos -> list of (to_aa1, sig) for P/LP missense in ClinVar subset."""
    idx = {}
    with open(os.path.join(BASE, "data/clinvar_grch38_small.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row["ClinSigSimple"] != "1" and row["ClinicalSignificance"] not in (
                "Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"):
                continue
            if row["ReviewStatus"] == "no assertion criteria provided":
                continue
            m = PROT_RE.search(row["Name"])
            if not m or not m.group(3) or m.group(3) == "=":
                continue
            f3, pos, t3 = m.group(1), int(m.group(2)), m.group(3)
            if f3 not in AA or t3 not in AA:
                continue
            idx.setdefault(row["GeneSymbol"], {}).setdefault(pos, []).append(
                (AA[f3] + str(pos) + AA[t3], row["ClinicalSignificance"],
                 row["Chromosome"] + "-" + row["PositionVCF"] + "-" +
                 row["ReferenceAlleleVCF"] + "-" + row["AlternateAlleleVCF"]))
    return idx

def same_codon(codon_idx, gene, protein_start, amino_acids, exclude_vid):
    """PS1: identical AA change known P/LP. PM5: different AA change at same codon
    known P/LP. Excludes the evaluated variant itself (G5)."""
    if not gene or not protein_start or not amino_acids:
        return False, False
    try:
        pos = int(protein_start)
    except (TypeError, ValueError):
        return False, False
    aas = str(amino_acids)
    if len(aas) != 3:
        return False, False
    my_change = aas[0] + str(pos) + aas[2]
    ps1 = pm5 = False
    for change, sig, other_vid in codon_idx.get(gene, {}).get(pos, []):
        if other_vid == exclude_vid:
            continue  # G5: the evaluated variant's own ClinVar record is never evidence
        if change == my_change:
            ps1 = True
        else:
            pm5 = True
    return ps1, pm5

if __name__ == "__main__":
    # smoke test on 3 variants
    tv = [("13", "32339662", "T", "A"), ("17", "7674220", "C", "T"), ("X", "154363127", "G", "A")]
    v = vep_batch(tv)
    g = gnomad_batch(tv)
    c = constraint_batch(["BRCA2", "TP53", "F8"])
    ci = build_codon_index()
    for (cc, p, r, a) in tv:
        k = vid(cc, p, r, a)
        f = v[k]
        ps1, pm5 = same_codon(ci, f["gene"], f["protein_start"], f["amino_acids"], k)
        print(k, f["most_severe"], f["gene"], f["sift"], f["polyphen"],
              "AF=", g[k]["af"], "PS1=", ps1, "PM5=", pm5)
    print("constraint BRCA2:", c.get("BRCA2"))
    print("codon index genes:", len(ci))

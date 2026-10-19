#!/usr/bin/env python3
"""Shared pieces for the Stage 3 prototypes: FASTA I/O, reference loaders, BLAST and best-hit calls.

Each assignment arm is the same operation on a different marker and reference:
  VT   SSU flank vs MaarjAM          unit = VTX code
  SH   full ITS  vs UNITE            unit = SH code
  LSU  LROR-FLR2 vs Delavaux (V18)   unit = reference species (nearest reference, Phase 1)

Loaders return {ref_id: info} and a FASTA copy of the reference with simple unique IDs
(r0000001...), so BLAST never parses pipes or spaces in headers. Every info dict has
  label  the unit assigned
  acc    INSDC accession of the reference sequence (used to drop self-hits)
  taxon  readable name for reports
  genus  genus in the reference's own taxonomy ('' if none)

assign() is vt_check.assign() generalised: best hit by bit score; a call needs identity,
coverage (of the shorter sequence) and e-value to pass. Coverage leaves out Ns in the reference:
Delavaux V18 pads many records with terminal Ns (and the gDAT MaarjAM set a few), which would
otherwise make exact matches to those records fail the coverage test; tied top scores across units give
'ambiguous'; margin = best identity minus the best identity to any other unit. Hits to the
accessions in drop_accs (normally every copy of the culture reference) are dropped first: once a
database has ingested the Stefani deposits, a copy would otherwise take its label from itself or
from its own culture-mates.
"""
import collections
import csv
import gzip
import hashlib
import os
import re
import subprocess

# genus -> order (Glomeromycota), from crosslineage.py; extend as needed. Changes from the Stage 1 dict:
# Viscospora and the other Glomerales genera of da Silva et al. 2024 (Taxonomy 4:41; Septoglomeraceae,
# Sclerocystaceae, Dominikiaceae, revised Glomeraceae) plus Microviscospora and Melanoglomus go to
# Glomerales; Albahypha and Alborhynchus go to Entrophosporales, as V18 files them (Entrophosporaceae).
ORDER = {
    **dict.fromkeys(["Rhizophagus", "Glomus", "Funneliformis", "Septoglomus", "Sclerocystis",
                     "Dominikia", "Microdominikia", "Kamienskia", "Nanoglomus", "Oehlia",
                     "Rhizoglomus", "Simiglomus", "Halonatospora", "Microkamienskia",
                     "Epigeocarpum", "Orientoglomus", "Viscospora", "Microviscospora", "Melanoglomus",
                     "Blaszkowskia", "Funneliglomus", "Complexispora", "Sclerocarpum", "Parvocarpum",
                     "Silvaspora", "Macrodominikia"], "Glomerales"),
    **dict.fromkeys(["Claroideoglomus", "Entrophospora", "Albahypha", "Alborhynchus"], "Entrophosporales"),
    **dict.fromkeys(["Diversispora", "Redeckera", "Otospora", "Acaulospora", "Kuklospora",
                     "Pacispora", "Gigaspora", "Scutellospora", "Racocetra", "Cetraspora",
                     "Dentiscutata", "Fuscutata", "Quatunica", "Intraornatospora",
                     "Paradentiscutata", "Bulbospora", "Corymbiglomus", "Sieverdingia"],
                    "Diversisporales"),
    **dict.fromkeys(["Archaeospora", "Ambispora", "Geosiphon", "Palaeospora"], "Archaeosporales"),
    **dict.fromkeys(["Paraglomus", "Innospora"], "Paraglomerales"),
}
# one genus under two names; compared as equal (and listed) where genera are checked
GENUS_SYNONYMS = {"Rhizoglomus": "Rhizophagus", "Claroideoglomus": "Entrophospora"}
UNNAMED = re.compile(r"\bsp\.?(\s|$)")  # 'Septoglomus sp.', 'Paraglomus sp. BR105' are not named species
ACC = re.compile(r"^([A-Z]{1,2}_?\d{5,8}(?:\.\d+)?)")


def genus_of(organism):
    return organism.split()[0].strip("[]") if organism else ""


def is_named(organism):
    return bool(organism) and not UNNAMED.search(organism)


def base_acc(acc):
    return acc.split(".")[0] if acc else ""


def open_text(path):
    """Plain or gzipped text."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    """{first header token: sequence (upper case)}; tolerates CRLF, wrapped lines and .gz."""
    seqs, name, chunks = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_fasta(path, seqs):
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")


def sha256(path, n=None):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()[:n] if n else h.hexdigest()


# ---------------------------------------------------------------- reference loaders

def _records(path):
    header, seq = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            elif line:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def parse_maarjam(header):
    """gDAT snapshot: 'gi|MOA01016|gb|AB555676| Glomeraceae Glomus sp. VTX00223'
    or the rewritten 2016 form 'VTX00223__AB555676__Glomeraceae_Glomus_sp.'"""
    m = re.match(r"^(VTX\d+)__([^_]+)__(\S+)", header)
    if m:
        toks = m.group(3).split("_")
        return {"label": m.group(1), "acc": m.group(2), "family": toks[0],
                "taxon": " ".join(toks[1:]), "genus": toks[1] if len(toks) > 1 else ""}
    idpart, _, rest = header.partition(" ")
    p, toks = idpart.split("|"), rest.split()
    if len(p) < 4 or not toks or not re.fullmatch(r"VTX\d+", toks[-1]):
        raise ValueError(f"unexpected MaarjAM header: {header}")
    return {"label": toks[-1], "acc": p[3], "maarjam_id": p[1], "family": toks[0],
            "taxon": " ".join(toks[1:-1]), "genus": toks[1] if len(toks) > 2 else ""}


def parse_unite(header):
    """UNITE general FASTA release, 'Name|ACCESSION|SHnnnnnnn.vvFU|reps|k__Fungi;p__...;s__Name'.
    Anything carrying an SH code works: the code, accession and taxonomy are found by pattern."""
    sh = re.search(r"SH\d+\.\d+FU", header)
    if not sh:
        raise ValueError(f"no SH code in UNITE header: {header}")
    fields = re.split(r"[|\s]", header)
    acc = next((f for f in fields if (ACC.match(f) or f.startswith("UDB")) and not re.fullmatch(r"SH\d+\.\d+FU", f)), "")
    tax = dict(re.findall(r"([kpcofgs])__([^;|]+)", header))
    species = tax.get("s", fields[0]).replace("_", " ")
    genus = tax.get("g", "")
    if genus.lower().startswith("unidentified") or "incertae" in genus.lower():  # e.g. Glomeraceae_gen_Incertae_sedis
        genus = ""
    return {"label": sh.group(0), "acc": acc, "taxon": species, "genus": genus,
            "phylum": tax.get("p", ""), "reftype": next((f for f in fields if f.startswith(("reps", "refs"))), "")}


def parse_delavaux(header):
    """Delavaux LSU DB (V18): 'MK348930.1_Septoglomeraceae_Funneliglomus_sanmartinense',
    'NG_060196.1_Gigasporaceae_Racocetra_castanea', outgroups 'Z19136.1_Schizosaccharomyces_pombe_Asco'.
    Ten V18 records put the database's genus clade before a binomial from another genus
    ('Septoglomus_Viscospora_viscosa', 'Entrophospora_Albahypha_drummondi'); the pipeline's clade
    lists use the first genus, so genus = first token and the label keeps the whole name."""
    m = re.match(r"^([A-Z]{1,2}_?\d+(?:\.\d+)?)_(.+)$", header.split()[0])
    if not m:
        raise ValueError(f"unexpected Delavaux header: {header}")
    toks = m.group(2).split("_")
    fam = toks[0] if toks[0].endswith("aceae") else ""
    rest = toks[1:] if fam else toks
    if not rest:
        raise ValueError(f"no genus in Delavaux header: {header}")
    return {"label": " ".join(rest), "acc": m.group(1), "family": fam, "genus": rest[0],
            "taxon": " ".join(rest), "amf": bool(fam)}


PARSERS = {"maarjam": parse_maarjam, "unite": parse_unite, "delavaux": parse_delavaux}


def load_reference(path, kind, keep=None):
    """Return ({ref_id: info}, fasta text with ref_ids, records read). Duplicate headers are fine.
    keep: regex; records whose header does not match are skipped before parsing."""
    parse = PARSERS[kind]
    pat = re.compile(keep) if keep else None
    refs, out, n = {}, [], 0
    for header, seq in _records(path):
        n += 1
        if pat and not pat.search(header):
            continue
        rid = f"r{len(refs) + 1:07d}"
        info = parse(header)
        info["header"] = header
        info["n_N"] = seq.upper().count("N")
        refs[rid] = info
        out.append(f">{rid}\n{seq}\n")
    if not refs:
        raise SystemExit(f"no sequences in {path}" + (f" match {keep!r}" if keep else ""))
    return refs, "".join(out), n


# ---------------------------------------------------------------- BLAST and calls

def run_blast(query, ref_fasta, workdir, name, task="megablast", threads=4, evalue=1e-10,
              max_targets=5000):
    """blastn every query against the reference; one HSP per subject. Returns {query: [hit]}."""
    os.makedirs(workdir, exist_ok=True)
    db = os.path.join(workdir, name)
    with open(db + ".fasta", "w") as fh:
        fh.write(ref_fasta)
    subprocess.run(["makeblastdb", "-in", db + ".fasta", "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL)
    cols = "qseqid sseqid pident length qlen slen evalue bitscore"
    res = subprocess.run(
        ["blastn", "-task", task, "-query", query, "-db", db, "-dust", "no",
         "-outfmt", f"6 {cols}", "-evalue", str(evalue), "-max_target_seqs", str(max_targets),
         "-max_hsps", "1", "-num_threads", str(threads)],
        check=True, capture_output=True, text=True)
    hits = collections.defaultdict(list)
    for line in res.stdout.splitlines():
        q, s, pid, ln, ql, sl, ev, bs = line.split("\t")
        hits[q].append({"ref": s, "pident": float(pid), "length": int(ln), "qlen": int(ql),
                        "slen": int(sl), "evalue": float(ev), "bitscore": float(bs)})
    return hits


EMPTY_CALL = {"call": "no_hit", "label": "", "pident": "", "cov": "", "tied": "", "best_label": "",
              "best_pident": "", "best_taxon": "", "best_genus": "", "margin": "", "dropped": 0}


def assign(hs, refs, min_id, min_cov, max_ev, drop_accs=None):
    """Best-hit call for one query (vt_check.assign, after dropping hits to drop_accs)."""
    row = dict(EMPTY_CALL)
    if drop_accs:
        kept = [h for h in hs if base_acc(refs[h["ref"]]["acc"]) not in drop_accs]
        row["dropped"] = len(hs) - len(kept)
        hs = kept
    if not hs:
        return row
    for h in hs:
        h["cov"] = h["length"] / min(h["qlen"], max(1, h["slen"] - refs[h["ref"]].get("n_N", 0)))
        h["label"] = refs[h["ref"]]["label"]
    best = max(hs, key=lambda h: (h["bitscore"], h["pident"]))
    info = refs[best["ref"]]
    best_by = {}
    for h in hs:
        if h["pident"] > best_by.get(h["label"], -1):
            best_by[h["label"]] = h["pident"]
    others = [p for u, p in best_by.items() if u != best["label"]]
    row.update(best_label=best["label"], best_pident=round(best["pident"], 2), best_taxon=info["taxon"],
               best_genus=info.get("genus", ""), cov=round(best["cov"], 3),
               margin=round(best["pident"] - max(others), 2) if others else "")
    passing = [h for h in hs if h["pident"] >= min_id and h["cov"] >= min_cov and h["evalue"] <= max_ev]
    if not passing:
        row["call"] = "unassigned"
        return row
    top = max(h["bitscore"] for h in passing)
    tied = sorted({h["label"] for h in passing if h["bitscore"] == top})
    pbest = max((h for h in passing if h["bitscore"] == top), key=lambda h: h["pident"])
    row.update(pident=round(pbest["pident"], 2), label=tied[0] if len(tied) == 1 else "",
               call="assigned" if len(tied) == 1 else "ambiguous",
               tied=",".join(tied) if len(tied) > 1 else "")
    return row

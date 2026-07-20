#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TRIEKY - Tandem Repeat Inspector for Emergent K-mer varietY

Post-processing tool for PacBio TRGT VCF files.

Author: Zoé Dmitrievsky
Neurogenetics Laboratory, Hospices Civils de Lyon

Main ideas
-----------
- Use a small TSV config to describe, per TRID/gene:
    * which motif sizes (k) to analyse,
    * an optional first-letter filter (e.g. "A"),
    * a minimal count threshold,
    * whether to keep only unexpected motifs (vs. expected ones in MOTIFS=),
    * whether to group rotationally equivalent motifs.

- For each TRGT VCF record:
    * split ALT into individual alleles (ALT1, ALT2, ...) on comma,
    * for each allele, extract a cleaned sequence (A/T/C/G only),
    * get expected motifs from the MOTIFS= field,
    * count k-mers according to the config,
    * apply filters and optional rotational grouping,
    * write a human-readable report, one block per allele.
"""

import sys
import csv
import re
from collections import Counter
from typing import Dict, List, Tuple

# Only accept canonical bases in cleaned ALT sequences
VALID = set("ATCG")


# ---------------------------------------------------------------------------
# Basic utility helpers
# ---------------------------------------------------------------------------

def clean_seq(s: str) -> str:
    """
    Uppercase the sequence and remove any character that is not A/T/C/G.
    This is applied to ALT subsequences (per allele) which may contain extra
    symbols or punctuation from TRGT.
    """
    return re.sub(r'[^ATCG]', '', s.upper())


def parse_info(info: str) -> Dict[str, str]:
    """
    Parse the INFO field of a VCF line into a dict.
    Example: "TRID=DMPK;MOTIFS=CAG,CAA" -> {"TRID":"DMPK", "MOTIFS":"CAG,CAA"}
    """
    out: Dict[str, str] = {}
    for item in info.split(';'):
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            out[k] = v
        else:
            # Standalone flag without value
            out[item] = ""
    return out


def rotations(s: str) -> List[str]:
    """
    Return all circular rotations of a string.
    E.g. "CAG" -> ["CAG", "AGC", "GCA"].
    """
    return [s[i:] + s[:i] for i in range(len(s))] if s else [s]


def is_rotation(a: str, b: str) -> bool:
    """
    Check if b is a circular rotation of a.
    """
    return len(a) == len(b) and (b in (a + a))


def expected_sets(expected: List[str]) -> List[set]:
    """
    For each expected motif, build the set of all its rotations.

    Example:
        expected = ["CAG"]
        -> [ {"CAG", "AGC", "GCA"} ]
    """
    sets: List[set] = []
    for m in expected:
        if set(m) <= VALID:
            sets.append(set(rotations(m)))
    return sets


def is_expected(m: str, exp_sets: List[set]) -> bool:
    """
    Return True if motif m belongs to any of the rotational equivalence sets
    derived from expected motifs.
    """
    return any(m in s for s in exp_sets)


# ---------------------------------------------------------------------------
# K-mer counting helpers
# ---------------------------------------------------------------------------

def count_sliding(seq: str, k: int) -> Counter:
    """
    Count k-mers in a sliding window fashion (step = 1 base).

    Example:
        seq = "AAAAAG"
        k = 5
        k-mers: "AAAAA", "AAAAG"
    """
    cnt: Counter = Counter()
    for i in range(0, len(seq) - k + 1):
        m = seq[i:i + k]
        if len(m) == k and set(m) <= VALID:
            cnt[m] += 1
    return cnt


def count_frame(seq: str, k: int, phase: int) -> Counter:
    """
    Count k-mers in a fixed reading frame:
    - Start at 'phase' (0..k-1),
    - Step by k each time (non-overlapping).

    Example:
        seq = "CAGCAGCAG"
        k = 3, phase = 0 -> CAG, CAG, CAG
    """
    cnt: Counter = Counter()
    for i in range(phase, len(seq) - k + 1, k):
        m = seq[i:i + k]
        if len(m) == k and set(m) <= VALID:
            cnt[m] += 1
    return cnt


# ---------------------------------------------------------------------------
# Rotational grouping
# ---------------------------------------------------------------------------

def group_by_rotations(counts: Dict[str, int]) -> List[Tuple[str, Dict[str, int]]]:
    """
    Group motifs that are rotationally equivalent.

    Input:
        counts: dict motif -> count

    Output:
        List of (rep, member_counts) where:
          - rep: representative motif (the most frequent in the group)
          - member_counts: dict motif -> count for all motifs in the group

    Groups are sorted by descending total count.
    """
    remaining = set(counts.keys())
    groups: List[Tuple[str, Dict[str, int]]] = []

    while remaining:
        seed = next(iter(remaining))
        eq = set(rotations(seed))
        members = sorted(
            [m for m in counts if m in eq],
            key=lambda x: (-counts[x], x)
        )
        member_counts = {m: counts[m] for m in members}
        rep = members[0]
        groups.append((rep, member_counts))
        remaining -= set(members)

    # sort groups by total count (descending)
    groups.sort(key=lambda g: -sum(g[1].values()))
    return groups


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def load_config(tsv_path: str):
    """
    Load a simple TSV config describing how to post-process each TRID.

    Required columns:
      - TRID
      - motif_sizes        (comma-separated list of k, may be empty)
      - startswith         ("ANY" or a single base)
      - count_min          (integer threshold)
      - keep_only_unexpected (true/false)
      - group_rotations    (true/false)

    Returns:
      cfg: dict TRID -> config dict
    """
    cfg: Dict[str, Dict] = {}
    with open(tsv_path, newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        required = {
            "TRID",
            "motif_sizes",
            "startswith",
            "count_min",
            "keep_only_unexpected",
            "group_rotations",
        }
        missing = required - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"Missing config columns: {sorted(missing)}")

        for row in r:
            trid = row["TRID"].strip()
            if not trid:
                continue

            # Parse motif sizes
            if row["motif_sizes"].strip():
                sizes = [int(x) for x in row["motif_sizes"].split(',') if x.strip()]
            else:
                sizes = []

            # First-letter constraint
            starts = row["startswith"].strip().upper() or "ANY"
            if starts == "ANY":
                starts = "ANY"

            cfg[trid] = {
                "sizes": sizes,
                "startswith": starts,
                "count_min": int(row["count_min"]),
                "keep_only_unexpected": row["keep_only_unexpected"].strip().lower() == "true",
                "group_rotations": row["group_rotations"].strip().lower() == "true",
            }
    return cfg


# ---------------------------------------------------------------------------
# TRGT VCF record iterator (per allele)
# ---------------------------------------------------------------------------

def iter_trgt_records(vcf_path: str):
    """
    Iterate over TRGT records in a VCF file, splitting ALT by allele.

    For each non-header line:
      - Parse INFO, retrieve TRID and MOTIFS
      - Split ALT on comma into ALT1, ALT2, ...
      - For each allele:
          * clean the allele sequence to A/T/C/G only
          * yield a per-allele record

    Yields (one per allele):
      {
        "chrom": chrom,
        "pos": pos,
        "trid": trid,
        "allele_index": 1-based index (1,2,...),
        "alt_seq": cleaned ALT sequence for that allele,
        "motifs_from_vcf": list of expected motifs (if any)
      }
    """
    with open(vcf_path, 'r') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue

            chrom, pos, _id, ref, alt, qual, flt, info = parts[:8]
            info_map = parse_info(info)
            trid = info_map.get("TRID")
            if not trid:
                # Not a TRGT locus or TRID missing
                continue

            motifs_from_vcf: List[str] = []
            if "MOTIFS" in info_map and info_map["MOTIFS"]:
                motifs_from_vcf = [
                    m.strip().upper()
                    for m in info_map["MOTIFS"].split(',')
                    if m.strip()
                ]

            # ALT may contain multiple alleles, separated by comma.
            # Example: "SEQALLELE1,SEQALLELE2"
            alt_alleles_raw = alt.split(',')

            for allele_idx, allele_raw in enumerate(alt_alleles_raw, start=1):
                alt_seq = clean_seq(allele_raw)
                yield {
                    "chrom": chrom,
                    "pos": pos,
                    "trid": trid,
                    "allele_index": allele_idx,
                    "alt_seq": alt_seq,
                    "motifs_from_vcf": motifs_from_vcf,
                }


# ---------------------------------------------------------------------------
# Junction filter: remove artificial k-mers from perfect repeats of expected motifs
# ---------------------------------------------------------------------------

def build_banned_kmers_from_expected(expected_list: List[str], k: int) -> set:
    """
    Build a set of k-mers to ignore because they are trivial junctions of
    a shorter expected motif repeated perfectly.

    For each expected motif m of length L < k:
      - build a repeated string rep = m * r long enough
      - ban all substrings of length k from rep

    Example:
      expected = ["AAAGG"], L=5, k=6
      rep = "AAAGGAAAGG..."
      banned k-mers include "AAAGGA", "AAGGAA", "AGGAAA", "GGAAAG", "GAAAGG"
    """
    banned: set = set()
    for m in expected_list or []:
        L = len(m)
        if L == 0 or L >= k:
            continue

        # length needed to see all k-mers crossing at least one junction:
        # we want at least k+L-1 bases; r is a ceil division
        r = (k + L - 1 + L - 1) // L
        rep = m * max(2, r)  # at least 2 repeats, usually more

        for i in range(0, len(rep) - k + 1):
            banned.add(rep[i:i + k])

    return banned


# ---------------------------------------------------------------------------
# Core analysis logic for a single allele of a TRGT record
# ---------------------------------------------------------------------------

def analyze(rec: Dict, cfg_row: Dict) -> Dict:
    """
    Analyze one TRGT allele according to its config row.

    Steps:
      - Determine motif sizes to use (from config or inferred from MOTIFS).
      - For each k:
          * If there is a single motif size and expected motifs:
              - use frame-based counting over all phases (0..k-1),
                accumulate only unexpected motifs.
          * Otherwise:
              - use sliding window.
          * Apply junction filter (ban k-mers from perfect repeats of expected motifs).
          * Apply simple filters (startswith, keep_only_unexpected, count_min).
          * Optionally group by rotational equivalence.
    """
    seq = rec["alt_seq"]
    if not seq:
        return {
            "chrom": rec["chrom"],
            "pos": rec["pos"],
            "trid": rec["trid"],
            "allele_index": rec["allele_index"],
            "blocks": [],
            "expected": rec.get("motifs_from_vcf", []),
            "note": "Empty ALT after cleaning",
        }

    # Determine motif sizes:
    sizes: List[int] = cfg_row["sizes"][:]
    if not sizes:
        # If no sizes in config, infer from expected motifs if present,
        # otherwise use a default small set.
        if rec["motifs_from_vcf"]:
            sizes = sorted({len(m) for m in rec["motifs_from_vcf"]})
        else:
            sizes = [3, 5, 6]

    exp_list: List[str] = rec["motifs_from_vcf"][:] if rec["motifs_from_vcf"] else []
    exp_sets = expected_sets(exp_list) if exp_list else []

    blocks: List[Dict] = []

    for k in sizes:
        # MODE 1: simple repeat with a single motif size and defined expected motifs.
        # We do frame-based counting over all phases, and accumulate only unexpected motifs.
        if len(sizes) == 1 and exp_sets:
            mode = "frame,phases=0..{}".format(k - 1)
            counts: Counter = Counter()
            for ph in range(k):
                cph = count_frame(seq, k, ph)
                for m, n in cph.items():
                    # Only accumulate motifs that are not rotational variants of expected motifs
                    if not is_expected(m, exp_sets):
                        counts[m] += n
        else:
            # MODE 2: more complex / mosaic repeat or no expected motifs.
            # Use sliding window.
            counts = count_sliding(seq, k)
            mode = "sliding"

        # Junction filter: remove k-mers that can be explained by
        # perfect concatenation of shorter expected motifs.
        if exp_list:
            banned = build_banned_kmers_from_expected(exp_list, k)
            if banned:
                for bm in list(counts.keys()):
                    if bm in banned:
                        del counts[bm]

        # Simple filters: startswith, keep_only_unexpected, count_min
        filtered: Dict[str, int] = {}
        for m, n in counts.items():
            # First-letter constraint, if any
            if cfg_row["startswith"] != "ANY" and not m.startswith(cfg_row["startswith"]):
                continue

            # Keep only motifs that are not rotational variants of expected motifs
            if cfg_row["keep_only_unexpected"] and exp_sets and is_expected(m, exp_sets):
                continue

            # Minimal count threshold
            if n >= cfg_row["count_min"]:
                filtered[m] = n

        # Build rows for UI/report, with or without rotational grouping
        if cfg_row["group_rotations"] and filtered:
            groups = group_by_rotations(filtered)
            rows = [
                {
                    "k": k,
                    "rep": rep,
                    "member_counts": member_counts,  # motif -> count
                }
                for (rep, member_counts) in groups
            ]
        else:
            # No grouping: still use the same structure (member_counts) for consistency
            rows = [
                {
                    "k": k,
                    "rep": m,
                    "member_counts": {m: n},
                }
                for m, n in sorted(filtered.items(), key=lambda x: (-x[1], x[0]))
            ]

        blocks.append({"k": k, "mode": mode, "rows": rows})

    return {
        "chrom": rec["chrom"],
        "pos": rec["pos"],
        "trid": rec["trid"],
        "allele_index": rec["allele_index"],
        "blocks": blocks,
        "expected": exp_list,
    }


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------

def all_reports_sorted(reports: List[Dict]) -> List[Dict]:
    """
    Sort reports by (chrom, numeric pos, allele_index, TRID) for stable output.
    """
    def key(rep: Dict):
        chrom = rep["chrom"]
        pos_str = rep["pos"]
        try:
            posi = int(pos_str)
        except ValueError:
            posi = 0
        allele_idx = rep.get("allele_index", 1)
        trid = rep["trid"]
        return (chrom, posi, allele_idx, trid)

    return sorted(reports, key=key)


def write_report(reports: List[Dict], out_path: str) -> None:
    """
    Write a human-readable text report to out_path.

    Each block looks like:
        === chr19:45770204 [DMPK] allele 1 ===
        Expected from VCF: CAG
          - k=3 (frame,phases=0..2):
              GCG <- GCG : 18, GGC : 17, CGG : 15
              CAC: 12
              ACA: 11
    """
    with open(out_path, 'w') as out:
        for rep in all_reports_sorted(reports):
            chrom = rep["chrom"]
            pos = rep["pos"]
            trid = rep["trid"]
            allele_idx = rep.get("allele_index", 1)

            header_loc = f"{chrom}:{pos} [${trid}]"
            # Slightly nicer header, with allele information
            out.write(f"=== {chrom}:{pos} [{trid}] allele {allele_idx} ===\n")

            exp = rep.get("expected") or []
            out.write(f"Expected from VCF: {', '.join(exp) if exp else '(none)'}\n")

            if "note" in rep:
                out.write(f"Note: {rep['note']}\n\n")
                continue

            for b in rep["blocks"]:
                out.write(f"  - k={b['k']} ({b['mode']}):\n")
                if not b["rows"]:
                    out.write("      (no motifs after filtering)\n")
                else:
                    for r in b["rows"]:
                        rep_motif = r["rep"]
                        mc: Dict[str, int] = r["member_counts"]

                        if len(mc) == 1:
                            # Single motif in the group -> simple line
                            m = next(iter(mc.keys()))
                            out.write(f"      {m}: {mc[m]}\n")
                        else:
                            # Multiple motifs in the rotational group:
                            # show detail per motif instead of a single total.
                            parts = [f"{m} : {mc[m]}" for m in mc.keys()]
                            out.write(f"      {rep_motif} <- " + ", ".join(parts) + "\n")
            out.write("\n")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    """
    Command-line entry point.

    Usage:
        python trieky.py <targets.tsv> <input.vcf> <output.txt>
    """
    if len(sys.argv) != 4:
        print("Usage: python trieky.py <targets.tsv> <input.vcf> <output.txt>")
        sys.exit(1)

    cfg_path, vcf_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    cfg = load_config(cfg_path)
    reports: List[Dict] = []

    for rec in iter_trgt_records(vcf_path):
        trid = rec["trid"]
        if trid not in cfg:
            # This TRID is not in the config, skip it
            continue
        reports.append(analyze(rec, cfg[trid]))

    if not reports:
        print("No matching TRID from config found in VCF.")
    write_report(reports, out_path)


if __name__ == "__main__":
    main()

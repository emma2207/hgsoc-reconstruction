#!/usr/bin/env python3

"""Normalize samtools mpileup lines for Conpair's parser expectations.

Conpair expects pileup rows with columns:
  CHROM POS REF BASES BASE_QUALS
and assumes BASES and BASE_QUALS are aligned character-by-character.

samtools mpileup emits an additional depth column and may include special
symbols in BASES (^, $, +/-indels, etc.) that are not quality-bearing.
This script strips/decodes those symbols and emits aligned 5-column rows.
"""

from __future__ import annotations

import sys
from typing import List, Optional, Tuple


def normalize_bases_and_quals(ref: str, bases: str, quals: str) -> Tuple[str, str]:
    ref_base = ref.upper()
    out_bases = []
    out_quals = []

    i = 0
    q = 0
    n = len(bases)

    while i < n:
        c = bases[i]

        if c == "^":
            # Start of read segment; next char is mapping quality marker.
            i += 2
            continue

        if c == "$":
            # End of read segment.
            i += 1
            continue

        if c in "+-":
            # Indel: +/-<len><seq>; does not consume a base-quality character.
            i += 1
            j = i
            while j < n and bases[j].isdigit():
                j += 1
            if j == i:
                # Malformed token; skip sign to avoid infinite loops.
                continue
            indel_len = int(bases[i:j])
            i = j + indel_len
            continue

        # For base-bearing tokens, consume one quality character if available.
        if q >= len(quals):
            break

        if c in ".,":
            b = ref_base
        else:
            b = c.upper()

        # Keep only canonical bases that Conpair's parser explicitly consumes.
        if b in {"A", "C", "G", "T"}:
            out_bases.append(b)
            out_quals.append(quals[q])

        q += 1
        i += 1

    return "".join(out_bases), "".join(out_quals)


def choose_bases_and_quals(fields: List[str]) -> Optional[Tuple[str, str]]:
    """Support multiple pileup layouts.

    1) Standard samtools mpileup:
       CHROM POS REF DEPTH BASES QUALS [...]
    2) Conpair example-like layout:
       CHROM POS REF BASES QUALS [0 VERBOSE ...]
    """

    if len(fields) < 5:
        return None

    # Standard mpileup has integer depth in column 4 (0-based index 3).
    if len(fields) >= 6 and fields[3].isdigit():
        return fields[4], fields[5]

    # Conpair-style: bases and qualities directly after REF.
    return fields[3], fields[4]


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: conpair_normalize_pileup.py <input.mpileup> <output.pileup>", file=sys.stderr)
        return 2

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path, "r", encoding="utf-8", errors="replace") as fin, open(
        out_path, "w", encoding="utf-8"
    ) as fout:
        for line in fin:
            if not line.strip() or line.startswith("[REDUCE RESULT]"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                # Fallback for whitespace-separated lines.
                fields = line.rstrip("\n").split()
            if len(fields) < 6:
                continue

            chrom, pos, ref = fields[0], fields[1], fields[2]
            selected = choose_bases_and_quals(fields)
            if selected is None:
                continue

            bases, quals = selected
            norm_bases, norm_quals = normalize_bases_and_quals(ref, bases, quals)

            # Extra safety: enforce strict alignment length.
            if len(norm_bases) != len(norm_quals):
                m = min(len(norm_bases), len(norm_quals))
                norm_bases = norm_bases[:m]
                norm_quals = norm_quals[:m]

            if not norm_bases:
                # Omit empty rows; Conpair treats absent markers as missing.
                continue

            fout.write(f"{chrom}\t{pos}\t{ref}\t{norm_bases}\t{norm_quals}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

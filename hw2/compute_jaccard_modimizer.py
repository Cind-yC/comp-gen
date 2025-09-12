#!/usr/bin/env python3
"""
compute_jaccard_modimizer.py

Compute Jaccard similarity between two sequences using only
k-mers that satisfy the "modimizer" condition:
    (zlib.crc32(kmer) & 0xffffffff) % m == 0

Outputs: jaccard, ani_exact, ani_approx, n_modimizers_a, n_modimizers_b
"""

import argparse
import math
import sys
import zlib

DNA = ("A", "C", "G", "T")

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute Jaccard/ANI using modimizer-sampled k-mers."
    )
    p.add_argument("-a", "--fasta_a", required=True, help="FASTA A")
    p.add_argument("-b", "--fasta_b", required=True, help="FASTA B")
    p.add_argument("-k", "--kmer", type=int, required=True, help="k-mer length")
    p.add_argument("-m", "--mod", type=int, required=True,
                   help="modimizer modulus (keep kmers with crc32(kmer) %% m == 0)")
    return p.parse_args()

def read_fasta(path):
    """
    Generator of (header, sequence_string) over all records in a FASTA.
    Returns sequence uppercased.
    """
    header = None
    chunks = []
    try:
        with open(path, "r") as fh:
            for line in fh:
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(chunks).upper()
                    header = line.strip()
                    chunks = []
                else:
                    chunks.append(line.strip())
            if header is not None:
                yield header, "".join(chunks).upper()
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: FASTA not found: {path}\n")
        sys.exit(1)

def clean_seq(seq: str) -> str:
    """Uppercase (already done) and convert non-ACGT to 'N'."""
    return "".join(ch if ch in DNA else "N" for ch in seq)

# new function 
def modimizer_kmers(seq: str, k: int, m: int) -> set:
    """
    Return the set of k-mers from seq that:
      - contain only A/C/G/T (no 'N')
      - satisfy (zlib.crc32(kmer) & 0xffffffff) % m == 0
    """
    if k <= 0 or len(seq) < k:
        return set()
    s = clean_seq(seq)
    out = set()
    for i in range(len(s) - k + 1):
        kmer = s[i:i+k]
        if "N" in kmer:
            continue
        h = zlib.crc32(kmer.encode("utf-8")) & 0xffffffff
        if (h % m) == 0:
            out.add(kmer)
    return out

def concat_seqs_from_fasta(path: str) -> str:
    """Concatenate all sequences from a FASTA file."""
    return "".join(seq for _h, seq in read_fasta(path))

def jaccard_sets(A: set, B: set) -> float:
    if not A and not B:
        return 0.0
    inter = len(A & B)
    union = len(A | B)
    return inter / union if union else 0.0

def ani_from_jaccard(j: float, k: int):
    """
    exact ANI = j ** (1/k)
    approx ANI ≈ 1 + (ln j)/k   (good when ANI ~ 1)
    """
    if j <= 0.0:
        return 0.0, 0.0
    exact = j ** (1.0 / k)
    approx = 1.0 + (math.log(j) / k)
    approx = max(0.0, min(1.0, approx))  # clamp
    return exact, approx

def main():
    args = parse_args()
    if args.kmer <= 0 or args.mod <= 0:
        sys.stderr.write("ERROR: k and m must be > 0\n")
        sys.exit(1)

    # load and modimizer-sample
    seq_a = concat_seqs_from_fasta(args.fasta_a)
    seq_b = concat_seqs_from_fasta(args.fasta_b)

    A = modimizer_kmers(seq_a, args.kmer, args.mod)
    B = modimizer_kmers(seq_b, args.kmer, args.mod)

    j = jaccard_sets(A, B)
    exact, approx = ani_from_jaccard(j, args.kmer)

    # Output single line (easy to redirect to .txt)
    print(
        f"a={args.fasta_a} b={args.fasta_b} k={args.kmer} m={args.mod} "
        f"jaccard={j:.6f} ani_exact={exact:.6f} ani_approx={approx:.6f} "
        f"n_modimizers_a={len(A)} n_modimizers_b={len(B)}"
    )

if __name__ == "__main__":
    main()

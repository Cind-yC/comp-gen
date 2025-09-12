#!/usr/bin/env python3
"""
Q2a: Build a graph from reads for k=3.

Modes:
  --nodes kminus1  -> classic de Bruijn: nodes=(k-1)-mers, edges=k-mers
  --nodes kmer     -> k-mer overlap graph: nodes=k-mers, edges when suffix/prefix overlap by k-1

Edges show multiplicity; in de Bruijn mode, edges also carry the k-mer label.

Render:
  dot   -T png debruijn_k3.dot -o debruijn_k3.png
  neato -T png debruijn_k3.dot -o debruijn_k3.png
"""

import argparse
from collections import Counter, defaultdict
import sys

DEFAULT_READS = [
    "ATTCA",
    "ATTGA",
    "CATTG",
    "CTTAT",
    "GATTG",
    "TATTT",
    "TCATT",
    "TCTTA",
    "TGATT",
    "TTATT",
    "TTCAT",
    "TTCTT",
    "TTGAT",
]

def parse_args():
    ap = argparse.ArgumentParser(description="Build de Bruijn / k-mer overlap graph to Graphviz DOT.")
    ap.add_argument("-f", "--reads_file", help="Optional file with one read per line. If omitted, uses built-in reads.")
    ap.add_argument("-k", "--kmer", type=int, default=3, help="k-mer length (default: 3)")
    ap.add_argument("-o", "--out", default="debruijn_k3.dot", help="Output DOT filename (default: debruijn_k3.dot)")
    ap.add_argument("--nodes", choices=["kminus1", "kmer"], default="kminus1",
                    help="Node type: 'kminus1' (classic de Bruijn) or 'kmer' (k-mer overlap graph). Default: kminus1")
    return ap.parse_args()

def load_reads(path=None):
    if not path:
        return DEFAULT_READS[:]
    reads = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s:
                reads.append(s.upper())
    return reads

def build_debruijn_kminus1(reads, k):
    """
    Classic de Bruijn:
      nodes: (k-1)-mers
      edges: u -> v for each k-mer, where u = kmer[:-1], v = kmer[1:]
    Returns (nodes, edges, edge_labels) where:
      nodes = set of strings
      edges = Counter((u,v) -> count)
      edge_labels = dict((u,v) -> set of kmers encountered)  (for labeling)
    """
    if k < 2:
        raise ValueError("k must be >= 2")
    nodes = set()
    edges = Counter()
    edge_labels = defaultdict(Counter)

    for read in reads:
        if len(read) < k: 
            continue
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            u = kmer[:-1]
            v = kmer[1:]
            nodes.update([u, v])
            edges[(u, v)] += 1
            edge_labels[(u, v)][kmer] += 1

    return nodes, edges, edge_labels

def build_kmer_overlap(reads, k):
    """
    k-mer overlap graph:
      nodes: k-mers
      edges: k1 -> k2 if suffix(k1, k-1) == prefix(k2, k-1) AND the pair occurs consecutively in a read.
    Returns (nodes, edges, edge_labels) where edge_labels can be empty or hold overlap info.
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    nodes = set()
    edges = Counter()
    edge_labels = defaultdict(Counter)

    for read in reads:
        if len(read) < k:
            continue
        kmers = [read[i:i+k] for i in range(len(read) - k + 1)]
        for km in kmers:
            nodes.add(km)
        for i in range(len(kmers) - 1):
            k1 = kmers[i]
            k2 = kmers[i+1]
            # These are consecutive in the read, so they necessarily overlap by k-1.
            edges[(k1, k2)] += 1
            # Optionally store the shared overlap (k-1-mer) if you want:
            overlap = k1[1:]  # == k2[:-1]
            edge_labels[(k1, k2)][overlap] += 1

    return nodes, edges, edge_labels

def write_dot(nodes, edges, edge_labels, out_path, mode):
    with open(out_path, "w") as out:
        out.write("digraph debruijn {\n")
        out.write('  graph [overlap=false, splines=true];\n')
        out.write('  node  [shape=box, style="rounded,filled", fillcolor="#eef5ff", fontname="Helvetica"];\n')
        out.write('  edge  [fontname="Helvetica"];\n')

        for n in sorted(nodes):
            out.write(f'  "{n}";\n')

        for (u, v), w in sorted(edges.items()):
            # label multiplicity; in de Bruijn mode also include the k-mer(s)
            if mode == "kminus1":
                # collect up to a few representative k-mers for readability
                kmers = ", ".join(list(edge_labels[(u, v)].keys())[:3])
                label_parts = []
                if kmers:
                    label_parts.append(kmers)
                if w > 1:
                    label_parts.append(f"x{w}")
                label = "\\n".join(label_parts) if label_parts else ""
            else:
                # kmer node mode: just multiplicity (edges are transitions)
                label = f"{w}" if w > 1 else ""

            if label:
                out.write(f'  "{u}" -> "{v}" [label="{label}"];\n')
            else:
                out.write(f'  "{u}" -> "{v}";\n')

        out.write("}\n")

def main():
    args = parse_args()
    reads = load_reads(args.reads_file)
    if not reads:
        sys.stderr.write("No reads found.\n")
        sys.exit(1)

    if args.nodes == "kminus1":
        nodes, edges, edge_labels = build_debruijn_kminus1(reads, args.kmer)
    else:  # args.nodes == "kmer"
        nodes, edges, edge_labels = build_kmer_overlap(reads, args.kmer)

    write_dot(nodes, edges, edge_labels, args.out, args.nodes)
    sys.stderr.write(
        f"Wrote DOT to {args.out} with {len(nodes)} nodes and {len(edges)} edges "
        f"(mode={args.nodes}, k={args.kmer}).\n"
    )

if __name__ == "__main__":
    main()

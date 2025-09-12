"""
read in fasta file
introduce specified level of mutations
output new fasta file

in command line: 
- input filename (-i)
- output filename (-o)
- level of mutations (-m)
- random seed (-s)

only introduce substitutions by changing the input base to a 
different base at random
substitutions should be inserted at random locations at a 
specified rate
"""
import argparse
import random
import sys
import math

DNA = ('A', 'C', 'G','T') # tuple in python - fixed, never changes, just for lookup

def parse_args(): 
    parser = argparse.ArgumentParser(
        description="Either mutate FASTA or compute Jaccard/ANI."
    )
    
    # mutate mode 
    parser.add_argument("-i", "--input", help="Input FASTA filename")
    parser.add_argument("-o", "--output", help="Output FASTA filename")
    parser.add_argument("-m", "--mutation_rate", type=float, 
                        help="Mutation rate (e.g., 0.015 for 1.5%% of bases)")
    parser.add_argument("-s", "--seed", type=int, 
                        help="Random seed for reproducibility")
    parser.add_argument("--wrap", type=int, default=60, 
                        help="Line wrap width for FASTA output (default: 60; use 0 for no wrap)") # for better formatting 
   
    # jaccard mode 
    parser.add_argument("-a", "--fasta_a", help="FASTA A (jaccard mode)")
    parser.add_argument("-b", "--fasta_b", help="FASTA B (jaccard mode)")
    parser.add_argument("-k", "--kmer", type=int, help="k-mer length (jaccard mode)")
   
    args = parser.parse_args() # collects flags/values & puts them into the args object 
    # e.g. args.input = "input.fa", args.mutation_rate = 0.01, etc. 

    # decide mode 
    jaccard_mode = (args.fasta_a is not None or args.fasta_b is not None or args.kmer is not None)

    # ignore error handling for now 
    # choose jaccard_mode
    if jaccard_mode:
        if not (args.fasta_a and args.fasta_b and args.kmer):
            parser.error("Jaccard mode requires 2 fasta files & kmer")
        if args.kmer <= 0:
            parser.error("kmer length must be > 0.")
        return args
    
    # choose mutate mode:
    if not (args.input and args.output and args.mutation_rate is not None and args.seed is not None):
        parser.error("mutate mode requires input, output, rate and seed.")
    if not(0.0 <= args.mutation_rate <= 1.0):
        parser.error("Mutation rate must be between 0.0 and 1.0")
    return args

def read_fasta(path):
    """
    generator - yields (header, sequence_string) for each FASTA record
    sequences translated to uppercase for simpler mutation logic
    """
    header=None
    chunks = []
    try: 
        with open(path, "r") as fh: 
            for line in fh:
                if not line:
                    continue
                if line.startswith(">"): # if hitting header line 
                    if header is not None: # if already have previous header + sequence stored, emit the record 
                        yield header, "".join(chunks).upper() # how to visualize this? 
                    header = line.strip() # store new header
                    chunks = [] # reset new sequence storage for new record
                else: 
                    chunks.append(line.strip()) # if it's not a header but part of DNA sequence, add to chunks
            if header is not None: # ensures last record in the file is also output -- but how to visualize this? 
                yield header, "".join(chunks).upper()
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: Input FASTA '{path}' not found.\n")
        sys.exit(1)

def wrap_seq(seq, width):
    """
    return seq wrapped to fixed line width; if width <= 0, return as a single line
    """
    if width and width > 0: 
        return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))
    else:
        return seq

def mutate_sequence(seq, rate, rng):
    """
    return mutated version of seq
    """
    seq_list = list(seq)

    mutable_idxs = [i for i, b in enumerate(seq_list) if b in DNA]
    n_mutable = len(mutable_idxs)

    # handles the case of mutation rate = 0
    if n_mutable == 0 or rate == 0.0:
        return seq
    
    n_to_mutate = int(round(rate * n_mutable))
    # handles the case of rate too low 
    if n_to_mutate == 0:
        return seq
    
    # randomly select n_to_mutate unique indices from list mutable_idxs
    chosen_positions = rng.sample(mutable_idxs, n_to_mutate) # random number generator 


    for i in chosen_positions:
        orig = seq_list[i]
        choices = [b for b in DNA if b != orig] # make a list of 3 bases diff from orig 
        seq_list[i] = rng.choice(choices) # randomly pick one base out of the 3 & replace
    
    return "".join(seq_list)

def clean_seq(seq: str) -> str:
    """
    upper case, map non-ACGT to N
    """
    s = seq.upper()
    return "".join(ch if ch in DNA else "N" for ch in s)

def kmer_set(seq:str, k:int) -> set:
    """
    return set of all k-mers (length k) from seq that contain only ACGT. any kmer w N is skipped
    """
    if k <= 0 or len(seq) < k:
        return set()
    s = clean_seq(seq)
    out = set()
    for i in range(len(s) - k + 1):
        kmer = s[i:i+k]
        if "N" in kmer:
            continue 
        out.add(kmer)
    return out

def jaccard(path_a:str, path_b:str, k:int) -> float:
    """
    compute jaccard between kmer sets
    """
    seq_a = "".join(seq for _h, seq in read_fasta(path_a))
    seq_b = "".join(seq for _h, seq in read_fasta(path_b))

    A = kmer_set(seq_a, k)
    B = kmer_set(seq_b, k)
    if not A and not B:
        return 0.0
    inter = len(A & B)
    union = len(A | B)
    return inter/union if union else 0.0

def ani(j:float,k:int) -> tuple[float, float]:
    """
    return exact ani and approx ani from jaccard and k 
    exact ANI = J**(1/k)
    approx ANI = 1 + (ln J)/k
    """
    if j <= 0.0:
        return 0.0, 0.0
    exact = j ** (1.0 / k)
    approx = 1.0 + (math.log(j) / k)
    approx = max(0.0, min(1.0, approx))
    return exact, approx


def main():
    args = parse_args()
    # jaccard mode
    if args.fasta_a and args.fasta_b and args.kmer:
        j = jaccard(args.fasta_a, args.fasta_b, args.kmer)
        exact, approx = ani(j, args.kmer)
        print(f"a={args.fasta_a} b={args.fasta_b} k={args.kmer} jaccard={j:.6f} ani_exact={exact:.6f} ani_approx={approx:.6f}")
        return

    # mutate mode 
    rng = random.Random(args.seed)
    records = list(read_fasta(args.input))
    if not records:
        sys.stderr.write("ERROR: No FASTA records found in input.\n")
        sys.exit(1)

    with open(args.output, "w") as out:
        for header, seq in records:
            mutated = mutate_sequence(seq, args.mutation_rate, rng)
            out.write(f"{header}\n")
            out.write(f"{wrap_seq(mutated, args.wrap)}\n")
    sys.stderr.write(f"Wrote mutated FASTA to {args.output}\n")

if __name__ == "__main__":
    main()
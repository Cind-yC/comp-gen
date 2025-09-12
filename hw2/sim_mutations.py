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

DNA = ('A', 'C', 'G','T') # tuple in python - fixed, never changes, just for lookup

def parse_args(): 
    parser = argparse.ArgumentParser(
        description="Introduce random substitution mutations into a FASTA file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA filename")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA filename")
    parser.add_argument("-m", "--mutation_rate", type=float, required=True,
                        help="Mutation rate (e.g., 0.015 for 1.5%% of bases)")
    parser.add_argument("-s", "--seed", type=int, required=True,
                        help="Random seed for reproducibility")
    parser.add_argument("--wrap", type=int, default=60, 
                        help="Line wrap width for FASTA output (default: 60; use 0 for no wrap)") # for better formatting 
    args = parser.parse_args() # collects flags/values & puts them into the args object 
    # e.g. args.input = "input.fa", args.mutation_rate = 0.01, etc. 

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

def main():
    args = parse_args()
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
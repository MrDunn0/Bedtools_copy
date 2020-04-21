import parse_args

args = parse_args()
if args.command == 'getfasta':
    if args.fi:
        print(args.fi)


# open and parse fasta file
def parse_fasta(file):
    sequences = {}
    with open(file, 'r') as in_f:
        sequence_name = None
        sequence = []
        for line in in_f:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences[sequence_name] = "".join(sequence)
                    sequence = []
                sequence_name = line[1:]
            else:
                sequence.append(line)
        sequences[sequence_name] = "".join(sequence)
    return sequences
from parse_args import parse_args
import getfasta
import sort

if __name__ == "__main__":
    args = parse_args()
    if args.command == 'sort':
        if args.sizeA:
            sort.write_for_size(sort.sizeA(args.i))
        elif args.sizeD:
            sort.write_for_size(sort.sizeD(args.i))
        elif args.chrThenSizeA:
            sort.write_for_chr_sorted(sort.chrThenSizeA(args.i))
        elif args.chrThenSizeD:
            sort.write_for_chr_sorted(sort.chrThenSizeD(args.i))
        elif args.chrThenScoreA:
            sort.write_for_chr_sorted(sort.chrThenScoreA(args.i))
        elif args.chrThenScoreD:
            sort.write_for_chr_sorted(sort.chrThenSizeD(args.i))
        else:
            sort.write_for_chr_sorted(sort.sort_by_default(args.i))

    elif args.command == 'subtract':
        pass

    elif args.command == 'merge':
        pass

    elif args.command == 'intersect':
        pass

    elif args.command == 'getfasta':
        input_bed = args.i
        input_fasta = getfasta.parse_sequences(args.fi)
        print(input_fasta)
        getfasta.do_bed(input_bed, input_fasta)
        pass


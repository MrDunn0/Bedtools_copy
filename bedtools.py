<<<<<<< HEAD
import parse_args
import getfasta
=======
from parse_args import parse_args
>>>>>>> sort_fun
import sort

if __name__ == "__main__":
    args = parse_args()
    if args.command == 'sort':
        if args.sizeA:
            sort.write_for_size(sort.sizeA(args.i))
        elif args.sizeD:
            sort.write_for_size(sort.sizeD(args.i))
        elif args.chrThenSizeA:
            args.write_for_chr_sorted(args.chrThenSizeA(args.i))
        elif args.chrThenSizeD:
            args.write_for_chr_sorted(args.chrThenSizeD(args.i))
        elif args.chrThenScoreA:
            args.write_for_chr_sorted(args.chrThenScoreA(args.i))
        elif args.chrThenScoreD:
            args.write_for_chr_sorted(args.chrThenSizeD(args.i))
        else:
            sort.write_for_chr_sorted(sort.sort_by_default(args.i))

    elif args.command == 'subtract':
        pass

    elif args.command == 'merge':
        pass

    elif args.command == 'intersect':
        pass

    elif args.command == 'getfasta':
        pass


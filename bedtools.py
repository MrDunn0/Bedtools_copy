import argparse
from main_funcs import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='bed file')
    subparsers = parser.add_subparsers(dest='command')

    # options for sort
    sort = subparsers.add_parser('sort',
                                 help='The bedtools sort tool sorts a feature file \
                                 by chromosome and other criteria.')
    sort.add_argument('--sizeA', action="store_true",
                      help='Sort by feature size in ascending order.')
    sort.add_argument('--sizeD', action="store_true",
                      help='Sort by feature size in descending order.')
    sort.add_argument('--chrThenSizeA', action="store_true",
                      help='Sort by chromosome (asc), then by feature size (asc).')
    sort.add_argument('--chrThenSizeD', action="store_true",
                      help='Sort by chromosome (asc), then by feature size (desc).')
    sort.add_argument('--chrThenScoreA', action="store_true",
                      help='Sort by chromosome (asc), then by score (asc).')
    sort.add_argument('--chrThenScoreD', action="store_true",
                      help='Sort by chromosome (asc), then by score (desc).')
    # sort.add_argument('-g', action="store_true",
    #                   help='Define sort order by order of tab-delimited file \
    #                     with chromosome names in the first column.')
    # sort.add_argument('--faidx', action="store_true",
    #                   help='Define sort order by order of tab-delimited file \
    #                     with chromosome names in the first column. Sort by specified chromosome order.')

    # options for subtract
    subtract = subparsers.add_parser('subtract',
                                     help='bedtools subtract bedtools` subtract \
                                     searches for features in B that overlap A by at least the number \
                                     of base pairs given by the -f option.')
    subtract._action_groups.pop()
    required = subtract.add_argument_group('required arguments')
    optional = subtract.add_argument_group('optional arguments')
    required.add_argument('-a', help='A file', required=True)
    required.add_argument('-b', help='B file', required=True)
    optional.add_argument('-f', type=float,
                          help='Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).')
    optional.add_argument('-s', action='store_true',
                          help='Force “strandedness”. That is, only report hits in B that overlap A \
                        on the same strand. By default, overlaps are reported without respect to strand.')
    optional.add_argument('-S', action='store_true',
                          help='Require different strandedness. That is, only report hits in B that overlap A\
                         on the _opposite_ strand. By default, overlaps are reported without respect to strand.')
    optional.add_argument('-A', action='store_true',
                          help='Remove entire feature if any overlap. That is, by default, only subtract \
                        the portion of A that overlaps B. Here, if any overlap is found (or -f amount), \
                        the entire feature is removed.')

    # options for merge

    merge = subparsers.add_parser('merge',
                                  help='bedtools merge combines overlapping or “book-ended” features\
                                  in an interval file into a single feature which spans \
                                  all of the combined features.')
    merge._action_groups.pop()
    required = merge.add_argument_group('required arguments')
    optional = merge.add_argument_group('optional arguments')
    required.add_argument('-i', help='Input BED file', required=True)
    optional.add_argument('-d', type=int,
                          help='Maximum distance between features allowed for features to be merged. \
                        Default is 0. That is, overlapping and/or book-ended features are merged.')

    # options for intersect

    intersect = subparsers.add_parser("intersect",
                                      help='bedtools intersect allows one to screen \
                                      for overlaps between two sets of genomic features. ')
    intersect._action_groups.pop()
    required = intersect.add_argument_group('required arguments')
    optional = intersect.add_argument_group('optional arguments')

    required.add_argument('-a', help='A file', required=True)
    required.add_argument('-b', help='B file', required=True)
    optional.add_argument('-s', action='store_true',
                          help='Force strandedness. That is, only merge features that are the same strand. \
                        By default, this is disabled.')
    optional.add_argument('-S', action='store_true',
                          help='Force merge for one specific strand only. Follow with + or - to force merge \
                        from only the forward or reverse strand, respectively. \
                        By default, merging is done without respect to strand.')

    # options for getfasta
    getfasta = subparsers.add_parser("getfasta",
                                     help='bedtools getfasta extracts sequences from a FASTA file'
                                          'for each of the intervals defined in a BED file.')
    getfasta.add_argument('--fi', '-fi',
                          help='Specify an output file name.')
    getfasta.add_argument('--fo', '-fo', default=False,
                          help='Specify an output file name.'
                               'By default, output goes to stdout.')
    getfasta.add_argument('--bed', '-bed',
                          help='BED file of ranges to extract from -fi')
    getfasta.add_argument('--tab', '-tab', action='store_true', default=False,
                          help='Report extract sequences in a tab-delimited format'
                               'instead of in FASTA format.')
    # getfasta.add_argument('--name', action='store_true',
    #                       help='Use the “name” column in the BED file \
    #                       for the FASTA headers in the output FASTA file.')
    # getfasta.add_argument('--bedOut', action='store_true',
    #                       help='Report extract sequences in a tab-delimited BED format \
    #                       instead of in FASTA format.')
    # getfasta.add_argument('--s', action='store_true',
    #                       help='Force strandedness. If the feature occupies the antisense strand, \
    #                       the sequence will be reverse complemented. \
    #                       Default: strand information is ignored.')
    # getfasta.add_argument('--split', action='store_true',
    #                       help='Given    BED12 input, extract and concatenate the sequences \
    #                       from the BED “blocks” (e.g., exons)')

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.command == 'sort':
        if args.sizeA:
            write_for_size(sizeA(args.i))
        elif args.sizeD:
            write_for_size(sizeD(args.i))
        elif args.chrThenSizeA:
            write_for_chr_sorted(chrThenSizeA(args.i))
        elif args.chrThenSizeD:
            write_for_chr_sorted(chrThenSizeD(args.i))
        elif args.chrThenScoreA:
            write_for_chr_sorted(chrThenScoreA(args.i))
        elif args.chrThenScoreD:
            write_for_chr_sorted(chrThenSizeD(args.i))
        else:
            write_for_chr_sorted(sort_by_default(args.i))

    elif args.command == 'getfasta':

        # define input bed
        if args.i:
            input_bed = args.i
        if args.bed:
            input_bed = args.bed
        # define input fasta
        input_fasta = parse_sequences(args.fi)
        # define output

        # our tuples with names and sequences to return
        names_and_seqs = list(get_sequences(input_bed, input_fasta))

        # fasta or tab, stdout or in file
        give_me_fasta(names_and_seqs, file_out=args.fo, tab_out=args.tab)

    else:
        files_are_sorted = [True, True]
        if args.i:
            if not sorting_check(args.i):
                files_are_sorted[1] = False
                print("Please, sort your BED file first, then try again.")
        else:
            if not sorting_check(args.a):
                files_are_sorted[1] = False
                print(f"Your file {args.a} is not sorted. Please, sort it, then try again.")
            if not sorting_check(args.b):
                files_are_sorted[1] = False
                print(f"Your file {args.b} is not sorted. Please, sort it, then try again.")

        if all(files_are_sorted):

            if args.command == "subtract":
                overlaps_collector(args.a, args.b, mode="subtract", s=args.s, S=args.S, A=args.A, f=args.f)

            elif args.command == 'intersect':
                overlaps_collector(args.a, args.b, mode="intersect", s=args.s, S=args.S)

            elif args.command == 'merge':
                main_merge(args.i, d=args.d)

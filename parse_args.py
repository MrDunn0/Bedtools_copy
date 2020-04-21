#!/usr/bin/python
# -*- coding: UTF-8 -*-


import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='bed file')
    subprasers = parser.add_subparsers(dest='command')

    # options for sort
    sort = subprasers.add_parser('sort',
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

    # optoins for subtract
    subtract = subprasers.add_parser('subtract',
                                     help='bedtools subtract bedtools` subtract \
                                     searches for features in B that overlap A by at least the number \
                                     of base pairs given by the -f option.')
    subtract.add_argument('-f', action='store_true',
                          help='Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).')
    subtract.add_argument('-F', action='store_true',
                          help='Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).')
    subtract.add_argument('-r', action='store_true',
                          help='Require that the fraction of overlap be reciprocal for A and B. \
                        In other words, if -f is 0.90 and -r is used, this requires that B \
                        overlap at least 90% of A and that A also overlaps at least 90% of B.')
    subtract.add_argument('-e', action='store_true',
                          help='Require that the minimum fraction be satisfied for A _OR_ B. \
                        In other words, if -e is used with -f 0.90 and -F 0.10 this requires \
                        that either 90% of A is covered OR 10% of B is covered. Without -e, \
                        both fractions would have to be satisfied.')
    subtract.add_argument('-s', action='store_true',
                          help='Force “strandedness”. That is, only report hits in B that overlap A \
                        on the same strand. By default, overlaps are reported without respect to strand.')
    subtract.add_argument('-S', action='store_true',
                          help='Require different strandedness. That is, only report hits in B that overlap A\
                         on the _opposite_ strand. By default, overlaps are reported without respect to strand.')
    subtract.add_argument('-A', action='store_true',
                          help='Remove entire feature if any overlap. That is, by default, only subtract \
                        the portion of A that overlaps B. Here, if any overlap is found (or -f amount), \
                        the entire feature is removed.')
    subtract.add_argument('-N', action='store_true',
                          help='Same as -A except when used with -f, the amount is the sum of all features \
                        (not any single feature).')

    # options for merge
    merge = subprasers.add_parser('merge',
                                  help='bedtools merge combines overlapping or “book-ended” features\
                                  in an interval file into a single feature which spans \
                                  all of the combined features.')
    merge.add_argument('-s', action='store_true',
                       help='Force strandedness. That is, only merge features that are the same strand. \
                        By default, this is disabled.')
    merge.add_argument('-S', action='store_true',
                       help='Force merge for one specific strand only. Follow with + or - to force merge \
                        from only the forward or reverse strand, respectively. \
                        By default, merging is done without respect to strand.')
    merge.add_argument('-d', action='store_true',
                       help='Maximum distance between features allowed for features to be merged. \
                        Default is 0. That is, overlapping and/or book-ended features are merged.')
    merge.add_argument('-c', action='store_true',
                       help='Specify columns from the input file to operate upon (see -o option, below). \
                        Multiple columns can be specified in a comma-delimited list.')
    merge.add_argument('--header', action='store_true',
                       help='Print the header from the A file prior to results.')

    # options for intersect
    intersect = subprasers.add_parser("intersect",
                                      help='bedtools intersect allows one to screen \
                                      for overlaps between two sets of genomic features. ')
    intersect.add_argument('--some_arg', action='store_true',
                           help='')

    # options for getfasta
    getfasta = subprasers.add_parser("getfasta",
                                     help='bedtools getfasta extracts sequences from a FASTA file\
                                     for each of the intervals defined in a BED file.')
    getfasta.add_argument('--fi', action='store_true',
                          help='Specify an output file name.')
    getfasta.add_argument('--fo', action='store_true',
                          help='Specify an output file name. \
                          By default, output goes to stdout.')

    getfasta.add_argument('--name', action='store_true',
                          help='Use the “name” column in the BED file \
                          for the FASTA headers in the output FASTA file.')
    getfasta.add_argument('--tab', action='store_true',
                          help='Report extract sequences in a tab-delimited format \
                          instead of in FASTA format.')
    getfasta.add_argument('--bedOut', action='store_true',
                          help='Report extract sequences in a tab-delimited BED format \
                          instead of in FASTA format.')
    getfasta.add_argument('--s', action='store_true',
                          help='Force strandedness. If the feature occupies the antisense strand, \
                          the sequence will be reverse complemented. \
                          Default: strand information is ignored.')
    getfasta.add_argument('--split', action='store_true',
                          help='Given BED12 input, extract and concatenate the sequences \
                          from the BED “blocks” (e.g., exons)')

    return parser.parse_args()

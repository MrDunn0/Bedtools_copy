# Bedtools_copy

### Descriprion
This is the skin-deep version of Bedtool realized in Python 3.

### Getting started and basic usage

#### Quick start

For calling short help message or check installation type:

`python bedtools.py -h` or `python bedtools.py --help` or `python bedtools.py <TOOLNAME> -h`

For getting started type:

`python bedtools.py <TOOLNAME> <TOOL_ARGUMENT1> <TOOL_ARGUMENT2>`

#### Implemented tools

List of tools, which are implemented in this Bedtools realization:

`sort` sorts a feature file by chromosome and other criteria

`getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED file

`intersect` allows one to screen for overlaps between two sets of genomic features

`subtract` allows you to subtract from the file A all intersections with the file B

`merge` combines overlapping or “book-ended” features in an interval file into a single feature which spans all of the combined features

Since this is a raw version, only a small part of the original parameters is available.

### Input and output files

### Other

Tool was written during the course of Python from Bioinformatic Institute in 2020.

#### Developers:
    
Tatiana Maslikova (https://github.com/Neonbird), Yulia Yakovleva (https://github.com/Yulia-Yakovleva) and Mikhail Ushakov (https://github.com/MrDunn0)

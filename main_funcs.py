b_intervals = list()


# getfasta code


class SeqObject():
    """
    Создаем класс, шоб было как в биопитоне
    """
    def __init__(self, seq_id=None, description=None, name=None, sequence=None):
        self.seq_id = seq_id
        self.description = description
        self.name = name
        self.sequence = sequence

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)


def iterate_on_fasta(file):
    header, sequence = None, []
    for line in file:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield header, ''.join(sequence)
            header, sequence = line, []
        else:
            sequence.append(line)
    if header:
        yield header, ''.join(sequence)


def fasta_tuple_to_object(tuple):
    name = tuple[0].lstrip('>')
    seq_id = name.split(' ')[0]
    description = ' '.join(name.split(' ')[1:])
    sequence = tuple[1]
    seq_obj = SeqObject(name=name, seq_id=seq_id,
                        description=description,
                        sequence=sequence)
    return seq_obj


def parse_sequences(filename):
    with open(filename, 'r') as text_file:
        list_of_objects = []
        fasta_tuples = iterate_on_fasta(text_file)
        for element in fasta_tuples:
            list_of_objects.append(fasta_tuple_to_object(element))
        return list_of_objects


def get_sequences(file, list_of_objects):
    with open(file, 'r') as in_f:
        for line in in_f:
            line = line.strip()
            if line.startswith("track") or line.startswith("#"):
                pass
            else:
                chrom = line.split('\t')[0]
                start = int(line.split('\t')[1])
                end = int(line.split('\t')[2])

                for obj in list_of_objects:
                    if obj.seq_id == chrom:
                        name = f"{chrom}:{start}-{end}"
                        """
                        наш скрипт включает все границы
                        первая позиция = первый нуклеотид
                        десятая позиция = десятый нуклеотид
                        """
                        seq = obj.sequence[start - 1:end]
                        yield name, seq


def give_me_fasta(tuples_to_return, file_out=False, tab_out=False):
    print(file_out)

    if file_out and not tab_out:
        with open(f"{file_out}", 'w') as out_file:
            for name, seq in tuples_to_return:
                out_file.write(f">{name}\n{seq}\n")

    if not file_out and not tab_out:
        for name, seq in tuples_to_return:
            print(f">{name}")
            print(seq)

    if file_out and tab_out:
        with open(f"{file_out}", 'w') as out_file:
            for name, seq in tuples_to_return:
                out_file.write(f"{name}\t{seq}\n")

    if not file_out and tab_out:
        for name, seq in tuples_to_return:
            print(f"{name}\t{seq}")


# Sort code


def get_header_and_start_line(file):
    header = []
    i = 0
    with open(file) as f:
        while True:
            line = f.readline()
            if line.startswith("track") or line.startswith("#"):
                header.append(line)
                i += 1
            else:
                return header, i


def chr_then_start_ordering(line, chr_dict_with_sorted_start, sign):
    chr, start = line[0], line[1]
    if chr not in chr_dict_with_sorted_start:
        chr_dict_with_sorted_start[chr] = [line]
        return chr_dict_with_sorted_start
    else:
        min_pos, max_pos = 0, len(chr_dict_with_sorted_start[chr])
        while max_pos - min_pos > 0:
            med = (min_pos + max_pos) // 2
            if chr_dict_with_sorted_start[chr][med][1] > start:
                max_pos = med
            else:
                min_pos = med + 1
        chr_dict_with_sorted_start[chr] = \
            chr_dict_with_sorted_start[chr][:max_pos] + [line] + chr_dict_with_sorted_start[chr][max_pos:]
        return chr_dict_with_sorted_start


def size_ordering(line, sorted_list, order):
    if not sorted_list:
        sorted_list.append(line)
        return sorted_list
    size = int(line[2]) - int(line[1])
    min_index, max_index = 0, len(sorted_list)
    if order == "asc":
        while max_index - min_index > 0:
            med = (min_index + max_index) // 2
            if int(sorted_list[med][2]) - int(sorted_list[med][1]) > size:
                max_index = med
            else:
                min_index = med + 1
        sorted_list = sorted_list[:max_index] + [line] + sorted_list[max_index:]
    elif order == "desc":
        while max_index - min_index > 0:
            med = (min_index + max_index) // 2
            if int(sorted_list[med][2]) - int(sorted_list[med][1]) < size:
                max_index = med
            else:
                min_index = med + 1
        sorted_list = sorted_list[:max_index] + [line] + sorted_list[max_index:]
    return sorted_list


def apply_fun_to_file(function, file, storage, order="asc"):
    header, index_start_line = get_header_and_start_line(file)
    with open(file) as f:
        for i in range(index_start_line):
            f.readline()
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip().split("\t")
            storage = function(line, storage, order)
    return header, storage


def chr_then_size_ordering(line, sorted_dict, order="asc"):
    chr, size = line[0], int(line[2]) - int(line[1])
    if chr not in sorted_dict:
        sorted_dict[chr] = [line]
        return sorted_dict
    else:
        min_pos, max_pos = 0, len(sorted_dict[chr])
        if order == "asc":
            while max_pos - min_pos > 0:
                med = (min_pos + max_pos) // 2
                if int(sorted_dict[chr][med][2]) - int(sorted_dict[chr][med][1]) > size:
                    max_pos = med
                else:
                    min_pos = med + 1
            sorted_dict[chr] = \
                sorted_dict[chr][:max_pos] + [line] + sorted_dict[chr][max_pos:]
        else:
            while max_pos - min_pos > 0:
                med = (min_pos + max_pos) // 2
                if int(sorted_dict[chr][med][2]) - int(sorted_dict[chr][med][1]) < size:
                    max_pos = med
                else:
                    min_pos = med + 1
            sorted_dict[chr] = \
                sorted_dict[chr][:max_pos] + [line] + sorted_dict[chr][max_pos:]
        return sorted_dict


def chr_then_score_sorting(line, chr_dict_with_sorted_score, order="asc"):
    chr, score = line[0], line[4]
    if chr not in chr_dict_with_sorted_score:
        chr_dict_with_sorted_score[chr] = [line]
        return chr_dict_with_sorted_score
    else:
        min_pos, max_pos = 0, len(chr_dict_with_sorted_score[chr])
        if order == "asc":
            while max_pos - min_pos > 0:
                med = (min_pos + max_pos) // 2
                if chr_dict_with_sorted_score[chr][med][4] > score:
                    max_pos = med
                else:
                    min_pos = med + 1
            chr_dict_with_sorted_score[chr] = \
                chr_dict_with_sorted_score[chr][:max_pos] + [line] + chr_dict_with_sorted_score[chr][max_pos:]
        elif order == "desc":
            while max_pos - min_pos > 0:
                med = (min_pos + max_pos) // 2
                if chr_dict_with_sorted_score[chr][med][4] < score:
                    max_pos = med
                else:
                    min_pos = med + 1
            chr_dict_with_sorted_score[chr] = \
                chr_dict_with_sorted_score[chr][:max_pos] + [line] + chr_dict_with_sorted_score[chr][max_pos:]
        return chr_dict_with_sorted_score


def sort_by_default(file):
    '''by chromosome then start position in ascending'''
    return apply_fun_to_file(chr_then_start_ordering, file, dict())


def sizeA(file):
    '''Sort by feature size in ascending order'''
    return apply_fun_to_file(size_ordering, file, list())


def sizeD(file):
    '''Sort by feature size in descending order'''
    return apply_fun_to_file(size_ordering, file, list(), "desc")


def chrThenSizeA(file):
    ''' Sort by chromosome (asc), then by feature size (asc).'''
    return apply_fun_to_file(chr_then_size_ordering, file, dict())


def chrThenSizeD(file):
    ''' Sort by chromosome (asc), then by feature size (asc).'''
    return apply_fun_to_file(chr_then_size_ordering, file, dict(), "desc")


def chrThenScoreA(file):
    '''Sort by chromosome (asc), then by score (asc).'''
    return apply_fun_to_file(chr_then_score_sorting, file, dict(), "asc")


def chrThenScoreD(file):
    '''Sort by chromosome (asc), then by score (asc).'''
    return apply_fun_to_file(chr_then_score_sorting, file, dict(), "desc")


def write_for_chr_sorted(tuple_header_and_chr_dict_with_sorted_start):
    header, chr_dict_with_sorted_start = tuple_header_and_chr_dict_with_sorted_start
    print("".join(header))
    for chr in sorted(chr_dict_with_sorted_start.keys()):
        for line in chr_dict_with_sorted_start[chr]:
            print("\t".join(line))


def write_for_size(tuple_header_sorted_list):
    header, sorted_list = tuple_header_sorted_list
    print("".join(header))
    for line in sorted_list:
        print("\t".join(line))


# Intersect and Subtract code

class ReadingError(Exception):
    """raised when reading stops both in A and B files, but the intervals in the files still remain"""
    pass


def bed_reader(file):
    with open(file, "r") as f:
        for line in f:
            if line.startswith("track") or line.startswith("#"):
                continue
            yield line.strip().split()


def sorting_check(bed_file):
    reader = bed_reader(bed_file)
    prev_int = None
    while True:
        try:
            if not prev_int:
                prev_int = next(reader)
            next_int = next(reader)
            if (next_int[0] == prev_int[0] and int(next_int[1]) < int(prev_int[1])) or next_int[0] < prev_int[0]:
                return False
            else:
                prev_int = next_int
        except StopIteration:
            return True


def chr_start_sort(file):
    with open(f"sorted_{file.split('/')[-1][:-4]}.bed", "w") as out_f:
        reader = bed_reader(file)
        intervals = [interval for interval in reader]
        intervals.sort(key=lambda x: (x[0], int(x[1])))
        for interval in intervals:
            out_f.write("\t".join(interval) + "\n")


def overlap_check(a_int, b_int):
    a_start, a_end = int(a_int[1]), int(a_int[2])
    b_start, b_end = int(b_int[1]), int(b_int[2])
    return not (b_end <= a_start or b_start >= a_end)


def get_overlap(a_int, b_int, mode=None, fraction=None):
    a_start, a_end = int(a_int[1]), int(a_int[2])
    b_start, b_end = int(b_int[1]), int(b_int[2])

    overlap = a_int[:]
    if a_start >= b_start:
        overlap[1] = a_start
        if a_end <= b_end:
            overlap[2] = a_end
        else:
            overlap[2] = b_end
    else:
        overlap[1] = b_start
        if b_end >= a_end:
            overlap[2] = a_end
        else:
            overlap[2] = b_end
    if mode == "subtract":
        overlap_len = overlap[2] - overlap[1] - 1
        a_len = a_end - a_start - 1
        return True if overlap_len / a_len >= fraction else False
    return "\t".join([str(i) for i in overlap])


def write_overlaps(a_int, opened_file):
    global b_intervals
    cut_i = 0
    for i, b_int in enumerate(b_intervals):
        if not overlap_check(a_int, b_int):
            if int(a_int[2]) <= int(b_int[1]):
                cut_i = i + 1
        else:
            overlap = get_overlap(a_int, b_int)
            opened_file.write(overlap + "\n")
    b_intervals = b_intervals[cut_i:]


def overlaps_collector(file_a, file_b, mode, **kwargs):
    file_postfix = "intersections" if mode == "intersect" else "subtractions"
    with open(f"{file_a.split('/')[-1][:-4]}_{file_postfix}.bed", "w") as out_f:
        a_parser = bed_reader(file_a)
        b_parser = bed_reader(file_b)
        read_a, read_b = False, True
        a_int = next(a_parser)
        warnings = 0
        try:
            while True:
                if read_a:
                    a_int = next(a_parser)
                    read_a = False
                elif read_b:
                    b_int = next(b_parser)
                    read_b = False
                else:
                    raise ReadingError(f"{a_int, b_int}")

                a_chr, b_chr = a_int[0].replace("chr", ""), b_int[0].replace("chr", "")
                if a_chr != b_chr:
                    if a_chr > b_chr:
                        b_intervals.clear()
                        read_b = True
                    else:
                        if mode == "subtract":
                            write_subtractions(a_int, out_f, non_overlapping=kwargs["A"], fraction=kwargs["f"])
                        elif mode == "intersect":
                            write_overlaps(a_int, out_f)
                        read_a = True
                    continue
                if kwargs["s"] or kwargs["S"]:
                    if len(a_int) >= 6 and len(b_int) >= 6:
                        if kwargs["s"] and a_int[5] != b_int[5]:
                            continue
                        elif kwargs["S"] and a_int[5] != b_int[5]:
                            continue
                    elif warnings == 0:
                        print(
                            "WARNING! Strandedness can't be taken into account since there are no six columns" \
                            "in one of the files")
                        warnings += 1

                if overlap_check(a_int, b_int):  # a_int overlaps b_int
                    b_intervals.append(b_int)
                    read_b = True
                elif int(a_int[1]) >= int(b_int[2]):  # a_start >= b_end
                    read_b = True
                else:  # a_start <= b_end
                    if mode == "subtract":
                        write_subtractions(a_int, out_f, non_overlapping=kwargs["A"], fraction=kwargs["f"])
                    elif mode == "intersect":
                        write_overlaps(a_int, out_f)
                    read_a = True

        except StopIteration:
            try:
                next(b_parser)
            except StopIteration:
                if mode == "subtract":
                    write_subtractions(a_int, out_f, non_overlapping=kwargs["A"], fraction=kwargs["f"])
                    for a_int in a_parser:
                        out_f.write("\t".join(a_int[0:3]) + "\n")
                elif mode == "intersect":
                    if b_intervals:
                        write_overlaps(a_int, out_f)
            return


def write_subtractions(a_int, opened_file, non_overlapping=None, fraction=None):
    global b_intervals
    cut_i = 0
    overlaps = None
    for i, b_int in enumerate(b_intervals):
        if not overlap_check(a_int, b_int):
            if int(a_int[2]) <= int(b_int[1]):
                cut_i = i + 1
        else:
            overlaps = True
            if non_overlapping:
                break
            else:
                if fraction and not get_overlap(a_int, b_int, mode="subtract", fraction=fraction):
                    opened_file.write("\t".join(a_int) + "\n")
                else:
                    subtraction = get_subtraction(a_int, b_intervals[i:], mode="multiple")
                    for interval in subtraction:
                        opened_file.write(interval + "\n")
                b_intervals = b_intervals[cut_i:]
                return
    if not overlaps:
        opened_file.write("\t".join(a_int) + "\n")
    b_intervals = b_intervals[cut_i:]


def get_subtraction(a_int, b_int, mode):
    def subtractor(a_int, b_int):
        a_start, a_end = int(a_int[1]), int(a_int[2])
        b_start, b_end = int(b_int[1]), int(b_int[2])
        left = a_int
        right = a_int
        if a_start >= b_start:
            left = None
            if a_end > b_end:
                right[1] = b_end
                right[2] = a_end
            else:
                right = None
        else:
            if a_end <= b_end:
                right = None
                left[1] = a_start
                left[2] = b_start
            else:
                left[1] = a_start
                left[2] = b_start
                right[1] = b_end
                right[2] = a_end
        return [left, right]

    if mode == "single":
        return subtractor(a_int, b_int)
    elif mode == "multiple":
        temp_1 = [a_int]
        for b in b_int:
            temp_2 = []
            for a in temp_1:
                if a:
                    if not overlap_check(a, b):
                        temp_2.append(a)
                    else:
                        subtraction = subtractor(a, b)
                        for el in subtraction:
                            if el:
                                temp_2.append(el)
            temp_1.clear()
            temp_1 = temp_2
        return ["\t".join([str(i) for i in interval]) for interval in sorted(temp_1, key=lambda x: x[1])]


# Merge code


def mergeable(prev_int, next_int, d=0):
    prev_start, prev_end = int(prev_int[1]), int(prev_int[2])
    next_start, next_end = int(next_int[1]), int(next_int[2])
    return prev_int[0] == next_int[0] and (next_start - prev_end) <= d


def merge(prev_int, next_int):
    #  next start is always >= prev start
    result = prev_int
    if int(next_int[2]) > int(prev_int[2]):
        result[2] = next_int[2]
    return result


def main_merge(file, d=0):
    with open(f"{file.split('/')[-1][:-4]}_merged.bed", "w") as out_f:
        reader = bed_reader(file)
        previous = None
        while True:
            try:
                if not previous:
                    previous = next(reader)
                interval = next(reader)
                if mergeable(previous, interval, d=d):
                    previous = merge(previous, interval)
                else:
                    out_f.write("\t".join(previous) + "\n")
                    previous = interval

            except StopIteration:
                out_f.write("\t".join(previous ) + "\n")
                return


# intersect("wes_71.final_sorted.bed", "bad_exons.bed")
# intersect(argv[1], argv[2])
# subtract("sorted_intervals.bed", "bad_exons.bed", A=None, f=0.1)
# print(sorting_check("test.bed"))
# chr_start_end_sort("tests/unsorted.bed")
# main_merge("wes_71.final_sorted_intersections.bed", d=100)
# overlaps_collector("sorted_intervals.bed", "bad_exons.bed", mode="subtract", s=None, S=None, A=None, f=None)
b_intervals = list()


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

    overlap = a_int
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


def write_subtractions(a_int, opened_file, non_overlapping=None, fraction=None):
    global b_intervals
    cut_i = 0
    for i, b_int in enumerate(b_intervals):
        if not overlap_check(a_int, b_int):
            if int(a_int[2]) <= int(b_int[1]):
                cut_i = i + 1
        else:
            if not non_overlapping:
                if fraction and not get_overlap(a_int, b_int, mode="subtract", fraction=fraction):
                    opened_file.write("\t".join(a_int) + "\n")
                else:
                    subtraction = get_subtraction(a_int, b_intervals[i:], mode="multiple")
                    for interval in subtraction:
                        opened_file.write(interval + "\n")
                b_intervals = b_intervals[cut_i:]
                return
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


def intersect(file_a, file_b, **kwargs):
    with open(f"{file_a.split('/')[-1][:-4]}_intersections.bed", "w") as out_f:
        a_parser = bed_reader(file_a)
        b_parser = bed_reader(file_b)
        read_a, read_b = False, True
        a_int = next(a_parser)
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
                        b_intervals.clear()  # a_int can't intersect with intervals on prev chr
                        read_b = True
                    else:
                        if b_intervals:
                            write_overlaps(a_int, out_f)
                        read_a = True
                    continue

                if kwargs["s"] or kwargs["S"]:
                    warnings = 0
                    if len(a_int) >= 6 and len(b_int) >= 6:
                        if kwargs["s"] and a_int[5] != b_int[5]:
                            continue
                        elif kwargs["S"] and a_int[5] != b_int[5]:
                            continue
                    elif warnings == 0:
                        print(
                            "WARNING! Strandedness can't be taken into account since one of the files does not have 6 columns")
                        warnings += 1

                if overlap_check(a_int, b_int):  # a_int overlaps b_int
                    b_intervals.append(b_int)
                    read_b = True
                elif int(a_int[1]) >= int(b_int[2]):  # a_start >= b_end
                    read_b = True
                else:  # a_start <= b_end
                    if b_intervals:
                        write_overlaps(a_int, out_f)
                    read_a = True

        except StopIteration:
            return


def subtract(file_a, file_b, **kwargs):
    with open(f"{file_a.split('/')[-1][:-4]}_subtractions.bed", "w") as out_f:
        a_parser = bed_reader(file_a)
        b_parser = bed_reader(file_b)
        read_a, read_b = False, True
        a_int = next(a_parser)
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
                        write_subtractions(a_int, out_f, non_overlapping=kwargs["A"], fraction=kwargs["f"])
                        read_a = True
                    continue
                if kwargs["s"] or kwargs["S"]:
                    warnings = 0
                    if len(a_int) >= 6 and len(b_int) >= 6:
                        if kwargs["s"] and a_int[5] != b_int[5]:
                            continue
                        elif kwargs["S"] and a_int[5] != b_int[5]:
                            continue
                    elif warnings == 0:
                        print(
                            "WARNING! Strandedness can't be taken into account since one of the files does not have 6 columns")
                        warnings += 1

                if overlap_check(a_int, b_int):  # a_int overlaps b_int
                    b_intervals.append(b_int)
                    read_b = True
                elif int(a_int[1]) >= int(b_int[2]):  # a_start >= b_end
                    read_b = True
                else:  # a_start <= b_end
                    write_subtractions(a_int, out_f, non_overlapping=kwargs["A"], fraction=kwargs["f"])
                    read_a = True

        except StopIteration:
            try:
                next(b_parser)
            except StopIteration:
                write_subtractions(a_int, out_f, non_overlapping=kwargs["A"], fraction=kwargs["f"])
            for a_int in a_parser:
                out_f.write("\t".join(a_int[0:3]) + "\n")
            return


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
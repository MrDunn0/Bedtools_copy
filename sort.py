import parse_args


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
    with open("sorted_file.txt", "w") as f:
        f.write("".join(header))
        for chr in sorted(chr_dict_with_sorted_start.keys()):
            for line in chr_dict_with_sorted_start[chr]:
                f.write("\t".join(line))
                f.write("\n")


def write_for_size(tuple_header_sorted_list):
    header, sorted_list = tuple_header_sorted_list
    with open("sorted_file.txt", "w") as f:
        f.write("".join(header))
        for line in sorted_list:
            f.write("\t".join(line))
            f.write("\n")


# write_for_chr_sorted(sort_by_default("unsorted.bed"))
# write_for_size(sizeA("unsorted.bed"))
# write_for_size(sizeD("unsorted.bed"))
# write_for_chr_sorted(chrThenSizeA("unsorted.bed"))
# write_for_chr_sorted(chrThenSizeD("unsorted.bed"))
# write_for_chr_sorted(chrThenScoreA("unsorted.bed"))

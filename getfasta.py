import parse_args


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
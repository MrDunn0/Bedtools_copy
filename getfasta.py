class SeqObject():

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
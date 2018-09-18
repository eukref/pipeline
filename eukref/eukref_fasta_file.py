from Bio import SeqIO

from eukref.formatted_header import Format


class EukrefFastaFile:
    def __init__(self, path):
        self.path = path
        self.entries = []

        self.reload()

    def save(self, path):
        SeqIO.write(self.entries, path, "fasta")

    def normalize_headers(self):
        self.entries = [self._normalize_header(e) for e in self.entries]

    def filter_short(self, min_len):
        self.entries = [e for e in self.entries if len(e.seq) >= min_len]

    def reload(self):
        self.entries = []

        with open(self.path, 'r+') as f:
            for record in SeqIO.parse(f, "fasta"):
                self.entries.append(record)

    @staticmethod
    def _normalize_header(entry):
        entry.id = Format.detect_format(entry.id).to_eukref_format()
        entry.description = ''
        entry.name = ''
        return entry

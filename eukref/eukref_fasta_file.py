import os
from Bio import SeqIO
from pathlib import Path
import shutil as sh

from eukref.formatted_header import Format
import eukref.utils as utils


class EukrefFastaFile:
    def __init__(self, path):
        self.path = path
        self.entries = []

        self.reload()

    def save_as(self, path=None):
        """Saves loaded cache to a file"""
        path = path or self.path
        SeqIO.write(self.entries, path, "fasta")
        return EukrefFastaFile(path)

    def remove_source(self):
        os.remove(self.path)

    def normalize_headers(self):
        """Normalize name for each entry in the loaded cache"""
        self.entries = [self._normalize_header(e) for e in self.entries]

    def filter_short(self, min_len):
        """Removes short entries from loaded cache"""
        self.entries = [e for e in self.entries if len(e.seq) >= min_len]

    def reload(self):
        """Reloads the entries from the original file (self.path variable)"""
        self.entries = []

        with open(self.path, 'r+') as f:
            for record in SeqIO.parse(f, "fasta"):
                self.entries.append(record)

    # utils

    def uchime_and_save(self, ref_path, quiet=False):
        """Performs the vsearch unchime and rewrites the original file"""
        unchimed_path = utils.append_to_path(self.path, "_uchimed")

        utils.uchime(self.path, unchimed_path, ref_path, quiet=quiet)

        sh.move(unchimed_path, self.path)
        self.reload()

    def cluster_and_sort(self, quiet=False):
        """Performs vsearch cluster&search and rewrites the original file"""
        cls_srt_path = utils.append_to_path(self.path, "_clustered_sorted")

        utils.cluster_and_sort(self.path, cls_srt_path, quiet=True)

        sh.move(cls_srt_path, self.path)
        self.reload()

    @staticmethod
    def _normalize_header(entry):
        entry.id = Format.detect_format(entry.id).to_eukref_format()
        entry.description = ''
        entry.name = ''
        return entry

import os
import shutil as sh
from pathlib import Path

from eukref.eukref_fasta_file import EukrefFastaFile


class IterationException(Exception):
    pass


class Iteration:
    """A representation of EukRef algorithm iteration
    :param session_path: A sting, the path to session folder
    :param idx: An int, the id of iteration

    """

    def __init__(self, session_path, idx):
        self.name = f"loop_{idx}"
        self.path = Path(session_path).joinpath(self.name)

        self.dataset = None
        self._create_folder()

    def add_dataset(self, path):
        """Copies the dataset specified by path from outside to its folder"""
        return sh.copy(path, self.datasets_path())

    def dataset(self):
        return EukrefFastaFile(self.datasets_path())

    def datasets_path(self):
        return self._iteration_path_for('datasets.fasta')

    def _iteration_path_for(self, *paths):
        return os.path.join(self.path, *paths)

    def _create_folder(self):
        if os.path.exists(self.path):
            raise IterationException(f"The folder {self.path} exists")
        else:
            os.makedirs(self.path)

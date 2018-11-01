import os
import datetime
import shutil as sh
from eukref.iteration import Iteration

from eukref.eukref_fasta_file import EukrefFastaFile


class SessionException(Exception):
    pass


class Session:
    """A representation of EukRef gbretrieve run
        :param session_path: A sting, the path to session folder

        """

    def __init__(self, session_path=None):
        self.name = Session._generate_session_name()
        self.iterations = []

        if not session_path:
            session_path = os.path.join(os.getcwd(), self.name)

        self.path = session_path

        self._create_folder()

    def add_starting_set(self, path):
        """Copies the starting set specified by path from outside to its folder"""
        return sh.copy(path, self.starting_set_path())

    def starting_set(self):
        return EukrefFastaFile(self.starting_set_path())

    def starting_set_path(self):
        return self._session_path_for('starting_set.fasta')

    def _session_path_for(self, *paths):
        return os.path.join(self.path, *paths)

    def _create_folder(self):
        if os.path.exists(self.path):
            raise SessionException(f"The folder {self.path} exists")
        else:
            os.makedirs(self.path)

    # Iteration logic

    def init_next_iteration(self):
        iteration = Iteration(self.path, len(self.iterations) + 1)
        self.iterations.append(iteration)
        return iteration

    @staticmethod
    def _generate_session_name():
        return datetime.datetime.now().strftime("gbretrieve_results_%Y_%m_%d_%H_%M_%S")

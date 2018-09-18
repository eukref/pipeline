import os
import datetime


class SessionException(Exception):
    pass


class Session:
    def __init__(self, session_path=None):
        self.name = Session._generate_session_name()

        if not session_path:
            session_path = os.path.join(os.getcwd(), self.name)

        self.session_path = session_path

        self._create_session_folder()

    def starting_set_path(self):
        return self._session_path_for('starting_set.fasta')

    def _session_path_for(self, *paths):
        return os.path.join(self.session_path, *paths)

    def _create_session_folder(self):
        if os.path.exists(self.session_path):
            raise SessionException(f"The folder {self.session_path} exists")
        else:
            os.makedirs(self.session_path)

    @staticmethod
    def _generate_session_name():
        return datetime.datetime.now().strftime("gbretrieve_results_%Y_%m_%d_%H_%M_%S")

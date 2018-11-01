import unittest
import random
import os
import shutil as sh

from eukref.eukref_fasta_file import *


class TestEukrefFastaFile(unittest.TestCase):
    ACC_1 = "AAD44166.1"
    ACC_2 = "AM711902.1"
    ACC_3 = "ABC123"

    HEADER_1 = f"gi|5524211|gb|{ACC_1}| cytochrome b [Elephas maximus maximus]"
    HEADER_2 = f"{ACC_2} Canis lupus lupus complete mitochondrial genome"
    HEADER_3 = f"{ACC_3} something here"

    TMP_FILE = 'tmp.fasta'
    TMP_FILE_2 = 'tmp2.fasta'

    def setUp(self):
        fasta_content = f">{self.HEADER_1}\n"
        fasta_content += f"{self._generate_random_seq()}\n"
        fasta_content += f">{self.HEADER_2}\n"
        fasta_content += f"{self._generate_random_seq()}\n"
        fasta_content += f">{self.HEADER_3}\n"
        fasta_content += f"{self._generate_random_seq()}\n"

        with open(self.TMP_FILE, 'w') as f:
            f.write(fasta_content)

    def tearDown(self):
        os.remove(self.TMP_FILE)

        if os.path.isfile(self.TMP_FILE_2):
            os.remove(self.TMP_FILE_2)

    def test_entries_load(self):
        f = EukrefFastaFile(self.TMP_FILE)
        self.assertEqual(3, len(f.entries))

    def test_saving_to_file(self):
        f = EukrefFastaFile(self.TMP_FILE)
        f.entries[0].id = 'test'
        f.entries[0].name = ''
        f.entries[0].description = ''
        f.save_as(self.TMP_FILE_2)

        self.assertTrue(os.path.isfile(self.TMP_FILE_2))
        with open(self.TMP_FILE_2) as file:
            self.assertEqual('>test', file.readline().strip())

    def test_headers_normalization(self):
        f = EukrefFastaFile(self.TMP_FILE)
        f.normalize_headers()
        f.save_as(self.TMP_FILE_2)

        headers = [self.ACC_1, self.ACC_2, self.ACC_3]

        with open(self.TMP_FILE_2) as file:
            hdrs = [f for f in file if f[0] == '>']
            for i, line in enumerate(hdrs):
                self.assertEqual(line.strip().replace('>', ''), headers[i])

    def test_filtering(self):
        f = EukrefFastaFile(self.TMP_FILE)
        f.entries[0].seq = 'GTA'

        f.filter_short(10)
        self.assertEqual(len(f.entries), 2)
        self.assertEqual(f.entries[0].description, self.HEADER_2)

    def test_reloading(self):
        f = EukrefFastaFile(self.TMP_FILE)

        with open(self.TMP_FILE, 'a+') as file:
            file.write('>test\n')
            file.write('GTAC\n')

        f.reload()
        self.assertEqual(4, len(f.entries))

    def test_uchime(self):
        ref_path = 'data/reference.fasta'

        f = EukrefFastaFile(self.TMP_FILE)
        f.uchime_and_save(ref_path, quiet=True)

        self.assertEqual(4, len(f.entries))
        self.assertFalse(os.path.isfile('tmp_uchimed.fasta'))

    @staticmethod
    def _generate_random_seq():
        return ''.join(random.choices("GTAC", k=32))


if __name__ == '__main__':
    unittest.main()

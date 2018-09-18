import unittest
from eukref.formatted_header import *


class TestFormat(unittest.TestCase):
    BAR_SEPARATED = 'gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]'
    BAR_SEPARATED_ACC = 'AAD44166.1'

    SIMPLE = 'AM711902.1 Canis lupus lupus complete mitochondrial genome'
    SIMPLE_ACC = 'AM711902.1'

    def test_bar_format_detection(self):
        f = Format.detect_format(self.BAR_SEPARATED)
        self.assertIsInstance(f, BarSeparatedFormat)

    def test_simple_format_detection(self):
        f = Format.detect_format(self.SIMPLE)
        self.assertIsInstance(f, SimpleFormat)

    def test_acc_number_in_bar_separated_format(self):
        f = Format.detect_format(self.BAR_SEPARATED)
        self.assertEqual(self.BAR_SEPARATED_ACC, f.accession_number())

    def test_acc_number_in_simple_format(self):
        f = Format.detect_format(self.SIMPLE)
        self.assertEqual(self.SIMPLE_ACC, f.accession_number())


if __name__ == '__main__':
    unittest.main()

import re


class HeadFormatterException(Exception):
    pass


class Format:
    acc_regexp = None

    @staticmethod
    def detect_format(header):
        klass = None

        for cls in Format.__subclasses__():
            if cls.is_compatible(header):
                klass = cls

        if not klass:
            raise HeadFormatterException(f'Cannot detect format for header "{header}"')

        return klass(header)

    def __init__(self, header):
        self.header = header
        self.acc_number = None

    def accession_number(self):
        if self.acc_number:
            return self.acc_number

        match = re.match(self.acc_regexp, self.header)
        self.acc_number = match.group('acc_number')

        if not self.acc_number:
            raise HeadFormatterException(f'There is no acc number for {self.header}')

        return self.acc_number

    def to_eukref_format(self):
        return f'{self.accession_number()}'

    @classmethod
    def is_compatible(cls, header):
        return re.match(cls.acc_regexp, header)


# gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
class BarSeparatedFormat(Format):
    acc_regexp = re.compile(r'^gi\|.+?\|gb\|(?P<acc_number>.+?)\|.*$')


# AM711902.1 Canis lupus lupus complete mitochondrial genome
class SimpleFormat(Format):
    acc_regexp = re.compile(r'^(?P<acc_number>[^|]+?)(\s.*)?$')

import os
from pathlib import Path
import binascii
import subprocess


def eukref_root(*extra_paths):
    root_path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    return root_path.joinpath(*extra_paths)


def append_to_path(path, suffix):
    folder, f_name = os.path.split(path)
    name, ext = os.path.splitext(f_name)
    f_name = '.'.join([f"{name}_{suffix}", ext])
    return Path(folder).joinpath(f_name)


def get_random_hex():
    return binascii.b2a_hex(os.urandom(15))


def uchime(input_path, output_path, ref_path, quiet=False):
    command = f"vsearch -uchime_ref {ref_path} -db {input_path} -nonchimeras {output_path} -strand plus"

    stdout = open(os.devnull, 'wb') if quiet is True else None
    subprocess.call(command, shell=True, stderr=stdout)

    if stdout:
        stdout.close()


def cluster_and_sort(input_path, output_path, quiet=False):
    tmp_path = append_to_path(input_path, get_random_hex())
    stdout = open(os.devnull, 'wb') if quiet is True else None

    command = f"vsearch --sortbylength {input_path} --output {tmp_path}"
    subprocess.call(command, shell=True, stderr=stdout)

    command = f"vsearch --cluster_smallmem {tmp_path} --id 0.97 --centroids {output_path}"
    subprocess.call(command, shell=True, stderr=stdout)

    os.remove(tmp_path)

    if stdout:
        stdout.close()

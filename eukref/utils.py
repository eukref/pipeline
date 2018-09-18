import os
from pathlib import Path


def eukref_root(*extra_paths):
    root_path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    return root_path.joinpath(*extra_paths)

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for sequence manipulation
"""

import gzip
import os
import shutil
import tempfile
from typing import List

from tqdm import tqdm  # type: ignore

from kmermaid.batch2 import Batch


def set_tempdir(path: str, create: bool = True) -> None:
    """Set new temporary directory and creates it if needed.

    :param path: path to new temporary directory
    :type path: str
    :param create: create if not found, defaults to True
    :type create: bool
    :raises AssertionError: if not found and create is False
    """
    if not os.path.isdir(path):
        if create:
            os.makedirs(path, exist_ok=True)
        else:
            raise AssertionError(f"folder not found: {path}")
    tempfile.tempdir = path


def input_file_exists(path: str) -> None:
    """Check if a file exists.

    :param path: to input file
    :type path: str
    :raises AssertionError: if file not found
    """
    if not os.path.isfile(path):
        raise AssertionError(f"input file not found: {path}")


def copy_batches(
    batches: List[Batch], output_path: str, compress: bool = False
) -> None:
    """Copy generated batches to output folder.

    :param batches: generated batches.
    :type batches: List[Batch]
    :param output_path: path to output folder.
    :type output_path: str
    :param compress: gzip batches, defaults to False
    :type compress: bool
    """
    batch_list = tqdm(
        (batch for batch in batches if os.path.isfile(batch.temp)),
        total=len(batches),
    )
    for current_batch in batch_list:
        if compress:
            with gzip.open(
                os.path.join(output_path, f"{os.path.basename(current_batch.tmp)}.gz"),
                "wb",
            ) as OH, open(current_batch.tmp, "rb") as IH:
                for line in IH:
                    OH.write(line)
        else:
            shutil.copy(current_batch.tmp, output_path)

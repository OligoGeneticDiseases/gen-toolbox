import os
import sys
import hail
from pathlib import Path
import hail.utils

from src.utils.general.string_operations import trim_prefix


def get_metadata(metadata_path):
    p = Path(metadata_path)
    metadata_dict = dict()
    assert p.exists()
    with p.open(encoding="latin-1") as f:
        for line in f.readlines():
            s = line.strip().split("\t")
            ecode = trim_prefix(s[0])
            if ecode not in metadata_dict:
                if len(s) >= 3:
                    metadata_dict[ecode] = [s[1].strip(), s[2].strip()]
                else:
                    metadata_dict[ecode] = ["NA", "NA"]
            else:
                hail.utils.warning(
                    "Metadata: found duplicate key {0} for line {1}. Existing object {2}.".format(
                        ecode, s, (ecode, metadata_dict[ecode])
                    )
                )
    return metadata_dict

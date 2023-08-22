import os

# Import a faster walker if it is installed
try:
    from scandir import walk
except ImportError:
    from os import walk

import shutil
import re

def trim_prefix(filename, seperator="."):
    prefix_clean = filename.rsplit(".")[0].rsplit("_")[0]
    return prefix_clean

def eval_regex(text, regex):
    """
    Will return the parsed string with the defined regex rule
    :param text: Text to be parsed
    :param regex: Regex ruleset
    :return: Parsed string or None if not able to parse.
    """
    result = regex.search(text)
    return result
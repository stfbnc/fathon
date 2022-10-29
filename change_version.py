import os
import sys
import re

INIT_FILE = os.path.join("fathon", "__init__.py")
DOCS_MAIN = os.path.join("docs", "index.rst")
SETUP_FILE = os.path.join(".", "setup.py")


def replace_version(file_name, replace_str):
    with open(file_name, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if re.search("\\d+\\.\\d+\\.\\d+", line) is not None:
            lines[i] = re.sub("\\d+\\.\\d+\\.\\d+", replace_str, line)
            break

    with open(file_name, "w") as f:
        f.write("".join(lines))


if __name__ == "__main__":
    new_version = sys.argv[1]
    replace_version(INIT_FILE, new_version)
    replace_version(DOCS_MAIN, new_version)
    replace_version(SETUP_FILE, new_version)

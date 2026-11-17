#!/usr/bin/env python3
"""Report definite and indirect leaks allocated by the Piernik executable."""

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


LEAK_KINDS = {"Leak_DefinitelyLost", "Leak_IndirectlyLost"}
ALLOCATOR_FUNCTIONS = {
    "malloc", "calloc", "realloc", "reallocarray", "free",
    "strdup", "strndup", "memalign", "posix_memalign", "aligned_alloc",
    "__libc_malloc", "__libc_calloc", "__libc_realloc",
    "operator new", "operator new[]",
}


def allocation_frame(error):
    """Return the first frame outside Valgrind and allocator wrappers."""
    for frame in error.findall("./stack/frame"):
        obj = frame.findtext("obj", default="")
        function = frame.findtext("fn", default="")
        if "vgpreload_" in obj or function in ALLOCATOR_FUNCTIONS:
            continue
        return frame
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("reports", nargs="+", type=Path)
    args = parser.parse_args()

    executable = args.executable.resolve()
    project_leaks = []
    for report in args.reports:
        try:
            root = ET.parse(report).getroot()
        except (OSError, ET.ParseError) as error:
            print(f"{report}: cannot read Valgrind XML: {error}", file=sys.stderr)
            return 1

        statuses = root.findall("status")
        if not statuses or statuses[-1].findtext("state") != "FINISHED":
            print(f"{report}: Valgrind report did not finish", file=sys.stderr)
            return 1

        for error in root.findall("error"):
            kind = error.findtext("kind", default="")
            if kind not in LEAK_KINDS:
                continue
            frame = allocation_frame(error)
            if frame is None:
                continue
            obj = frame.findtext("obj", default="")
            if not obj or Path(obj).resolve() != executable:
                continue

            location = frame.findtext("file", default="unknown source")
            line = frame.findtext("line", default="")
            if line:
                location = f"{location}:{line}"
            description = error.findtext("xwhat/text", default=kind)
            project_leaks.append(f"{report}: {description} ({location})")

    if project_leaks:
        print("\n".join(project_leaks), file=sys.stderr)
        return 1

    print(f"No Piernik-originating definite/indirect leaks in {len(args.reports)} reports.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

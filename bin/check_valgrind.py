#!/usr/bin/env python3
"""Report Valgrind findings and detect leaks allocated by Piernik."""

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


def frame_location(frame):
    if frame is None:
        return "allocation site unavailable"

    function = frame.findtext("fn", default="unknown function")
    source = frame.findtext("file", default="")
    line = frame.findtext("line", default="")
    if source:
        location = f"{source}:{line}" if line else source
        return f"{function} at {location}"
    return function


def print_stack(stack):
    frames = stack.findall("frame")
    if not frames:
        print("    (no frames available)")
    for frame in frames:
        address = frame.findtext("ip", default="?")
        function = frame.findtext("fn", default="unknown function")
        directory = frame.findtext("dir", default="")
        filename = frame.findtext("file", default="")
        line = frame.findtext("line", default="")
        obj = frame.findtext("obj", default="")
        location = str(Path(directory) / filename) if filename else obj
        if line and filename:
            location = f"{location}:{line}"
        print(f"    {address}: {function}" + (f" ({location})" if location else ""))


def print_error_stacks(element):
    for child in element:
        if child.tag == "auxwhat" and child.text:
            print(f"  {child.text}")
        elif child.tag == "stack":
            print_stack(child)
        elif child.tag == "origin":
            description = child.findtext("what")
            if description:
                print(f"  {description}")
            print_error_stacks(child)


def is_project_leak(error, executable):
    frame = allocation_frame(error)
    if frame is None:
        return False, frame

    obj = frame.findtext("obj", default="")
    if not obj:
        return False, frame
    if executable.is_absolute():
        return Path(obj).resolve() == executable, frame
    return Path(obj).name == executable.name, frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-e", "--executable",
        type=Path,
        help="binary path to use instead of the executable recorded in the report",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="show full Valgrind call stacks for selected findings",
    )
    parser.add_argument(
        "-a", "--all",
        action="store_true",
        help="show every unsuppressed Valgrind finding, including MPI/runtime",
    )
    parser.add_argument("reports", nargs="+", type=Path)
    args = parser.parse_args()

    findings = []
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

        executable = args.executable
        if executable is None:
            executable_name = root.findtext("./args/argv/exe", default="")
            if not executable_name:
                print(
                    f"{report}: executable is missing from Valgrind XML; "
                    "pass --executable",
                    file=sys.stderr,
                )
                return 1
            executable = Path(executable_name)
        elif executable.is_absolute():
            executable = executable.resolve()

        for error in root.findall("error"):
            kind = error.findtext("kind", default="")
            frame = None
            project_leak = False
            if kind in LEAK_KINDS:
                project_leak, frame = is_project_leak(error, executable)
            elif args.all:
                stack = error.find("./stack")
                if stack is not None:
                    frame = stack.find("./frame")

            if kind not in LEAK_KINDS and not args.all:
                continue
            if kind in LEAK_KINDS or args.all:
                description = error.findtext("xwhat/text")
                if not description:
                    description = error.findtext("what")
                if not description:
                    description = kind
                findings.append({
                    "report": report,
                    "pid": root.findtext("pid", default="?"),
                    "kind": kind,
                    "description": description,
                    "frame": frame,
                    "project_leak": project_leak,
                    "error": error,
                })

    project_leaks = [finding for finding in findings
                     if finding["project_leak"]]
    display_all = args.all
    displayed = findings if display_all else project_leaks

    if args.verbose:
        for finding in displayed:
            owner = "Piernik" if finding["project_leak"] else "external"
            print(
                f"{finding['report']} (PID {finding['pid']}): "
                f"{owner} {finding['kind']}: {finding['description']}"
            )
            print_error_stacks(finding["error"])
    else:
        output = sys.stdout if display_all else sys.stderr
        for finding in displayed:
            owner = "Piernik" if finding["project_leak"] else "external"
            print(
                f"{finding['report']}: {owner}: {finding['description']} "
                f"({frame_location(finding['frame'])})",
                file=output,
            )

    if project_leaks:
        return 1

    if display_all and not displayed:
        print(f"Displayed {len(findings)} unsuppressed Valgrind findings.")
    elif args.verbose:
        print(
            "No Piernik-originating definite/indirect leaks in "
            f"{len(args.reports)} reports."
        )
    else:
        print(f"No Piernik-originating definite/indirect leaks in {len(args.reports)} reports.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Tests for Valgrind report presentation."""

import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


SCRIPT = Path(__file__).with_name("check_valgrind.py")


class CheckValgrindTest(unittest.TestCase):
    def test_detects_project_leak_from_executable_recorded_in_report(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "report.xml"
            deleted_executable = Path(directory) / "piernik"
            root = ET.Element("valgrindoutput")
            status = ET.SubElement(root, "status")
            ET.SubElement(status, "state").text = "FINISHED"
            args = ET.SubElement(root, "args")
            argv = ET.SubElement(args, "argv")
            ET.SubElement(argv, "exe").text = "./piernik"
            error = ET.SubElement(root, "error")
            ET.SubElement(error, "kind").text = "Leak_DefinitelyLost"
            xwhat = ET.SubElement(error, "xwhat")
            ET.SubElement(xwhat, "text").text = "16 bytes definitely lost"
            stack = ET.SubElement(error, "stack")
            frame = ET.SubElement(stack, "frame")
            ET.SubElement(frame, "fn").text = "test_memory_leak"
            ET.SubElement(frame, "obj").text = str(deleted_executable)
            ET.ElementTree(root).write(report)

            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(report)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn("test_memory_leak", result.stderr)

    def test_verbose_shows_all_stacks_without_changing_leak_status(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "piernik"
            report = Path(directory) / "report.xml"
            root = ET.Element("valgrindoutput")
            ET.SubElement(root, "pid").text = "123"
            status = ET.SubElement(root, "status")
            ET.SubElement(status, "state").text = "FINISHED"
            error = ET.SubElement(root, "error")
            ET.SubElement(error, "kind").text = "Leak_DefinitelyLost"
            xwhat = ET.SubElement(error, "xwhat")
            ET.SubElement(xwhat, "text").text = "8 bytes definitely lost"
            stack = ET.SubElement(error, "stack")
            for name, obj in (("malloc", "/lib/vgpreload_memcheck.so"),
                              ("allocate", str(executable)),
                              ("caller", str(executable))):
                frame = ET.SubElement(stack, "frame")
                ET.SubElement(frame, "ip").text = "0x123"
                ET.SubElement(frame, "fn").text = name
                ET.SubElement(frame, "obj").text = obj
                ET.SubElement(frame, "dir").text = "/source"
                ET.SubElement(frame, "file").text = "sample.F90"
                ET.SubElement(frame, "line").text = "42"
            ET.SubElement(error, "auxwhat").text = "Block was allocated here"
            aux_stack = ET.SubElement(error, "stack")
            aux_frame = ET.SubElement(aux_stack, "frame")
            ET.SubElement(aux_frame, "ip").text = "0x456"
            ET.SubElement(aux_frame, "obj").text = "/lib/unknown.so"
            origin = ET.SubElement(error, "origin")
            ET.SubElement(origin, "what").text = "Created by a heap allocation"
            origin_stack = ET.SubElement(origin, "stack")
            origin_frame = ET.SubElement(origin_stack, "frame")
            ET.SubElement(origin_frame, "fn").text = "origin_function"
            ET.ElementTree(root).write(report)

            result = subprocess.run(
                [sys.executable, str(SCRIPT), "-e", str(executable), "-v", str(report)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 1)
            for text in ("malloc", "allocate", "caller", "/source/sample.F90:42",
                         "Block was allocated here", "0x456", "/lib/unknown.so",
                         "Created by a heap allocation", "origin_function"):
                self.assertIn(text, result.stdout)
            self.assertLess(result.stdout.index("malloc"), result.stdout.index("allocate"))
            self.assertLess(result.stdout.index("allocate"), result.stdout.index("caller"))

    def test_all_shows_external_summary_without_stack(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "report.xml"
            report.write_text(
                "<valgrindoutput><status><state>FINISHED</state></status>"
                "<error><kind>InvalidRead</kind><what>Invalid read</what>"
                "<stack><frame><fn>library_call</fn></frame>"
                "<frame><fn>caller</fn></frame></stack></error></valgrindoutput>"
            )
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "-e", "/piernik", "-a", str(report)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0)
            self.assertIn("Invalid read", result.stdout)
            self.assertIn("library_call", result.stdout)
            self.assertNotIn("caller", result.stdout)

    def test_all_and_verbose_show_external_stacks(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "report.xml"
            report.write_text(
                "<valgrindoutput><status><state>FINISHED</state></status>"
                "<error><kind>InvalidRead</kind><what>Invalid read</what>"
                "<stack><frame><fn>library_call</fn></frame>"
                "<frame><fn>caller</fn></frame></stack></error></valgrindoutput>"
            )
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "-e", "/piernik", "-a", "-v", str(report)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0)
            self.assertIn("Invalid read", result.stdout)
            self.assertIn("library_call", result.stdout)
            self.assertIn("caller", result.stdout)

    def test_color_always_emits_ansi_and_color_never_does_not(self):
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "report.xml"
            report.write_text(
                "<valgrindoutput><status><state>FINISHED</state></status>"
                "<error><kind>InvalidRead</kind><what>Invalid read</what>"
                "</error></valgrindoutput>"
            )
            base = [sys.executable, str(SCRIPT), "-e", "/piernik", "-a"]
            colored = subprocess.run(
                base + ["--color", "always", str(report)],
                capture_output=True, text=True, check=False,
            )
            plain = subprocess.run(
                base + ["--color", "never", str(report)],
                capture_output=True, text=True, check=False,
            )
            self.assertIn("\033[", colored.stdout)
            self.assertNotIn("\033[", plain.stdout)


if __name__ == "__main__":
    unittest.main()

import runpy
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


compare_time = runpy.run_path(str(Path(__file__).with_name("gdf_distance")))["compare_time"]


def time_file(time):
    return SimpleNamespace(attrs={"time": [time]})


class CompareTimeTest(unittest.TestCase):
    def test_opposite_signs_warn_and_return_one(self):
        for t1, t2 in ((-2., 3.), (3., -2.)):
            with self.subTest(t1=t1, t2=t2), self.assertLogs(level="WARNING") as logs:
                self.assertEqual(compare_time(time_file(t1), time_file(t2)), 1.)
            self.assertIn("negative simulation time", logs.output[0])

    def test_other_sign_combinations_do_not_warn(self):
        for t1, t2, expected in ((-2., -3., .2), (-2., 0., 1.), (0., 0., 0.), (2., 3., .2)):
            with self.subTest(t1=t1, t2=t2), patch("logging.warning") as warning:
                self.assertAlmostEqual(compare_time(time_file(t1), time_file(t2)), expected)
                warning.assert_not_called()


if __name__ == "__main__":
    unittest.main()

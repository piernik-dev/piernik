import contextlib
import io
import itertools
import runpy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np


SCRIPT = Path(__file__).with_name("gdf_distance")
module = runpy.run_path(str(SCRIPT))
compare_time = module["compare_time"]


def time_file(time):
    return SimpleNamespace(attrs={"time": [time]})


class GdfTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        print(file=sys.stderr, flush=True)

    def subTest(self, msg=None, **params):
        description = ", ".join(f"{key}={value!r}" for key, value in params.items())
        print(f"  Subtest: {description or msg}", file=sys.stderr, flush=True)
        if msg is None:
            return super().subTest(**params)
        return super().subTest(msg, **params)


class CompareTimeTest(GdfTestCase):
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

    def test_large_times_do_not_overflow(self):
        self.assertAlmostEqual(compare_time(time_file(1.e308), time_file(5.e307)), 1. / 3.)


class CompareDataTest(GdfTestCase):
    def test_small_and_large_values_do_not_underflow_or_overflow(self):
        for value in (1.e-200, 1.e200):
            with self.subTest(value=value):
                result = module["compare_data"](np.array([value]), np.array([2. * value]))
                self.assertAlmostEqual(result, 1. / 3.)

    def test_unrepresentable_difference_is_reported(self):
        with self.assertRaisesRegex(ValueError, "underflowed"):
            module["compare_data"](np.array([1., 0.]), np.array([1., 1.e-300]))


def make_file(path, fields=None, grids=None, particles=None, time=1., edge=2., id_dtype="int64"):
    fields = {"dens": 1.} if fields is None else fields
    grids = [(0, (0, 0, 0), (2, 1, 1))] if grids is None else grids
    with h5py.File(path, "w") as f:
        f.attrs["piernik"] = [2.]
        f.attrs["time"] = [time]
        params = f.create_group("simulation_parameters")
        for key, value in {
            "dimensionality": [1], "geometry": [0], "refine_by": [2],
            "domain_left_edge": [0., 0., 0.], "domain_right_edge": [edge, 1., 1.],
            "domain_dimensions": [2, 1, 1],
        }.items():
            params.attrs[key] = value
        f.create_group("domains/base").attrs["n_d"] = [2, 1, 1]
        f.create_dataset("grid_level", data=[level for level, _, _ in grids])
        f.create_dataset("grid_dimensions", data=[size for _, _, size in grids])
        f.create_dataset("grid_particle_count", data=[len(particles or {})] + [0] * (len(grids) - 1))
        data = f.create_group("data")
        data.attrs["cg_count"] = [len(grids)]
        for i, (level, off, size) in enumerate(grids):
            grid = data.create_group(f"grid_{i:010d}")
            grid.attrs["level"] = [level]
            grid.attrs["off"] = off
            grid.attrs["n_b"] = size
            for name, value in fields.items():
                grid.create_dataset(name, data=np.full(tuple(reversed(size)), value))
            if i == 0 and particles:
                stars = grid.create_group("particles/stars")
                stars.attrs["n_part"] = [len(particles)]
                ids = list(particles)
                stars.create_dataset("id", data=np.array(ids, dtype=id_dtype))
                for key in module["pkeys_flat"]:
                    values = [particles[pid] if key == "position_x" else 0. for pid in ids]
                    stars.create_dataset(key, data=values)


def total_distance(first, second, levels=None):
    with contextlib.redirect_stdout(io.StringIO()):
        norms, by_level, level_range = module["piernik_gdf_compare"](str(first), str(second), levels or [])
    values = list(norms.values()) + [v for level in level_range for v in by_level[level].values()]
    return 1. - np.prod([1. - value for value in values])


class DatasetMetricTest(GdfTestCase):
    def setUp(self):
        super().setUp()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.paths = [Path(self.tmp.name) / f"{i}.h5" for i in range(3)]

    def test_self_comparison_with_extreme_domain_sizes(self):
        for edge in (1.e-200, 1., 1.e200):
            with self.subTest(edge=edge):
                make_file(self.paths[0])
                with h5py.File(self.paths[0], "a") as f:
                    f["simulation_parameters"].attrs["domain_right_edge"] = [edge] * 3
                self.assertEqual(total_distance(self.paths[0], self.paths[0]), 0.)

    def test_domain_overlap_with_extreme_values(self):
        for edge in (1.e-200, 1., 1.e200):
            with self.subTest(edge=edge):
                make_file(self.paths[0])
                make_file(self.paths[1])
                with h5py.File(self.paths[0], "a") as f:
                    f["simulation_parameters"].attrs["domain_right_edge"] = [edge] * 3
                with h5py.File(self.paths[1], "a") as f:
                    f["simulation_parameters"].attrs["domain_right_edge"] = [edge, edge, 2. * edge]
                with h5py.File(self.paths[0], "r") as f1, h5py.File(self.paths[1], "r") as f2:
                    distance, _, _ = module["compare_domains"](f1, f2)
                    self.assertAlmostEqual(distance, .5)

    def test_domain_width_larger_than_float_range(self):
        make_file(self.paths[0])
        make_file(self.paths[1])
        with h5py.File(self.paths[0], "a") as f:
            f["simulation_parameters"].attrs["domain_left_edge"] = [-1.e308, 0., 0.]
            f["simulation_parameters"].attrs["domain_right_edge"] = [1.e308, 1., 1.]
        with h5py.File(self.paths[1], "a") as f:
            f["simulation_parameters"].attrs["domain_right_edge"] = [1.e308, 1., 1.]
        self.assertEqual(total_distance(self.paths[0], self.paths[0]), 0.)
        with h5py.File(self.paths[0], "r") as f1, h5py.File(self.paths[1], "r") as f2:
            distance, _, _ = module["compare_domains"](f1, f2)
            self.assertAlmostEqual(distance, .5)

    def test_particle_id_triangle_counterexample(self):
        make_file(self.paths[0], particles={1: 1., 2: 0., 3: 0.})
        make_file(self.paths[1], particles={2: 0., 3: 0., 4: 0.})
        make_file(self.paths[2], particles={1: -1., 2: 0., 3: 0.})
        ab = total_distance(self.paths[0], self.paths[1])
        bc = total_distance(self.paths[1], self.paths[2])
        ac = total_distance(self.paths[0], self.paths[2])
        self.assertEqual((ab, bc, ac), (1., 1., 1.))
        self.assertLessEqual(ac, ab + bc)

    def test_field_and_geometry_structure_differences(self):
        make_file(self.paths[0])
        make_file(self.paths[1], fields={"other": 1.})
        make_file(self.paths[2], edge=3.)
        for path in self.paths[1:]:
            self.assertEqual(total_distance(self.paths[0], path), 1.)
            self.assertEqual(total_distance(path, self.paths[0]), 1.)
        self.assertEqual(total_distance(self.paths[0], self.paths[0]), 0.)
        output = subprocess.run([sys.executable, str(SCRIPT), str(self.paths[0]), str(self.paths[2])],
                                capture_output=True, text=True, check=True).stdout
        self.assertIn("Total difference between", output)
        self.assertEqual(float(output.splitlines()[-1].split()[-1]), 1.)
        subprocess.run([sys.executable, str(SCRIPT), *map(str, self.paths)],
                       capture_output=True, text=True, check=True)

    def test_refinement_coverage_and_grid_partition(self):
        base = (0, (0, 0, 0), (2, 1, 1))
        make_file(self.paths[0], grids=[base, (1, (0, 0, 0), (2, 2, 2))])
        make_file(self.paths[1], grids=[base, (1, (2, 0, 0), (2, 2, 2))])
        make_file(self.paths[2], grids=[(0, (0, 0, 0), (1, 1, 1)), (0, (1, 0, 0), (1, 1, 1)),
                                        (1, (0, 0, 0), (2, 2, 2))])
        self.assertEqual(total_distance(self.paths[0], self.paths[1]), 1.)
        self.assertEqual(total_distance(self.paths[0], self.paths[2]), 0.)

    def test_matching_structure_and_old_field_name(self):
        make_file(self.paths[0], fields={"cree01": 1.})
        make_file(self.paths[1], fields={"cr_e-e01": 2.})
        self.assertAlmostEqual(total_distance(self.paths[0], self.paths[1]), 1. / 3.)
        self.assertAlmostEqual(total_distance(self.paths[1], self.paths[0]), 1. / 3.)

    def test_duplicate_particle_ids_rejected(self):
        make_file(self.paths[0], particles={1: 1., 2: 0.})
        with h5py.File(self.paths[0], "a") as f:
            ids = f["data/grid_0000000000/particles/stars/id"]
            assert isinstance(ids, h5py.Dataset)
            ids[:] = [1, 1]
        with self.assertRaisesRegex(ValueError, "unique integers"):
            total_distance(self.paths[0], self.paths[0])

    def test_large_particle_ids_remain_distinct(self):
        make_file(self.paths[0], particles={2**53: 0.})
        make_file(self.paths[1], particles={2**53 + 1: 0.})
        self.assertEqual(total_distance(self.paths[0], self.paths[1]), 1.)

    def test_float_encoded_integer_ids(self):
        make_file(self.paths[0], particles={1: 1., 2: 2.}, id_dtype="float64")
        make_file(self.paths[1], particles={2: 2., 1: 1.})
        self.assertEqual(total_distance(self.paths[0], self.paths[1]), 0.)
        with h5py.File(self.paths[0], "a") as f:
            ids = f["data/grid_0000000000/particles/stars/id"]
            assert isinstance(ids, h5py.Dataset)
            ids[0] = 1.5
        with self.assertRaisesRegex(ValueError, "unique integers"):
            total_distance(self.paths[0], self.paths[1])

    def test_triangle_and_symmetry_with_matching_structure(self):
        descriptions = {}
        for i, path in enumerate(self.paths):
            make_file(path, fields={"dens": i + 1., "velx": i + .5, "vely": 2.},
                      particles={1: i + 1., 2: 2.}, time=i + 1.)
            descriptions[path] = (f"{path.name}: time={i + 1}, dens={i + 1}, "
                                  f"velx={i + .5}, particle 1 position_x={i + 1}")
        print("  Shared: vely=2, particle 2 position_x=2, other physical particle values=0",
              file=sys.stderr, flush=True)
        for a, b, c in itertools.permutations(self.paths):
            with self.subTest(a=descriptions[a], b=descriptions[b], c=descriptions[c]):
                ab = total_distance(a, b)
                self.assertAlmostEqual(ab, total_distance(b, a))
                self.assertLessEqual(total_distance(a, c), ab + total_distance(b, c) + 1.e-14)


if __name__ == "__main__":
    unittest.main(verbosity=2)

"""Integration test for the cylindrical discretization at r=0 and at walls."""

import argparse
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def test_uniform_equilibrium(executable):
    with tempfile.TemporaryDirectory(prefix="ugks2d-axis-") as directory:
        work = Path(directory)
        (work / "equilibrium.nml").write_text(
            """&axisym
 nx=8, nr=6, length=1.0, radius=0.6,
 nu=24, ns=10, na=8, cfl=0.45, end_time=0.0008,
 rho_initial=0.001, gas_r=403.5744316182839, uniform_test=.true.,
/
""",
            encoding="ascii",
        )
        run = subprocess.run(
            [str(executable), "--axisymmetric", "equilibrium.nml"],
            cwd=work,
            capture_output=True,
            text=True,
            timeout=60,
        )
        result = work / "axisymmetric_final.dat"
        if run.returncode or not result.is_file():
            raise AssertionError("axisymmetric run did not produce a solution\n" + run.stdout + run.stderr)
        data = np.loadtxt(result, skiprows=1)
        assert data.shape == (48, 8)
        np.testing.assert_allclose(data[:, 2], 0.001, atol=2e-13, rtol=0)
        np.testing.assert_allclose(data[:, 3:5], 0.0, atol=2e-7, rtol=0)
        np.testing.assert_allclose(data[:, 5], 300.0, atol=2e-6, rtol=0)
        np.testing.assert_allclose(data[:, 6], 121.072329485485, atol=2e-7, rtol=0)
        assert np.all(data[:, 7] == 1)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("executable", type=Path)
    options = cli.parse_args()
    test_uniform_equilibrium(options.executable.resolve())
    print("PASS: stationary Maxwellian, cylindrical axis and wall mass balance")

"""Check the five-species axisymmetric entry point in a closed 300 K chamber."""

import argparse
import math
import subprocess
import tempfile
from pathlib import Path


def check_uniform_gas(binary, accuracy=False):
    input_text = """&axisym
 nx=12, nr=6, length=0.6, radius=0.3,
 nu=20, ns=8, na=6, u_scale=1.4, u_center=0, s_max=7,
 gas_r=403.5744316182838, cold_wall=300, aft_wall=300,
 rho_initial=0.0003, cfl=0.45, end_time=0.000035, average_start=0.00001,
 uniform_test=.true., species_file='mixture_species.dat',
/
"""
    if accuracy:
        input_text = input_text.replace(
            "\n/\n",
            """
 x_ref_lo=0.2, x_ref_hi=0.4, x_ref_split=2,
 r_ref_lo=0.0, r_ref_hi=0.1, r_ref_split=2,
 axial_half_range=.true., spatial_order=2, transport_matching=.true.,
 translational_prandtl=0.6666666666666667,
 rotation_tstar=91.5,0.0,100.0,91.5,0.0,
 enable_vibration=.true., vibration_probability=5*0.02,
 vibration_bands(1,1)=3393.5, vibration_bands(1,2)=2294.8,
 vibration_bands(1:2,3)=959.7,959.7,
 vibration_bands(1,4)=3121.9, vibration_bands(1,5)=6332.4,
/
""",
        )
        input_text = input_text.replace("end_time=0.000035, average_start=0.00001", "end_time=0.000006, average_start=0.000002")
    species = """# Synthetic fixture: name molar_mass_g_per_mol mass_fraction diameter_m omega rot_dof rot_probability
N2 28.013 0.35 4.00e-10 0.74 2 0.2
H2O 18.016 0.25 3.00e-10 0.80 3 0.2
CO2 44.011 0.20 4.40e-10 0.85 2 0.2
CO 28.011 0.12 3.80e-10 0.76 2 0.2
H2 2.016 0.08 2.70e-10 0.70 2 0.2
"""
    with tempfile.TemporaryDirectory() as temp:
        directory = Path(temp)
        (directory / "chamber.nml").write_text(input_text, encoding="ascii")
        (directory / "mixture_species.dat").write_text(species, encoding="ascii")
        result = subprocess.run(
            [str(binary), "--axisymmetric", "chamber.nml"],
            cwd=directory, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60, check=False,
        )
        assert result.returncode == 0, result.stdout
        lines = (directory / "axisymmetric_mixture_average.dat").read_text().splitlines()
        assert len(lines) == 1 + (16 * 8 if accuracy else 12 * 6)
        fractions = (0.35, 0.25, 0.20, 0.12, 0.08)
        molar_mass = (28.013, 18.016, 44.011, 28.011, 2.016)
        expected_n = 0.0003 * 6.02214076e26 * sum(y / m for y, m in zip(fractions, molar_mass))
        pressure = expected_n * 1.380649e-23 * 300
        if accuracy:
            bands = ((3393.5,), (2294.8,), (959.7, 959.7), (3121.9,), (6332.4,))
            evib = 0.0003 * 6.02214076e26 * 1.380649e-23 * sum(
                y / m * sum(theta / math.expm1(theta/300) for theta in species_bands)
                for y, m, species_bands in zip(fractions, molar_mass, bands)
            )
        for line in lines[1:]:
            _, _, number, density, ux, ur, trans_t, rot_t, p, *species_columns = map(float, line.split())
            assert len(species_columns) == (7 if accuracy else 6), line
            assert species_columns[-1] == 1, line
            assert max(abs(ux), abs(ur)) < 1e-6, line
            for actual, expected in zip(
                (number, density, trans_t, rot_t, p, *species_columns[:5]),
                (expected_n, 0.0003, 300, 300, pressure, *(0.0003 * y for y in fractions)),
            ):
                assert math.isclose(actual, expected, rel_tol=5e-8), (line, actual, expected)
            if accuracy:
                assert math.isclose(species_columns[-2], evib, rel_tol=5e-8), line
        assert "Completed axisymmetric mixture run." in result.stdout, result.stdout


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("executable", type=Path)
    cli.add_argument("--accuracy", action="store_true", help="Exercise local refinement and all optional closures")
    options = cli.parse_args()
    binary = options.executable.resolve()
    check_uniform_gas(binary)
    if options.accuracy:
        check_uniform_gas(binary, accuracy=True)
        print("PASS axisymmetric nonuniform chamber with all accuracy options enabled")
    print("PASS axisymmetric five-species 300 K equilibrium and wall reflection")

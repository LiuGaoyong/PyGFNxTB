from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional

import numpy as np
from ase import Atoms
from ase import units as U
from ase.calculators import calculator as ase_calc

from pygfnxtb.exe import run_xtb as _run_xtb

__all__ = ["XTB"]


class XTB(ase_calc.Calculator):
    """Base class for XTB calculators."""

    implemented_properties = [
        "energy",
        "forces",
        "stress",
    ]

    default_parameters = {
        "method": "GFN0-xTB",
        "charge": None,
        "multiplicity": None,
        "accuracy": 1.0,
        "guess": "sad",
        "max_iterations": 250,
        "mixer_damping": 0.4,
        "electric_field": None,
        "spin_polarization": None,
        "electronic_temperature": 300.0,
    }

    def __init__(
        self,
        restart=None,
        ignore_bad_restart_file=ase_calc.BaseCalculator._deprecated,
        label=None,
        atoms: Optional[Atoms] = None,
        directory=TemporaryDirectory(),
        **kwargs,
    ) -> None:
        super().__init__(
            atoms=atoms,
            label=label,
            restart=restart,
            ignore_bad_restart_file=ignore_bad_restart_file,
            directory=Path(str(directory)).absolute().__fspath__(),
            **kwargs,
        )
        assert isinstance(self.parameters, ase_calc.Parameters)
        self._method: str = self.parameters.get("method", "GFNFF")
        assert self._method.upper()[:3] == "GFN", f"Invalid: {self._method}."
        if self._method.upper() != "GFNFF":
            assert self._method[3] in "012", f"Invalid: {self._method}."
            self._method = self._method.upper()[:4]
            self.__method = f"--gfn {self._method[3]}"
        else:
            self.__method = "--gfnff"
        assert self._method in ["GFNFF", "GFN0", "GFN1", "GFN2"], (
            f"Invalid: {self._method.upper()}."
        )

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: Optional[list[str]] = None,
        system_changes: list[str] = ase_calc.all_changes,
    ) -> None:
        """Perform actual calculation with by calling the XTB exec."""
        properties = ["energy", "forces"] if not properties else properties
        ase_calc.Calculator.calculate(self, atoms, properties, system_changes)

        # Write coordinates to xyz or POSCAR format.
        assert isinstance(self.atoms, Atoms), "No atoms object set."
        Path(self.directory).mkdir(parents=True, exist_ok=True)
        if np.all(self.atoms.cell.lengths() < 1e-3):
            fname, format = "xtb.xyz", "xyz"
            pbc = False
        else:
            pbc = True
            fname, format = "POSCAR", "vasp"
            if self._method[3] == "2":
                raise ase_calc.CalculatorError(
                    f"The method of {self._method} isn't available"
                    " with periodic boundary conditions."
                )
        self.atoms.write(Path(self.directory) / fname, format=format)

        # Call xtb.exe and chech output files
        args: list[str] = [fname, self.__method, "--grad -I xtb.inp"]
        outs: list[str] = ["energy", "gradient"]
        if pbc:
            outs.append("gradlatt")
        with open(Path(self.directory) / "xtb.inp", "w") as f:
            f.write("$write\n")
            f.write("  dipole=true\n")
            f.write("  charges=true\n")
        content, _, err, is_success, filesexist = _run_xtb(
            *args,
            outputfiles=[f for f in outs],
            workdir=Path(self.directory),
        )
        if not is_success:
            raise ase_calc.CalculatorError(
                f"XTB calculation failed: {err}. "
                f"The XTB command by CLI is {content}."
            )
        for f, exist in zip(outs, filesexist):
            if not exist:
                raise ase_calc.CalculatorError(
                    f"XTB calculation failed: {f} not exist."
                )

        # Parse output files into results
        with Path(self.directory).joinpath("energy").open() as f:
            e = float(f.readlines()[1].split()[1]) * U.Hartree / U.eV
            self.results["energy"] = self.results["free_energy"] = e
        with Path(self.directory).joinpath("gradient").open() as f:
            v = f.readlines()[2 + len(self.atoms) : 2 + 2 * len(self.atoms)]
            f = np.loadtxt(v) * (U.Hartree / U.Bohr) / (U.eV / U.Angstrom)
            assert f.shape == (len(self.atoms), 3), ""
            self.results["forces"] = f
        if pbc:
            with Path(self.directory).joinpath("gradlatt").open() as f:
                virial = np.loadtxt(f.readlines()[2 + 3 : 2 + 2 * 3])
                virial *= (U.Hartree / U.Bohr) / (U.eV / U.Angstrom)
            self.results["stress"] = virial.flat[[0, 4, 8, 5, 2, 1]]


if __name__ == "__main__":
    from ase.build import bulk

    atoms = bulk("Cu", "fcc", 4.3)
    atoms *= (2, 2, 2)
    atoms.calc = XTB()
    print(atoms.get_potential_energy())
    print(atoms.get_forces())
    print(atoms.get_stress())

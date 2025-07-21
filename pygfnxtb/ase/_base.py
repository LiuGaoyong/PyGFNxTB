from pathlib import Path
from tempfile import mkdtemp
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
        "method": "GFN1-xTB",
        "charge": None,
        "multiplicity": None,
        "accuracy": 1.0,
        "guess": "sad",
        "max_iterations": 250,
        "mixer_damping": 0.4,
        "electronic_temperature": 300.0,
    }

    def __init__(
        self,
        restart=None,
        ignore_bad_restart_file=ase_calc.BaseCalculator._deprecated,
        label=None,
        atoms: Optional[Atoms] = None,
        directory: str = mkdtemp(),  # type: ignore
        **kwargs,
    ) -> None:
        super().__init__(
            atoms=atoms,
            label=label,
            restart=restart,
            ignore_bad_restart_file=ignore_bad_restart_file,
            directory=Path(directory).absolute().__fspath__(),
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
            self._method = self._method.upper()
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
        assert isinstance(self.parameters, ase_calc.Parameters)

        # Write coordinates to xyz or POSCAR format.
        assert isinstance(self.atoms, Atoms), "No atoms object set."
        Path(self.directory).mkdir(parents=True, exist_ok=True)
        if np.all(self.atoms.cell.lengths() < 1e-3):
            fname, format = "xtb.xyz", "xyz"
            pbc = False
        else:
            pbc = True
            fname, format = "POSCAR", "vasp"
            if self._method[3] != "0":
                raise ase_calc.CalculatorError(
                    f"The method of {self._method} isn't available with "
                    "periodic boundary conditions. Only GFN0 supports it."
                )
        self.atoms.write(Path(self.directory) / fname, format=format)

        # Call xtb.exe and chech output files
        args: list[str] = [fname, self.__method, "--grad -I xtb.inp"]
        args.append(f"--acc {self.parameters.get('accuracy', 1.0)}")
        args.append("--pop --verbose --norestart")
        with open(Path(self.directory) / "xtb.inp", "w") as f:
            chg: Optional[int] = self.parameters.get("charge", None)
            if chg is not None and int(chg) != 0:
                f.write(f"$chrg {str(int(chg))}")
            mtp: Optional[int] = self.parameters.get("multiplicity", None)
            if mtp is not None:
                f.write(f"$spin {str(int(2 * (mtp + 1) + 1))}")
            f.write("$scc\n")
            guess = self.parameters.get("guess", "sad")
            mixer = self.parameters.get("mixer_damping", 0.4)
            maxiter = self.parameters.get("max_iterations", 250)
            etemp = self.parameters.get("electronic_temperature", 300.0)
            f.write(f"  temp={etemp:.5f}\n")
            f.write(f"  iterations={maxiter:d}\n")
            f.write(f"  broydamp={mixer:.5f}\n")
            f.write(f"  guess={guess:s}\n")
            f.write("$write\n")
            f.write("  esp=true\n")
            f.write("  dipole=true\n")
            f.write("  charges=true\n")
            f.write("  mulliken=true\n")
        outs: list[str] = ["energy", "gradient"]
        if pbc:
            outs.append("gradlatt")
        content, out, err, is_success, filesexist = _run_xtb(
            *args,
            outputfiles=[f for f in outs],
            workdir=Path(self.directory),
        )
        if not is_success:
            raise ase_calc.CalculatorError(
                f"XTB calculation failed:\n{err}\n"
                f"XTB calculation outpt:\n{out}\n"
                f"The XTB command by CLI is {content}\n"
            )
        for f, exist in zip(outs, filesexist):
            if not exist:
                raise ase_calc.CalculatorError(
                    f"XTB calculation failed: {f} not exist."
                )
        # outdata = out.splitlines()
        # TODO: check accuracy of calculation for outdata

        # Parse output files into results
        with Path(self.directory).joinpath("energy").open() as f:
            e = float(f.readlines()[1].split()[1]) * U.Hartree / U.eV
            self.results["energy"] = self.results["free_energy"] = e
        with Path(self.directory).joinpath("gradient").open() as f:
            v = f.readlines()[2 + len(self.atoms) : 2 + 2 * len(self.atoms)]
            f = -np.loadtxt(v) * (U.Hartree / U.Bohr) / (U.eV / U.Angstrom)
            if len(self.atoms) > 1:
                assert f.shape == (len(self.atoms), 3), f"{f.shape}: {f}"
            elif len(self.atoms) == 1:
                assert f.shape == (3,), f"{f.shape}: {f}"
            else:
                pass
            self.results["forces"] = f
        if pbc:
            with Path(self.directory).joinpath("gradlatt").open() as f:
                virial = np.loadtxt(f.readlines()[2 + 3 : 2 + 2 * 3])
                _stress = virial * U.Hartree / self.atoms.get_volume()
            self.results["stress"] = _stress.flat[[0, 4, 8, 5, 2, 1]]


if __name__ == "__main__":
    from ase.build import bulk

    atoms = bulk("Cu", "fcc", 4.3)
    atoms *= (2, 2, 2)
    atoms.calc = XTB(directory=".", method="GFN0-xTB")
    print(atoms.get_potential_energy())
    print(atoms.get_forces())
    print(atoms.get_stress())

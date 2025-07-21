from pathlib import Path
from typing import Optional

import numpy as np
from ase import Atoms
from ase.calculators import calculator as ase_calc

from pygfnxtb.exe import run_xtb as _run_xtb

__all__ = ["XTB"]


class XTB(ase_calc.Calculator):
    """Base class for XTB calculators."""

    implemented_properties = [
        "energy",
        "forces",
        "charges",
        "dipole",
        # "stress", # Cannot support stress calculation for now.
    ]

    default_parameters = {
        "method": "GFN2-xTB",
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
        directory=".",
        **kwargs,
    ):
        super().__init__(
            restart,
            ignore_bad_restart_file,
            label,
            atoms,
            directory,
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

        assert isinstance(self.atoms, Atoms), "No atoms object set."
        Path(self.directory).mkdir(parents=True, exist_ok=True)
        if np.all(self.atoms.cell.lengths() < 1e-3):
            fname, format = "xtb.xyz", "xyz"
        else:
            fname, format = "POSCAR", "vasp"
        self.atoms.write(Path(self.directory) / fname, format=format)

        # xtb [options] <geometry> [options]
        args: list[str] = [fname, self.__method, "--grad"]
        outs: list[str] = ["energy", "gradient", "gradlatt"]
        if self._method.upper() == "GFNFF":
            outs.append("gfnff_charges")
        else:
            outs.append("charges")
        content, out, err, is_success, filesexist = _run_xtb(
            *args,
            outputfiles=[str(i) for i in outs],
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

        print(out)
        assert False, "Not implemented."


if __name__ == "__main__":
    from ase.build import bulk

    atoms = bulk("Cu", "fcc", 4.3)
    atoms *= (2, 2, 2)
    atoms.calc = XTB()
    print(atoms.get_potential_energy())

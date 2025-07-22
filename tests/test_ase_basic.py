from platform import system

import numpy as np
import pytest
from ase.atoms import Atoms
from ase.build import bulk, molecule
from ase.calculators.calculator import InputError
from ase.cluster import Octahedron

from pygfnxtb.ase import XTB


@pytest.mark.parametrize(
    "atoms",
    [
        molecule("H2O"),
        bulk("Cu", "fcc", 3.8, cubic=True) * (2, 2, 2),
    ],
)
@pytest.mark.parametrize("method", ["GFN1-xTB", "gfn0xTB", "GFN2xTB", "gfnff"])
def test_XTB_available(atoms: Atoms, method: str):
    atoms.calc = XTB(method=method)
    if atoms.cell.rank != 0:
        if not method.upper().startswith("GFN0"):
            pytest.skip(f"{method} does not support periodic systems.")
        else:
            print(atoms.get_stress())
    print(atoms.get_potential_energy())
    print(atoms.get_forces())


@pytest.mark.skipif(
    system().lower() == "windows",
    reason="The tblite is not supported on Windows.",
)
@pytest.mark.parametrize("method", ["GFN1-xTB", "GFN2-xTB"])
@pytest.mark.parametrize(
    "atoms",
    [molecule(i) for i in ("CH3CH2OCH3", "C2H6SO", "C6H6", "CH2SCH2")]
    + [Octahedron("Cu", 3, 1)],
)
def test_XTB_by_tblite(atoms: Atoms, method: str):
    from tblite.ase import TBLite

    atoms.calc = XTB(method=method)
    e0 = atoms.get_potential_energy()
    f0 = atoms.get_forces()

    atoms.calc = TBLite(method=method, verbosity=False, cache=False)
    e1 = atoms.get_potential_energy()
    f1 = atoms.get_forces()

    assert np.all(np.abs(e0 - e1) < 0.005)
    assert np.all(np.abs(f0 - f1) < 0.005)


def test_invalid_method():
    """GFN-xTB without method number is invalid, should raise an input error."""
    with pytest.raises(InputError):
        XTB(method="GFN-xTB")

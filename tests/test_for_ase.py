import pytest
from ase.atoms import Atoms
from ase.build import bulk, molecule

from pygfnxtb.ase import XTB


@pytest.mark.parametrize("atoms", [molecule("H2O"), bulk("Cu", "fcc", 3.6)])
def test_XTB_available(atoms: Atoms):
    atoms.calc = XTB()
    atoms.get_potential_energy()

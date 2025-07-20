# ruff: noqa: D100, D103
from pyxtb.exe import XTB_AVAILABLE


def test_available():
    assert XTB_AVAILABLE

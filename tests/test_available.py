# ruff: noqa: D100, D103
from pygfnxtb.exe import XTB_AVAILABLE, XTB_BIN, XTB_EXE


def test_available():
    assert XTB_BIN.exists(), f"XTB_BIN_DIR: {XTB_BIN} does not exist."
    assert XTB_BIN.is_dir(), f"XTB_BIN_DIR: {XTB_BIN} is not a directory."
    assert XTB_EXE.exists(), f"XTB_EXE: {XTB_EXE} does not exist."
    assert XTB_EXE.is_file(), f"XTB_EXE: {XTB_EXE} is not a file."
    assert XTB_AVAILABLE, "..."

"""The XTB Binary Executable Files."""

import subprocess
from os import environ as ENV
from pathlib import Path
from platform import machine, system

__all__ = [
    "XTB_BIN",
    "XTB_EXE",
    "run_xtb",
    "XTB_AVAILABLE",
    "THIS_DIR",
]

if subprocess._mswindows:  # type: ignore
    DEFAULT_EXE = "cmd"
else:
    DEFAULT_EXE = "bash"
THIS_DIR = Path(__file__).parent
XTB_BIN = THIS_DIR.joinpath(f"{system().lower()}_{machine().lower()}")
ENV["XTBPATH"] = THIS_DIR.joinpath("0_params").__fspath__()
XTB_EXE = XTB_BIN.joinpath("xtb.exe")


def run_script(
    content: str,
    exe: str = DEFAULT_EXE,
    outputfiles: list[str | Path] = [],
    workdir: str | Path = Path("."),
    timeout: float | None = None,
) -> tuple[str, str, str, bool, list[bool]]:
    """Run script by CLI based on `subprocess.Popen`.

        https://docs.python.org/3/library/subprocess.html

    Args:
        content (str): The content of script which to be run.
        exe: (str): The executor which to be run.
            If platform is Windows, it defaults to CMD.
            If platform is Mac or Linux, it defaults to Bash.
        outputfiles (list[str | Path], optional): The output files that
            will be check exist of not after script run. Defaults to [].
        workdir (str | Path, optional): The workdir. Defaults to Path(".").
        timeout (float, optional): ...

    Returns:
        tuple[str, bool, list[bool]]:
            The first item is input content before script run.
            The second item is output content after script run.
            The third item is error content after script run.
            The 4th item is whether script run successfully.
            The 5th item is whether exist or not for each output file.
    """
    if not isinstance(workdir, Path):
        workdir = Path(workdir)
    assert isinstance(workdir, Path)
    if workdir.exists():
        assert workdir.is_dir()
    else:
        workdir.mkdir(parents=True)
    kwargs = dict(shell=True, cwd=workdir, env=ENV)
    kwargs["stdin"] = subprocess.PIPE  # type: ignore
    kwargs["stdout"] = subprocess.PIPE  # type: ignore
    kwargs["stderr"] = subprocess.PIPE  # type: ignore
    with subprocess.Popen(exe, **kwargs) as p:  # type: ignore
        content_splitlines = str(content).splitlines()
        for line in content_splitlines:
            bline = bytes(f"{line}\n", encoding="utf-8")
            p.stdin.write(bline)  # type: ignore
        try:
            out, err = p.communicate(timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            p.kill()
            if subprocess._mswindows:  # type: ignore
                # Windows accumulates the output in a single blocking
                # read() call run on child threads, with the timeout
                # being done in a join() on those threads.  communicate()
                # _after_ kill() is required to collect that and add it
                # to the exception.
                exc.stdout, exc.stderr = p.communicate()
            else:
                # POSIX _communicate already populated the output so
                # far into the TimeoutExpired exception.
                p.wait()
            raise
        except:  # Including KeyboardInterrupt, communicate handled that.
            p.kill()
            # We don't call process.wait() as .__exit__ does that for us.
            raise
        is_success = bool(p.poll() == 0)
    filesexist = [workdir.joinpath(f).exists() for f in outputfiles]

    # decode for out & err
    result = dict(out=out, err=err)
    for k, v in result.items():
        assert isinstance(v, (str, bytes)), f"TYPE={type(v)}"
        if isinstance(v, bytes):
            has_decode = False
            encoding_lst = ["UTF-8", "GB2312", "GBK", "ISO-8859-1", "UTF-16"]
            for encoding in encoding_lst:
                try:
                    result[k] = v.decode(encoding=encoding)
                    has_decode = True
                    break
                except Exception:
                    has_decode = False
            if not has_decode:
                raise
    out, err = result["out"], result["err"]
    if is_success and exe == "cmd":
        out, nPS = out.strip().splitlines(), []
        for i, line in enumerate(out):
            if line.lower().startswith(f"{workdir.resolve()}>".lower()):
                nPS.append(i)
        if len(nPS) == 1:
            nPS.append(len(out))
        assert len(nPS) == 2
        out = "\n".join(out[nPS[0] + 1 : nPS[1]])

    return (content, out, err, is_success, filesexist)


def run_xtb(
    *args,
    outputfiles: list[str | Path] = [],
    workdir: str | Path = Path("."),
    timeout: float | None = None,
) -> tuple[str, str, str, bool, list[bool]]:
    """Run XTB by subprocess."""
    arguments = " ".join([str(arg) for arg in args])
    if subprocess._mswindows:  # type: ignore
        content = f'"{XTB_EXE.__fspath__()}" {arguments}'
    else:
        content = f"{XTB_EXE.__fspath__()} {arguments}"
    return run_script(content, DEFAULT_EXE, outputfiles, workdir, timeout)


def __available() -> bool:
    """Check if XTB is available."""
    try:
        assert XTB_BIN.exists(), f"XTB_BIN_DIR: {XTB_BIN} does not exist."
        assert XTB_BIN.is_dir(), f"XTB_BIN_DIR: {XTB_BIN} is not a directory."
        assert XTB_EXE.exists(), f"XTB_EXE: {XTB_EXE} does not exist."
        assert XTB_EXE.is_file(), f"XTB_EXE: {XTB_EXE} is not a file."
        _, out, _, is_success, _ = run_xtb("--version")
        if not is_success:
            return False
        else:
            return "xtb" in out and "Grimme" in out
    except Exception:
        return False


XTB_AVAILABLE = __available()

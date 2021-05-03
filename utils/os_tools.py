import os
import socket
import subprocess
from pathlib import Path
from sys import argv

from setproctitle import setproctitle


def filename_from_argv(argument: str = None) -> Path:
    if len(argv) < 2:
        raise ValueError("specify filename")
    if argument:
        fn = argument
    else:
        fn = argv[1]
    fn = fn.replace(".bin", "").replace(".meta.json", "")
    if fn.endswith("."):
        fn = fn[:-1]
    return Path(fn.replace(".bin", "").replace(".meta.json", ""))


def git_hash() -> str:
    output = subprocess.run(["git", "rev-parse", "--verify", "HEAD"], capture_output=True)
    return output.stdout.decode()


def check_heartbeat_needs_recompile() -> None:
    library = Path("heartbeat/heartbeat.so")
    code = library.with_suffix(".c")
    if code.stat().st_mtime > library.stat().st_mtime:
        raise RuntimeError("heartbeat.so is out of date. Please recompile it from source.")


def process_friendlyness(fn: Path) -> None:
    if socket.gethostname() == "standpc":
        # only handle other computers specially
        return
    setproctitle(f"[{fn.stem}] [rebound-watersim] read /home/winklerl23/sim-info.txt for more information")
    os.nice(5)


def is_ci() -> bool:
    return "CI" in os.environ

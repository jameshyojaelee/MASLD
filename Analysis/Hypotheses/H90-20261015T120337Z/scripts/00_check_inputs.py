"""H90 step 0: verify prespec sha256 and all 24 input hashes; write versions."""
import subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import h90_common as C

C.OUT.mkdir(exist_ok=True)
lines = []
C.check_inputs(lambda s: (print(s, flush=True), lines.append(s)))
(C.OUT / "inputs_check.txt").write_text("\n".join(lines) + "\n")
frz = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout
p = C.RUN / "pip_freeze.txt"
if not p.exists():
    p.write_text(frz)

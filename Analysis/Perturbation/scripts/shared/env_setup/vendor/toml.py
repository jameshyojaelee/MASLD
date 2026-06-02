# Lightweight shim: cell-load v0.8.7 imports `toml` (PyPI), which is not
# present in perturbation_state. Python 3.11's stdlib `tomllib` has the same
# semantics for parse-from-file when opened in binary mode. We keep the shim
# here in the vendored layer so we never have to pip install into
# perturbation_state. (Per project policy: no pip installs in this env.)
import tomllib as _tomllib


def load(fp):
    """Match `toml.load(open(path, 'r'))` semantics by re-opening in binary."""
    name = getattr(fp, "name", None)
    if name is not None:
        with open(name, "rb") as bf:
            return _tomllib.load(bf)
    # Fallback: read text, encode, parse.
    data = fp.read()
    if isinstance(data, str):
        data = data.encode("utf-8")
    return _tomllib.loads(data.decode("utf-8"))


def loads(s):
    return _tomllib.loads(s)

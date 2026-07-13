#!/usr/bin/env python3
"""Run IDR 2.0.4.2 under NumPy >=1.24 without mutating the environment."""

import numpy as np

# IDR 2.0.4.2 uses the removed NumPy scalar aliases internally. Mapping them
# to the behaviorally equivalent builtin types is the upstream-recommended
# compatibility adaptation; inputs, model, ranking, and IDR statistics remain
# unchanged.
if "int" not in np.__dict__:
    np.int = int
if "float" not in np.__dict__:
    np.float = float

import idr.idr

idr.idr.main()


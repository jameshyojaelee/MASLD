"""Python 3.10 compatibility surface for the Python 3.11 ``tomllib`` API.

The cluster's pinned pyBigWig module uses Python 3.10 and already provides
``tomli==2.0.1``.  Re-exporting only the standard-library API keeps the I/O
runtime isolated without adding or installing a dependency.
"""

from tomli import TOMLDecodeError, load, loads

__all__ = ["TOMLDecodeError", "load", "loads"]

"""Climate extreme indices.

Indices are declared in ``conf/mappings/indices.yaml`` rather than written in
code: each entry names the callable to invoke — usually an :mod:`xclim`
indicator — the dataset variables to feed it, and any literal arguments. Adding
an index is therefore a configuration change.

Example:
    >>> from climdata.extremes import extreme_index
    >>> indices = extreme_index(cfg, ds)              # doctest: +SKIP
    >>> tn10p = indices.calculate("tn10p").compute()  # doctest: +SKIP
"""

from .indices import extreme_index

__all__ = ["extreme_index"]

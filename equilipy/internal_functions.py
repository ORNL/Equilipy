"""Legacy internal helper functions."""

from __future__ import annotations

import numpy as np


def _dict2np(dictionary):
    header = list(dictionary.keys())
    vals = list([])
    for i, h in enumerate(header):
        val = dictionary[h]
        try:
            len(val)
        except TypeError:
            val = [val]

        if i == 0:
            vals.append(val)
            column_length = len(val)
        elif len(val) == column_length:
            vals.append(val)
        else:
            print("Error: different length of colums ")
            raise ValueError()
    res = np.squeeze(np.asarray(vals).T)

    try:
        L, n = res.shape
    except ValueError:
        # When only one condition is given
        res = res.reshape((1, len(res)))

    return header, res

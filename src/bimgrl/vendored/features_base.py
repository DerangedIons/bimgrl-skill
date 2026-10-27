# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (block_rush_larsen/features.py), adapted for bimgrl:
# imports made package-relative; no other changes unless noted below.
"""Action potential and calcium transient features of one beat

Shared by ``make_reference.py`` and ``tests/test_block_rush_larsen_accuracy.py``,
so that the reference and the scheme are measured in exactly the same way.
"""

from __future__ import annotations

import numpy as np


def features(t: np.ndarray, V: np.ndarray, Ca: np.ndarray) -> dict[str, float]:
    """Features of one beat, sampled on a (0.1 ms) grid

    Parameters
    ----------
    t : np.ndarray
        The time points
    V : np.ndarray
        The membrane potential at the time points, starting from rest
    Ca : np.ndarray
        The cytosolic calcium at the time points, starting from rest

    Returns
    -------
    dict[str, float]
        ``vpeak``: the maximum of `V`. ``apd90``: the time from the upstroke
        (the first crossing of the midpoint between rest and peak) to the
        first crossing after the peak of ``vpeak - 0.9 * (vpeak - V[0])``.
        ``cat``: the amplitude of the calcium transient, the maximum of `Ca`
        less its starting value. ``tact``: the time of the upstroke, linearly
        interpolated between the grid points.
    """
    i_peak = int(np.nanargmax(V))
    vpeak = V[i_peak]
    vrest = V[0]
    v90 = vpeak - 0.9 * (vpeak - vrest)
    vmid = 0.5 * (vpeak + vrest)
    up = int(np.argmax(V > vmid))
    after = np.where(V[i_peak:] < v90)[0]
    apd90 = t[i_peak + after[0]] - t[up] if len(after) else np.nan
    cat = np.nanmax(Ca) - Ca[0]
    if up > 0:
        tact = t[up - 1] + (vmid - V[up - 1]) * (t[up] - t[up - 1]) / (V[up] - V[up - 1])
    else:
        tact = t[0]
    return dict(vpeak=float(vpeak), apd90=float(apd90), cat=float(cat), tact=float(tact))

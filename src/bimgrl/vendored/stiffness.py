# Vendored from DerangedIons/blocked-rush-larsen @ 58af805 (cell_model_survey/survey/stiffness.py), adapted for bimgrl:
# imports made package-relative; no other changes unless noted below.
"""Stiffness of a model along its reference solution

The exact Jacobian (``gotranx.block_detection.jacobian_function``) is evaluated
every 1 ms along the reference, as Spiteri & Dean (2010) did, and
eigen-decomposed, ``J = V Λ V⁻¹``. The participation of state k in mode j is
``p_kj = V_kj (V⁻¹)_jk`` (Pérez-Arriaga 1982), reported as ``|p|`` normalized
to sum to one over the states, as in the Methods of the draft. A mode with
``Re λ < -20 /ms`` is stiff, and is owned by the states with ``p > 0.1``.

Two predictions of the stability limit of GRL:

* **P1** (the draft's): ``2 |Re λ| / |λ|²`` for the stiffest mode that is owned
  by states updated with forward Euler, i.e. whose diagonal entry of the
  Jacobian vanishes, which gotranx's GRL then steps with forward Euler. A mode
  counts as FE-owned when those states carry at least half of its participation.
* **P2**: the largest ``h`` for which the linearized step
  ``A(h) = I + h diag(φ(h J_ii)) J``, ``φ(z) = (eᶻ - 1)/z``, does not grow faster
  than the linearized model, ``ln ρ(A(h)) / h ≤ max(0, max_j Re λ_j) + 0.1 /ms``,
  at every checkpoint. The rate tolerance lets the step follow the growth of
  the model in the upstroke and in regenerative repolarization, and the slight
  growth of nearly neutral modes, which is not an instability. Beyond the limit
  ``h*`` of a mode, the step grows at about ``(|λ|² h - 2 |Re λ|) / 2``, so the
  tolerance moves the limit to ``h* (1 + 0.1 / |Re λ|)``: by at most 0.5 % for
  a stiff mode (``Re λ < -20 /ms``). P2 also sees the modes of strongly
  coupled states that GRL updates exponentially one by one, which P1 misses.
  The forward Euler analogue uses ``A = I + h J``.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

STIFF = -20.0  # /ms
OWNER = 0.1
FE_DIAGONAL = 1e-8  # the delta of gotranx's GRL, below which it takes a forward Euler step


def participation(J: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The eigenvalues of `J` and the normalized participation ``|p|``, column j for mode j"""
    lam, V = np.linalg.eig(J)
    W = np.linalg.inv(V)
    p = np.abs(V * W.T)
    p /= p.sum(axis=0, keepdims=True)
    return lam, p


def extremes(eigs: Iterable[np.ndarray]) -> dict[str, float]:
    """The extremes of the eigenvalues over the checkpoints, and the percentage of checkpoints
    with a complex pair (as in Table 2 of Spiteri & Dean)"""
    eigs = list(eigs)
    re = [e.real for e in eigs]
    im = [e.imag for e in eigs]
    complex_ = [bool(np.any(np.abs(e.imag) > 1e-9 * max(np.max(np.abs(e)), 1e-300))) for e in eigs]
    return dict(
        min_re=float(min(r.min() for r in re)),
        max_re=float(max(r.max() for r in re)),
        min_im=float(min(i.min() for i in im)),
        max_im=float(max(i.max() for i in im)),
        pct_complex=100.0 * sum(complex_) / len(eigs),
    )


def fe_limit(lam: complex) -> float:
    """The largest h with ``|1 + h λ| ≤ 1``: ``2 |Re λ| / |λ|²``"""
    lam = complex(lam)
    return 2 * abs(lam.real) / abs(lam) ** 2


def stiff_modes(
    J: np.ndarray, names: Sequence[str], fe_states: set[str], threshold: float = STIFF
) -> list[dict]:
    """The modes of `J` with ``Re λ < threshold``, stiffest first, with their owners, the update of
    each owner (``exp`` or ``FE``) and the share of the participation of the `fe_states`"""
    lam, p = participation(J)
    modes = []
    for j in np.argsort(lam.real):
        if lam[j].real >= threshold:
            break
        order = np.argsort(-p[:, j])
        owners = [(names[k], float(p[k, j])) for k in order if p[k, j] > OWNER]
        modes.append(
            dict(
                lam=complex(lam[j]),
                owners=owners,
                updates={name: ("FE" if name in fe_states else "exp") for name, _ in owners},
                fe_share=float(sum(p[k, j] for k in range(len(names)) if names[k] in fe_states)),
                fe_limit=fe_limit(lam[j]),
            )
        )
    return modes


def p1_limit(modes_over_time: Iterable[list[dict]]) -> float:
    """The smallest forward Euler limit of an FE-owned stiff mode, inf if there is none"""
    limits = [m["fe_limit"] for modes in modes_over_time for m in modes if m["fe_share"] >= 0.5]
    return min(limits, default=np.inf)


def _phi(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    small = np.abs(z) < 1e-8
    out = np.empty_like(z)
    out[small] = 1.0 + 0.5 * z[small]
    out[~small] = np.expm1(z[~small]) / z[~small]
    return out


#: How much faster than the linearized model the linearized step may grow, per ms
RATE_TOL = 0.1


def _stable(J: np.ndarray, lam: np.ndarray, h: float, scheme: str) -> bool:
    if scheme == "FE":
        A = np.eye(len(J)) + h * J
    elif scheme == "GRL":
        A = np.eye(len(J)) + h * _phi(h * np.diag(J))[:, None] * J
    else:
        raise ValueError(f"Unknown scheme {scheme!r}")
    rho = np.max(np.abs(np.linalg.eigvals(A)))
    return np.log(max(rho, 1e-300)) / h <= max(0.0, float(lam.real.max())) + RATE_TOL


def growing(lam: np.ndarray) -> bool:
    """Whether the linearized model grows: an eigenvalue with ``Re λ > 1e-6 max |λ|``"""
    return bool(lam.real.max() > 1e-6 * np.max(np.abs(lam)))


def p2_checkpoint(
    J: np.ndarray,
    scheme: str,
    h_lo: float = 1e-6,
    h_hi: float = 10.0,
    n_scan: int = 60,
    rel: float = 1e-3,
) -> float:
    """The largest h at which the linearized step of `scheme` (``"GRL"`` or ``"FE"``) is stable at
    this checkpoint; `h_hi` if it always is, `h_lo` if it never is. A scan on a log grid finds the
    first unstable h, then bisection"""
    lam = np.linalg.eigvals(J)
    hs = np.geomspace(h_lo, h_hi, n_scan)
    bad = next((i for i, h in enumerate(hs) if not _stable(J, lam, h, scheme)), None)
    if bad is None:
        return h_hi
    if bad == 0:
        return h_lo
    lo, hi = hs[bad - 1], hs[bad]
    while hi / lo > 1 + rel:
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if _stable(J, lam, mid, scheme) else (lo, mid)
    return float(lo)


def p2_limit(Js: Iterable[np.ndarray], scheme: str, h_hi: float = 10.0, **kwargs) -> float:
    """The smallest :func:`p2_checkpoint` over the checkpoints, `h_hi` if none limits"""
    return min((p2_checkpoint(J, scheme, h_hi=h_hi, **kwargs) for J in Js), default=h_hi)


def analyse(
    jacobian,
    states: np.ndarray,
    times: np.ndarray,
    parameters: np.ndarray,
    names: Sequence[str],
    symbolic_fe: set[str],
) -> tuple[dict, list[dict]]:
    """The stiffness of a model at the checkpoints

    Parameters
    ----------
    jacobian : Callable
        ``jacobian(t, y, p)``, ordered as the columns of `states`
    states : np.ndarray
        The reference states at the checkpoints, one row per time
    times : np.ndarray
        The checkpoints in ms
    parameters : np.ndarray
        The parameters, with the stimulus of the protocol
    names : Sequence[str]
        The state names, in the order of the columns
    symbolic_fe : set[str]
        The states whose diagonal entry of the Jacobian is zero symbolically

    Returns
    -------
    tuple[dict, list[dict]]
        One row of ``stiffness.csv`` (without the model and protocol) and the rows
        of ``stiffness_modes.csv``, one per checkpoint and stiff mode
    """
    eigs, Js, all_modes, rows, checkpoint_times = [], [], [], [], []
    with np.errstate(all="ignore"):
        for t, y in zip(times, states):
            J = np.asarray(jacobian(t, y, parameters), dtype=float)
            if not np.all(np.isfinite(J)):
                continue
            Js.append(J)
            checkpoint_times.append(float(t))
            fe = set(symbolic_fe) | {
                names[k] for k in range(len(names)) if abs(J[k, k]) < FE_DIAGONAL
            }
            modes = stiff_modes(J, names, fe)
            eigs.append(np.linalg.eigvals(J))
            all_modes.append(modes)
            for m in modes:
                rows.append(
                    dict(
                        t=float(t),
                        lam_re=m["lam"].real,
                        lam_im=m["lam"].imag,
                        owners=";".join(f"{n}:{p:.3f}" for n, p in m["owners"]),
                        updates=";".join(f"{n}:{u}" for n, u in m["updates"].items()),
                        fe_share=m["fe_share"],
                        fe_limit_ms=m["fe_limit"],
                    )
                )
    summary = dict(n_checkpoints=len(Js), n_skipped=len(times) - len(Js))
    summary.update(extremes(eigs))
    stiffest = min(rows, key=lambda r: r["lam_re"], default=None)
    fe_owned = [r for r in rows if r["fe_share"] >= 0.5]
    stiffest_fe = min(fe_owned, key=lambda r: r["fe_limit_ms"], default=None)
    summary.update(
        stiffest_lam=stiffest["lam_re"] if stiffest else np.nan,
        stiffest_t=stiffest["t"] if stiffest else np.nan,
        stiffest_owners=stiffest["owners"] if stiffest else "",
        stiffest_updates=stiffest["updates"] if stiffest else "",
        fe_owned_lam=stiffest_fe["lam_re"] if stiffest_fe else np.nan,
        fe_owned_owners=stiffest_fe["owners"] if stiffest_fe else "",
        p1_grl_ms=p1_limit(all_modes),
        n_growing=sum(growing(e) for e in eigs),
    )
    for scheme in ("GRL", "FE"):
        limits = [p2_checkpoint(J, scheme) for J in Js]
        h, t = min(zip(limits, checkpoint_times), default=(10.0, np.nan))
        summary[f"p2_{scheme.lower()}_ms"], summary[f"p2_{scheme.lower()}_t"] = h, t
    summary.update(
        fe_limit_ms=fe_limit(complex(summary["min_re"])) if eigs else np.nan,
    )
    return summary, rows


def summarize_modes(rows: Iterable[dict]) -> list[dict]:
    """One row per model, protocol, beats and set of owners of the stiff modes: the most negative
    ``Re λ``, its time, imaginary part and updates, the number of checkpoints, the largest FE share
    and the smallest forward Euler limit. The tables use this; the full per-checkpoint table is
    kept compressed and is not committed"""
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        owners = ";".join(sorted(o.split(":")[0] for o in str(r["owners"]).split(";") if o))
        groups.setdefault((r["model"], r["protocol"], str(r["beats"]), owners), []).append(r)
    out = []
    for (model, protocol, beats, owners), group in sorted(groups.items()):
        stiffest = min(group, key=lambda r: float(r["lam_re"]))
        out.append(
            dict(
                model=model,
                protocol=protocol,
                beats=beats,
                owners=owners,
                n_checkpoints=len(group),
                min_lam_re=float(stiffest["lam_re"]),
                lam_im_at_min=float(stiffest["lam_im"]),
                t_at_min=float(stiffest["t"]),
                participation_at_min=stiffest["owners"],
                updates_at_min=stiffest["updates"],
                max_fe_share=max(float(r["fe_share"]) for r in group),
                min_fe_limit_ms=min(float(r["fe_limit_ms"]) for r in group),
            )
        )
    return out

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

ArrayF = NDArray[np.float64]
ArrayC = NDArray[np.complex128]


@dataclass(slots=True)
class ChebyshevGrid:
    """Chebyshev–Gauss–Lobatto 网格，映射到 z in [0, zc]。

    节点顺序：
    - 第 0 个点对应 z = zc（匹配点）
    - 最后一个点对应 z = 0（无穷远）
    """

    n: int
    zc: float
    x: ArrayF
    z: ArrayF
    D: ArrayC


def chebyshev_lobatto_grid(n: int, zc: float) -> ChebyshevGrid:
    if n < 2:
        raise ValueError(f"n 必须至少为 2，当前 n={n}")
    if zc <= 0.0:
        raise ValueError(f"zc 必须为正，当前 zc={zc}")

    k = np.arange(n, dtype=float)
    x = np.cos(np.pi * k / (n - 1))

    c = np.ones(n, dtype=float)
    c[0] = 2.0
    c[-1] = 2.0
    c = c * ((-1.0) ** k)

    # 关键修正：必须是 x_i - x_j
    dX = x[:, None] - x[None, :]

    D_x = np.outer(c, 1.0 / c) / (dX + np.eye(n))
    D_x = D_x - np.diag(np.sum(D_x, axis=1))

    z = 0.5 * zc * (1.0 + x)
    D_z = (2.0 / zc) * D_x

    return ChebyshevGrid(
        n=n,
        zc=zc,
        x=x,
        z=z,
        D=np.asarray(D_z, dtype=complex),
    )

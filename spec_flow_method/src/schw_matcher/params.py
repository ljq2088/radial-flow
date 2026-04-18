from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(slots=True, frozen=True)
class SolverParams:
    """参数容器。

    约定：
    - z = 1/r
    - 视界位置 z_h = 1/r0
    - 匹配点 z_match 需满足 0 < z_match < z_h
    """

    r0: float = 2.0
    omega: float = 0.2
    ell: int = 0
    z_match: float = 0.25
    n_cheb: int = 25
    flow_eps: float = 1.0e-6
    flow_rtol: float = 1.0e-10
    flow_atol: float = 1.0e-12

    def __post_init__(self) -> None:
        if self.r0 <= 0.0:
            raise ValueError(f"r0 必须为正，当前 r0={self.r0}")
        if self.omega == 0.0:
            raise ValueError("当前实现要求 omega != 0")
        if self.ell < 0:
            raise ValueError(f"ell 必须为非负整数，当前 ell={self.ell}")
        if self.n_cheb < 4:
            raise ValueError(f"n_cheb 至少取 4，当前 n_cheb={self.n_cheb}")
        if self.flow_eps <= 0.0:
            raise ValueError(f"flow_eps 必须为正，当前 flow_eps={self.flow_eps}")
        if not (0.0 < self.z_match < self.z_horizon):
            raise ValueError(
                f"需要满足 0 < z_match < z_h=1/r0={self.z_horizon}, 当前 z_match={self.z_match}"
            )
        if self.flow_eps >= (self.z_horizon - self.z_match):
            raise ValueError(
                "flow_eps 过大：从 z_h - flow_eps 出发后会越过匹配点，"
                f"当前 z_h-z_match={self.z_horizon - self.z_match}, flow_eps={self.flow_eps}"
            )

    @property
    def z_horizon(self) -> float:
        return 1.0 / self.r0

    @property
    def r_match(self) -> float:
        return 1.0 / self.z_match

    @property
    def z_start(self) -> float:
        return self.z_horizon - self.flow_eps

    def with_updates(self, **kwargs: object) -> "SolverParams":
        return replace(self, **kwargs)

    def to_dict(self) -> dict[str, float | int]:
        return {
            "r0": self.r0,
            "omega": self.omega,
            "ell": self.ell,
            "z_match": self.z_match,
            "n_cheb": self.n_cheb,
            "flow_eps": self.flow_eps,
            "flow_rtol": self.flow_rtol,
            "flow_atol": self.flow_atol,
        }

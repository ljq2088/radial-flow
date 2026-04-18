from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(slots=True, frozen=True)
class SolverParams:
    """Kerr s=-2 Teukolsky 径向 Riccati + 外侧谱拼接参数。

    约定：
    - 内侧 flow 用 r 变量，从 r_+ + flow_eps 积到 r_match。
    - 外侧谱方法用 z = 1/r，区间为 z \in [0, z_match]。
    - 若未显式提供 lambda_sep，则默认使用 Schwarzschild 极限近似
      lambda = ell (ell + 1) - 2.
    """

    M: float = 1.0
    a: float = 0.3
    omega: float = 0.2
    ell: int = 2
    m_mode: int = 2
    lambda_sep: complex | None = None

    r_match: float = 8.0
    n_cheb: int = 32
    flow_eps: float = 1.0e-6
    flow_rtol: float = 1.0e-10
    flow_atol: float = 1.0e-12

    def __post_init__(self) -> None:
        if self.M <= 0.0:
            raise ValueError(f"M 必须为正，当前 M={self.M}")
        if abs(self.a) >= self.M:
            raise ValueError(f"当前只支持非极端 Kerr，要求 |a| < M，当前 a={self.a}, M={self.M}")
        if self.omega == 0.0:
            raise ValueError("当前实现要求 omega != 0")
        if self.ell < 2:
            raise ValueError(f"对于 s=-2，建议 ell >= 2，当前 ell={self.ell}")
        if abs(self.m_mode) > self.ell:
            raise ValueError(f"需满足 |m_mode| <= ell，当前 m_mode={self.m_mode}, ell={self.ell}")
        if self.n_cheb < 8:
            raise ValueError(f"n_cheb 至少取 8，当前 n_cheb={self.n_cheb}")
        if self.flow_eps <= 0.0:
            raise ValueError(f"flow_eps 必须为正，当前 flow_eps={self.flow_eps}")
        if self.r_match <= self.r_plus + self.flow_eps:
            raise ValueError(
                f"需要满足 r_match > r_plus + flow_eps，当前 r_match={self.r_match}, "
                f"r_plus={self.r_plus}, flow_eps={self.flow_eps}"
            )
        if abs(self.k_horizon) == 0.0:
            raise ValueError("当前 Riccati 变量要求 k_H = omega - m Omega_H != 0")

    @property
    def r_plus(self) -> float:
        return self.M + (self.M**2 - self.a**2) ** 0.5

    @property
    def r_minus(self) -> float:
        return self.M - (self.M**2 - self.a**2) ** 0.5

    @property
    def omega_horizon(self) -> float:
        return self.a / (self.r_plus**2 + self.a**2)

    @property
    def k_horizon(self) -> float:
        return self.omega - self.m_mode * self.omega_horizon

    @property
    def z_match(self) -> float:
        return 1.0 / self.r_match

    @property
    def r_start(self) -> float:
        return self.r_plus + self.flow_eps

    @property
    def lambda_value(self) -> complex:
        if self.lambda_sep is not None:
            return complex(self.lambda_sep)
        return complex(self.ell * (self.ell + 1) - 2.0)

    def with_updates(self, **kwargs: object) -> "SolverParams":
        return replace(self, **kwargs)

    def to_dict(self) -> dict[str, float | int | complex | None]:
        return {
            "M": self.M,
            "a": self.a,
            "omega": self.omega,
            "ell": self.ell,
            "m_mode": self.m_mode,
            "lambda_sep": self.lambda_sep,
            "r_match": self.r_match,
            "n_cheb": self.n_cheb,
            "flow_eps": self.flow_eps,
            "flow_rtol": self.flow_rtol,
            "flow_atol": self.flow_atol,
        }

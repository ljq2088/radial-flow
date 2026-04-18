from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wlexpr


REPO_ROOT = Path(__file__).resolve().parents[3]
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
PINN_ROOT = REPO_ROOT.parent / "PINN" / "SolvingTeukolsky"

for path in (SRC_ROOT, REPO_ROOT, PINN_ROOT):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from domain.safe_eval import try_compute_lambda_status, try_compute_ramp_status


@dataclass
class ScanSettings:
    l: int
    m: int
    s: int
    a_min: float
    a_max: float
    n_a: int
    omega_min: float
    n_omega: int
    omega_scale: float
    eta_max: float
    k_margin: float
    r_match: float
    n_cheb: int
    mma_ref_mode: str
    mma_kernel_path: str
    mma_wl_path_win: str
    mma_timeout_sec: float


def update_progress(prefix: str, done: int, total: int) -> None:
    width = 32
    frac = 0.0 if total <= 0 else done / total
    n_fill = min(width, max(0, int(round(width * frac))))
    bar = "#" * n_fill + "-" * (width - n_fill)
    sys.stdout.write(f"\r{prefix} [{bar}] {done}/{total} ({100.0 * frac:5.1f}%)")
    sys.stdout.flush()
    if done >= total:
        sys.stdout.write("\n")
        sys.stdout.flush()


def _wl_num(x: float) -> str:
    x = float(x)
    if x == 0.0:
        return "0."
    s = f"{x:.17e}"
    mant, exp = s.split("e")
    exp_i = int(exp)
    mant = mant.rstrip("0").rstrip(".")
    if "." not in mant:
        mant = mant + "."
    return f"{mant}*^{exp_i}"


def _parse_rin_result(result) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(result, dtype=object)
    if arr.ndim != 2 or arr.shape[1] < 3:
        raise ValueError(f"Unexpected Mathematica output shape: {arr.shape}")

    r_eval = np.asarray(arr[:, 0], dtype=float)
    re = np.asarray(arr[:, 1], dtype=float)
    im = np.asarray(arr[:, 2], dtype=float)
    return r_eval, re + 1j * im


class MathematicaRinSampler:
    def __init__(self, kernel_path: str, wl_path_win: str, timeout_sec: float = 8.0):
        self.kernel_path = str(kernel_path)
        self.wl_path_win = str(wl_path_win)
        self.timeout_sec = float(timeout_sec)
        self._session: WolframLanguageSession | None = None

    def _get_session(self) -> WolframLanguageSession:
        if self._session is None:
            self._session = WolframLanguageSession(kernel=self.kernel_path)
            self._session.evaluate(wlexpr(rf'Get["{self.wl_path_win}"]'))
        return self._session

    def close(self) -> None:
        if self._session is not None:
            self._session.terminate()
            self._session = None

    def evaluate_rin_at_points_direct(
        self,
        s: int,
        l: int,
        m: int,
        a: float,
        omega: float,
        r_query: np.ndarray,
    ) -> np.ndarray:
        r_query = np.asarray(r_query, dtype=float).reshape(-1)
        if r_query.size == 0:
            raise ValueError("r_query is empty")

        rp = r_plus_of_a(float(a))
        if np.any(r_query <= rp):
            bad = r_query[r_query <= rp][:5]
            raise ValueError(f"Some query radii are not outside the horizon. r_+={rp}, bad={bad}")

        r_list = ", ".join(_wl_num(float(r)) for r in r_query)
        expr = (
            "Quiet[Check[TimeConstrained["
            f"SampleRinAtPoints[{int(s)}, {int(l)}, {int(m)}, "
            f"{_wl_num(a)}, {_wl_num(omega)}, "
            f"{{{r_list}}}], {_wl_num(self.timeout_sec)}, $Aborted], $Failed]]"
        )
        result = self._get_session().evaluate(wlexpr(expr))
        if str(result) == "$Failed":
            raise RuntimeError("Mathematica evaluation returned $Failed")
        if str(result) == "$Aborted":
            raise TimeoutError(f"Mathematica evaluation timed out after {self.timeout_sec} s")
        r_eval, rin = _parse_rin_result(result)

        if len(r_eval) != len(r_query):
            raise ValueError(
                f"Length mismatch: Mathematica returned {len(r_eval)} points, requested {len(r_query)}."
            )
        max_diff = float(np.max(np.abs(r_eval - r_query)))
        if max_diff > 1.0e-9:
            raise ValueError(f"r_points mismatch between Python and Mathematica, max diff = {max_diff:.3e}")
        return rin


def r_plus_of_a(a: float, M: float = 1.0) -> float:
    return M + math.sqrt(max(M * M - a * a, 0.0))


def r_minus_of_a(a: float, M: float = 1.0) -> float:
    return M - math.sqrt(max(M * M - a * a, 0.0))


def omega_horizon(a: float, M: float = 1.0) -> float:
    rp = r_plus_of_a(a, M=M)
    return a / (2.0 * M * rp)


def k_horizon(a: float, omega: float, m: int, M: float = 1.0) -> float:
    return omega - m * omega_horizon(a, M=M)


def build_a_grid(a_min: float, a_max: float, n_a: int) -> np.ndarray:
    eps_min = 1.0 - a_max
    eps_max = 1.0 - a_min
    eps_grid = np.geomspace(eps_min, eps_max, n_a)
    a_grid = 1.0 - eps_grid[::-1]
    a_grid[0] = a_min
    a_grid[-1] = a_max
    return a_grid


def eta_from_omega(omega: np.ndarray | float, omega_scale: float) -> np.ndarray | float:
    return (2.0 / math.pi) * np.arctan(np.asarray(omega) / omega_scale)


def omega_from_eta(eta: np.ndarray | float, omega_scale: float) -> np.ndarray | float:
    return omega_scale * np.tan(0.5 * math.pi * np.asarray(eta))


def build_eta_grid(omega_min: float, n_omega: int, omega_scale: float, eta_max: float) -> np.ndarray:
    eta_min = float(eta_from_omega(omega_min, omega_scale))
    return np.linspace(eta_min, eta_max, n_omega)


def centers_to_edges(vals: np.ndarray) -> np.ndarray:
    vals = np.asarray(vals, dtype=float)
    if vals.ndim != 1 or len(vals) < 2:
        raise ValueError("Need at least two grid centers to build edges.")
    mids = 0.5 * (vals[:-1] + vals[1:])
    left = vals[0] - 0.5 * (vals[1] - vals[0])
    right = vals[-1] + 0.5 * (vals[-1] - vals[-2])
    edges = np.concatenate([[left], mids, [right]])
    return edges


def build_reference_radii(a: float, mode: str) -> np.ndarray:
    rp = r_plus_of_a(a)
    rm = r_minus_of_a(a)
    gap = max(rp - rm, 1.0e-8)

    if mode == "relative":
        raw = [
            rp + max(1.0e-4, 1.0e-2 * gap),
            rp + max(5.0e-4, 5.0e-2 * gap),
            1.05 * rp,
            1.5 * rp,
            max(8.0, 2.0 * rp),
            max(20.0, 5.0 * rp),
        ]
    elif mode == "fixed":
        raw = [
            rp + max(1.0e-4, 1.0e-2 * gap),
            max(2.2, 1.05 * rp),
            4.0,
            8.0,
            20.0,
            50.0,
        ]
    else:
        raise ValueError(f"Unknown mma_ref_mode={mode}")

    refs = sorted({float(r) for r in raw if r > rp})
    return np.asarray(refs, dtype=float)


def classify_lambda(a: float, omega: float, l: int, m: int, s: int) -> tuple[bool, str]:
    st = try_compute_lambda_status(a=a, omega=omega, l=l, m=m, s=s)
    return st.ok, st.code


def classify_ramp(
    a: float,
    omega: float,
    l: int,
    m: int,
    s: int,
    r_match: float,
    n_cheb: int,
) -> tuple[bool, str]:
    lam_status = try_compute_lambda_status(a=a, omega=omega, l=l, m=m, s=s)
    if not lam_status.ok:
        return False, lam_status.code
    st = try_compute_ramp_status(
        a=a,
        omega=omega,
        l=l,
        m=m,
        s=s,
        lambda_sep=lam_status.value,
        r_match=r_match,
        n_cheb=n_cheb,
    )
    return st.ok, st.code


def classify_mma(
    sampler,
    a: float,
    omega: float,
    l: int,
    m: int,
    s: int,
    mma_ref_mode: str,
) -> tuple[bool, str]:
    try:
        if abs(k_horizon(a, omega, m)) < 1.0e-2:
            return False, "mma_k_horizon_margin"
        r_refs = build_reference_radii(a, mode=mma_ref_mode)
        vals = sampler.evaluate_rin_at_points_direct(
            s=s,
            l=l,
            m=m,
            a=a,
            omega=omega,
            r_query=r_refs,
        )
        vals = np.asarray(vals, dtype=np.complex128)
        if not np.all(np.isfinite(vals.real) & np.isfinite(vals.imag)):
            return False, "mma_nonfinite"
        return True, "ok"
    except Exception as exc:
        return False, f"mma:{type(exc).__name__}"


def scan_mode(
    mode_name: str,
    classify_fn: Callable[[float, float], tuple[bool, str]],
    a_vals: np.ndarray,
    omega_vals: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, Counter, dict[str, list[dict[str, float | str]]]]:
    safe = np.zeros((len(a_vals), len(omega_vals)), dtype=bool)
    codes = np.empty((len(a_vals), len(omega_vals)), dtype=object)
    counts: Counter = Counter()
    failure_examples: dict[str, list[dict[str, float | str]]] = {}
    total = len(a_vals) * len(omega_vals)
    done = 0

    for j, omega in enumerate(omega_vals):
        for i, a in enumerate(a_vals):
            ok, code = classify_fn(float(a), float(omega))
            safe[i, j] = ok
            codes[i, j] = code
            counts[code] += 1
            if not ok:
                lst = failure_examples.setdefault(str(code), [])
                if len(lst) < 12:
                    lst.append({"a": float(a), "omega": float(omega), "code": str(code)})
            done += 1
            update_progress(f"[{mode_name}]", done, total)

    return safe, codes, counts, failure_examples


def setup_omega_ticks(ax: plt.Axes, omega_scale: float) -> None:
    tick_omegas = np.array([1.0e-6, 1.0e-4, 1.0e-2, 1.0, 10.0, 100.0], dtype=float)
    tick_etas = eta_from_omega(tick_omegas, omega_scale)
    mask = (tick_etas >= ax.get_ylim()[0]) & (tick_etas <= ax.get_ylim()[1])
    ax.set_yticks(tick_etas[mask])
    labels = []
    for w in tick_omegas[mask]:
        if w < 1.0e-3:
            labels.append(f"{w:.0e}")
        elif w < 1.0:
            labels.append(f"{w:.2g}")
        elif w < 10.0:
            labels.append(f"{w:.0f}")
        else:
            labels.append(f"{w:.0f}")
    ax.set_yticklabels(labels)


def plot_triptych(
    out_png: Path,
    settings: ScanSettings,
    a_vals: np.ndarray,
    eta_vals: np.ndarray,
    mode_data: list[tuple[str, np.ndarray, Counter]],
) -> None:
    a_edges = centers_to_edges(a_vals)
    eta_edges = centers_to_edges(eta_vals)

    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(18, 6.5),
        sharey=True,
        constrained_layout=True,
    )

    cmap = ListedColormap(["#d9d9d9", "#1b9e77"])

    k0 = np.array([float(settings.m) * omega_horizon(float(a)) for a in a_vals], dtype=float)
    k0_eta = eta_from_omega(k0, settings.omega_scale)
    km_minus = np.maximum(k0 - settings.k_margin, settings.omega_min)
    km_plus = k0 + settings.k_margin
    km_minus_eta = eta_from_omega(km_minus, settings.omega_scale)
    km_plus_eta = eta_from_omega(km_plus, settings.omega_scale)

    titles = {
        "lambda": r"$\lambda(a,\omega)$ finite",
        "ramp": r"$R_{\rm amp}$ finite",
        "mma": r"MMA $R_{\rm in}(r_{\rm ref})$ finite",
    }

    for ax, (name, safe_mask, counts) in zip(axes, mode_data):
        mesh = ax.pcolormesh(
            a_edges,
            eta_edges,
            safe_mask.T.astype(int),
            cmap=cmap,
            shading="auto",
            vmin=0,
            vmax=1,
        )
        ax.plot(a_vals, k0_eta, color="#2c3e50", lw=1.4, label=r"$k_H=0$")
        ax.plot(a_vals, km_minus_eta, color="#2c3e50", lw=0.9, ls="--", alpha=0.7, label=rf"$|k_H|={settings.k_margin:g}$")
        ax.plot(a_vals, km_plus_eta, color="#2c3e50", lw=0.9, ls="--", alpha=0.7)
        ax.set_title(titles[name], fontsize=12)
        ax.set_xlabel("a")
        ax.set_xlim(a_vals.min(), a_vals.max())
        ax.set_ylim(eta_vals.min(), eta_vals.max())
        setup_omega_ticks(ax, settings.omega_scale)
        safe_n = int(np.count_nonzero(safe_mask))
        total_n = safe_mask.size
        ax.text(
            0.02,
            0.02,
            f"safe={safe_n}/{total_n}\nfail={total_n - safe_n}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
        )

    axes[0].set_ylabel(r"compactified $\omega$: $\eta=\frac{2}{\pi}\arctan(\omega/\omega_s)$")
    axes[1].legend(loc="upper left", fontsize=9, framealpha=0.9)

    fig.suptitle(
        (
            f"Teukolsky safe-domain scan  "
            f"(l,m,s)=({settings.l},{settings.m},{settings.s}),  "
            f"a∈[{settings.a_min:.0e},{settings.a_max:.6f}],  "
            f"ω∈[{settings.omega_min:.0e},∞)"
        ),
        fontsize=14,
    )
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Scan and plot Teukolsky safe domains.")
    parser.add_argument("--l", type=int, default=2)
    parser.add_argument("--m", type=int, default=2)
    parser.add_argument("--s", type=int, default=-2)
    parser.add_argument("--a-min", type=float, default=1.0e-5)
    parser.add_argument("--a-max", type=float, default=0.999999)
    parser.add_argument("--n-a", type=int, default=48)
    parser.add_argument("--omega-min", type=float, default=1.0e-6)
    parser.add_argument("--n-omega", type=int, default=72)
    parser.add_argument("--omega-scale", type=float, default=1.0)
    parser.add_argument("--eta-max", type=float, default=0.995)
    parser.add_argument("--k-margin", type=float, default=1.0e-2)
    parser.add_argument("--r-match", type=float, default=8.0)
    parser.add_argument("--n-cheb", type=int, default=32)
    parser.add_argument("--mma-ref-mode", choices=["relative", "fixed"], default="relative")
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=["lambda", "ramp", "mma"],
        default=["lambda", "ramp", "mma"],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "tests" / "outputs" / "safe_domain_scan",
    )
    parser.add_argument(
        "--mma-kernel-path",
        type=str,
        default="/mnt/f/mma/WolframKernel.exe",
    )
    parser.add_argument(
        "--mma-wl-path-win",
        type=str,
        default="F:/EMRI/Radial_flow/Radial_Function.wl",
    )
    parser.add_argument("--mma-timeout-sec", type=float, default=8.0)
    args = parser.parse_args()

    settings = ScanSettings(
        l=args.l,
        m=args.m,
        s=args.s,
        a_min=args.a_min,
        a_max=args.a_max,
        n_a=args.n_a,
        omega_min=args.omega_min,
        n_omega=args.n_omega,
        omega_scale=args.omega_scale,
        eta_max=args.eta_max,
        k_margin=args.k_margin,
        r_match=args.r_match,
        n_cheb=args.n_cheb,
        mma_ref_mode=args.mma_ref_mode,
        mma_kernel_path=args.mma_kernel_path,
        mma_wl_path_win=args.mma_wl_path_win,
        mma_timeout_sec=args.mma_timeout_sec,
    )

    a_vals = build_a_grid(args.a_min, args.a_max, args.n_a)
    eta_vals = build_eta_grid(args.omega_min, args.n_omega, args.omega_scale, args.eta_max)
    omega_vals = np.asarray(omega_from_eta(eta_vals, args.omega_scale), dtype=float)

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    mode_results: list[tuple[str, np.ndarray, Counter]] = []
    mode_code_grids: dict[str, np.ndarray] = {}
    mode_fail_examples: dict[str, dict[str, list[dict[str, float | str]]]] = {}

    mma_sampler = None
    try:
        if "mma" in args.modes:
            mma_sampler = MathematicaRinSampler(
                kernel_path=args.mma_kernel_path,
                wl_path_win=args.mma_wl_path_win,
                timeout_sec=args.mma_timeout_sec,
            )

        for mode_name in args.modes:
            if mode_name == "lambda":
                classify_fn = lambda a, omega: classify_lambda(a, omega, args.l, args.m, args.s)
            elif mode_name == "ramp":
                classify_fn = lambda a, omega: classify_ramp(
                    a,
                    omega,
                    args.l,
                    args.m,
                    args.s,
                    r_match=args.r_match,
                    n_cheb=args.n_cheb,
                )
            elif mode_name == "mma":
                if mma_sampler is None:
                    raise RuntimeError("Mathematica sampler was not initialized.")
                classify_fn = lambda a, omega: classify_mma(
                    mma_sampler,
                    a,
                    omega,
                    args.l,
                    args.m,
                    args.s,
                    mma_ref_mode=args.mma_ref_mode,
                )
            else:
                raise ValueError(f"Unknown mode {mode_name}")

            safe_mask, code_grid, counts, failure_examples = scan_mode(
                mode_name=mode_name,
                classify_fn=classify_fn,
                a_vals=a_vals,
                omega_vals=omega_vals,
            )
            mode_results.append((mode_name, safe_mask, counts))
            mode_code_grids[mode_name] = code_grid
            mode_fail_examples[mode_name] = failure_examples

    finally:
        if mma_sampler is not None:
            mma_sampler.close()

    plot_modes = list(mode_results)
    if len(plot_modes) != 3:
        name_to_result = {name: (mask, counts) for name, mask, counts in plot_modes}
        full_modes = []
        for name in ("lambda", "ramp", "mma"):
            if name in name_to_result:
                full_modes.append((name, *name_to_result[name]))
            else:
                full_modes.append((name, np.zeros((len(a_vals), len(omega_vals)), dtype=bool), Counter({"not_scanned": len(a_vals) * len(omega_vals)})))
        plot_modes = full_modes

    png_path = out_dir / "safe_domain_triptych.png"
    plot_triptych(
        out_png=png_path,
        settings=settings,
        a_vals=a_vals,
        eta_vals=eta_vals,
        mode_data=plot_modes,
    )

    np.savez_compressed(
        out_dir / "safe_domain_scan_data.npz",
        a_vals=a_vals,
        eta_vals=eta_vals,
        omega_vals=omega_vals,
        **{f"{name}_safe": mask for name, mask, _ in mode_results},
    )

    summary = {
        "settings": asdict(settings),
        "output_png": str(png_path),
        "counts": {
            name: {str(k): int(v) for k, v in counts.items()}
            for name, _, counts in mode_results
        },
        "failure_examples": mode_fail_examples,
        "omega_tick_examples": [1.0e-6, 1.0e-4, 1.0e-2, 1.0, 10.0, 100.0],
    }
    with open(out_dir / "safe_domain_scan_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"saved figure -> {png_path}")
    print(f"saved summary -> {out_dir / 'safe_domain_scan_summary.json'}")
    for name, mask, counts in mode_results:
        print(f"[{name}] safe={int(mask.sum())}/{mask.size}, counts={dict(counts)}")


if __name__ == "__main__":
    main()

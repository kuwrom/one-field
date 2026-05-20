"""
Display helpers for consistent terminal output.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

SEP = "═" * 82
THIN = "─" * 82


def H(title: str) -> None:
    """Print a major section header."""
    print(f"\n{SEP}\n  {title}\n{SEP}")


def S(title: str) -> None:
    """Print a subsection header."""
    print(f"\n  ── {title} " + "─" * max(0, 70 - len(title)))


def box(lines: list[str]) -> None:
    """Print a framed box."""
    w = max(len(l) for l in lines) + 4
    print("  ┌" + "─" * w + "┐")
    for l in lines:
        print(f"  │  {l:<{w-2}}│")
    print("  └" + "─" * w + "┘")


def pct(pred: float, ref: float) -> float:
    """Percentage error."""
    return 100.0 * (pred - ref) / ref

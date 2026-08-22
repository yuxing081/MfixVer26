"""Plot standardized JKR velocity-sweep metrics."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence


REQUIRED_SWEEP_FIELDS = (
    "impact_velocity",
    "rebound_velocity",
    "effective_restitution",
    "captured",
    "maximum_overlap",
    "contact_duration",
)


def sweep_series(
    rows: Sequence[Mapping[str, float]],
) -> dict[str, tuple[float, ...]]:
    if not rows:
        raise ValueError("velocity sweep contains no rows")
    missing = tuple(
        field
        for field in REQUIRED_SWEEP_FIELDS
        if any(field not in row for row in rows)
    )
    if missing:
        raise ValueError(f"velocity sweep is missing fields: {missing}")
    ordered = sorted(rows, key=lambda row: float(row["impact_velocity"]))
    impact = tuple(float(row["impact_velocity"]) for row in ordered)
    if len(set(impact)) != len(impact):
        raise ValueError("impact velocities must be unique")
    return {
        field: tuple(float(row[field]) for row in ordered)
        for field in REQUIRED_SWEEP_FIELDS
    }


def plot_velocity_sweep(
    rows: Sequence[Mapping[str, float]],
    output_directory: Path,
    file_prefix: str,
) -> tuple[Path, ...]:
    """Write the five standard velocity-sweep plots.

    Matplotlib is imported only when plotting is requested, so computation and
    CSV comparison do not depend on a graphics stack.
    """
    series = sweep_series(rows)
    output_directory.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("matplotlib is required to create plots") from error

    definitions = (
        (
            "rebound_velocity",
            series["rebound_velocity"],
            "Rebound velocity",
            "Rebound velocity (m/s)",
            "line",
        ),
        (
            "effective_restitution",
            series["effective_restitution"],
            "Effective restitution",
            "Effective restitution (-)",
            "line",
        ),
        (
            "capture_state",
            series["captured"],
            "Capture/rebound state",
            "Captured (1=yes)",
            "step",
        ),
        (
            "maximum_overlap",
            series["maximum_overlap"],
            "Maximum overlap",
            "Maximum overlap (m)",
            "line",
        ),
        (
            "contact_duration",
            series["contact_duration"],
            "Contact duration",
            "Contact duration (s)",
            "line",
        ),
    )
    paths = []
    for suffix, values, title, ylabel, plot_kind in definitions:
        figure, axis = plt.subplots(figsize=(6.4, 4.0))
        if plot_kind == "step":
            axis.step(series["impact_velocity"], values, where="mid")
            axis.set_yticks((0, 1))
        else:
            axis.plot(series["impact_velocity"], values, marker="o")
        axis.set_xlabel("Impact velocity (m/s)")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(True, alpha=0.3)
        figure.tight_layout()
        path = output_directory / f"{file_prefix}_{suffix}.png"
        figure.savefig(path, dpi=160)
        plt.close(figure)
        paths.append(path)
    return tuple(paths)

def plot_contact_validation(
    mfix_rows: Sequence[Mapping[str, object]],
    reference_rows: Sequence[Mapping[str, object]],
    output_directory: Path,
    file_prefix: str,
) -> tuple[Path, ...]:
    """Plot MFIX contact evidence against an independent Python trajectory."""
    if not mfix_rows or not reference_rows:
        raise ValueError("contact plotting requires MFIX and reference rows")
    output_directory.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("matplotlib is required to create plots") from error

    paths = []

    def save(figure, suffix: str) -> None:
        figure.tight_layout()
        path = output_directory / f"{file_prefix}_{suffix}.png"
        figure.savefig(path, dpi=180)
        plt.close(figure)
        paths.append(path)

    figure, axis = plt.subplots(figsize=(6.6, 4.2))
    axis.plot(
        [float(row["time"]) for row in mfix_rows],
        [float(row["normal_speed"]) for row in mfix_rows],
        label="MFIX26",
    )
    axis.plot(
        [float(row["time"]) for row in reference_rows],
        [float(row["normal_speed"]) for row in reference_rows],
        "--",
        label="Python JKR reference",
    )
    axis.set(xlabel="Time (s)", ylabel="Normal separation velocity (m/s)")
    axis.grid(True, alpha=0.3)
    axis.legend()
    save(figure, "velocity_time")

    figure, axis = plt.subplots(figsize=(6.6, 4.2))
    axis.plot(
        [float(row["time"]) for row in mfix_rows],
        [float(row["angular_speed_1"]) for row in mfix_rows],
        label="MFIX26 particle 1",
    )
    if any(row.get("angular_speed_2") is not None for row in mfix_rows):
        axis.plot(
            [float(row["time"]) for row in mfix_rows],
            [float(row["angular_speed_2"]) for row in mfix_rows],
            label="MFIX26 particle 2",
        )
    axis.plot(
        [float(row["time"]) for row in reference_rows],
        [0.0 for _ in reference_rows],
        "--",
        label="Python normal reference",
    )
    axis.set(xlabel="Time (s)", ylabel="Angular speed (rad/s)")
    axis.grid(True, alpha=0.3)
    axis.legend()
    save(figure, "angular_velocity_time")

    mfix_active = tuple(
        row for row in mfix_rows if row.get("active")
    )
    reference_active = tuple(
        row for row in reference_rows if row.get("active")
    )
    figure, axis = plt.subplots(figsize=(6.6, 4.2))
    axis.plot(
        [float(row["overlap"]) for row in mfix_active],
        [float(row["observed_normal_force"]) for row in mfix_active],
        label="MFIX26",
    )
    axis.plot(
        [float(row["overlap"]) for row in reference_active],
        [float(row["normal_force"]) for row in reference_active],
        "--",
        label="Python JKR reference",
    )
    axis.set(xlabel="Overlap (m)", ylabel="Total normal force (N)")
    axis.grid(True, alpha=0.3)
    axis.legend()
    save(figure, "normal_force_overlap")

    mfix_radius = tuple(
        row for row in mfix_rows if row.get("contact_radius") is not None
    )
    reference_radius = tuple(
        row for row in reference_rows if row.get("contact_radius") is not None
    )
    figure, axis = plt.subplots(figsize=(6.6, 4.2))
    axis.plot(
        [float(row["time"]) for row in mfix_radius],
        [float(row["contact_radius"]) for row in mfix_radius],
        label="MFIX26 overlap-inferred",
    )
    axis.plot(
        [float(row["time"]) for row in reference_radius],
        [float(row["contact_radius"]) for row in reference_radius],
        "--",
        label="Python JKR reference",
    )
    axis.set(xlabel="Time (s)", ylabel="Contact radius (m)")
    axis.grid(True, alpha=0.3)
    axis.legend()
    save(figure, "contact_radius_time")
    return tuple(paths)

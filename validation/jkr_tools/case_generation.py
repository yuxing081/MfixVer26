"""Generate and validate MFIX 26.1.2 JKR cases."""

from __future__ import annotations

import math
import re
from pathlib import Path
import shutil

from .config import CaseConfig, MaterialConfig, NumericalConfig, ParticleConfig


def hertz_contact_time(
    radius: float,
    density: float,
    young_modulus: float,
    poisson_ratio: float,
    wall_young_modulus: float,
    wall_poisson_ratio: float,
    normal_speed: float,
    wall_contact: bool,
) -> float:
    if min(radius, density, young_modulus, normal_speed) <= 0.0:
        raise ValueError("contact-time inputs must be positive")
    mass = (4.0 / 3.0) * math.pi * radius**3 * density
    effective_mass = mass if wall_contact else mass / 2.0
    effective_radius = radius if wall_contact else radius / 2.0
    if wall_contact:
        effective_young = 1.0 / (
            (1.0 - poisson_ratio**2) / young_modulus
            + (1.0 - wall_poisson_ratio**2) / wall_young_modulus
        )
    else:
        effective_young = young_modulus / (2.0 * (1.0 - poisson_ratio**2))
    stiffness = (4.0 / 3.0) * math.sqrt(effective_radius) * effective_young
    return 3.21 * (effective_mass / stiffness) ** 0.4 * normal_speed ** (-0.2)


def case_contact_time(case: CaseConfig, material: MaterialConfig) -> float:
    return hertz_contact_time(
        material.particle_radius,
        material.particle_density,
        material.young_modulus,
        material.poisson_ratio,
        material.wall_young_modulus,
        material.wall_poisson_ratio,
        case.contact_time_speed,
        case.wall_contact,
    )

def boolean_text(value: bool) -> str:
    return ".True." if value else ".False."


def particle_input_text(particles: tuple[ParticleConfig, ...]) -> str:
    if not particles:
        raise ValueError("a case must contain particles")
    user_count = 1 + len(particles[0].diagnostic_values)
    if any(1 + len(item.diagnostic_values) != user_count for item in particles):
        raise ValueError("all particles must use equal user-scalar columns")
    user_columns = " ".join(
        f"| usr_{index}" for index in range(1, user_count + 1)
    )
    rows = []
    for particle in particles:
        values = (
            *particle.position,
            1,
            particle.diameter,
            particle.density,
            *particle.velocity,
            particle.surface_energy,
            *particle.diagnostic_values,
        )
        rows.append(" ".join(
            str(value) if isinstance(value, int) else f"{value:.12e}"
            for value in values
        ))
    header = f"""Version 2.0
========================================================================
Instructions:
Dimension_is_2_or_3
Particle_count_follows
Read_flags_are_T_or_F
Coordinates_are_required
Phase_ID_is_scalar
Diameter_is_scalar
Density_is_scalar
Velocity_has_three_components
Temperature_is_optional
Species_are_optional
User_scalars_are_optional
Unused_values_use_defaults
All_particles_share_columns
Data_starts_after_header
========================================================================
Dimension:    3
Particles:    {len(particles)}
========================================================================
Variable      Read(T/F)     Default value (when Read=F)
========================================================================
Coordinates:     T          Must always be T
Phase_ID:        T
Diameter:        T
Density:         T
Velocity:        T
Temperature:     F          293.15
Species:         F
User_Scalar:     T
========================================================================
| x | y | z | Phase_ID | Diameter | Density | u | v | w {user_columns} |
========================================================================
"""
    return header + "\n".join(rows) + "\n"


def _add_domain(lines: list[str], prefix: str, domain: float) -> None:
    for axis, low, high in (("x", "w", "e"), ("y", "s", "n"), ("z", "b", "t")):
        lines.extend((
            f"{prefix}_{axis}_{low}(1) = 0.0",
            f"{prefix}_{axis}_{high}(1) = {domain:.12e}",
        ))


def _monitor_lines(case: CaseConfig, domain: float) -> list[str]:
    lines = [
        "monitor_name(1) = 'PARTICLE_1'",
        "monitor_type(1) = 104",
        f"monitor_dt(1) = {case.monitor_interval:.12e}",
    ]
    _add_domain(lines, "monitor", domain)
    lines.extend(f"monitor_{field}(1) = .True." for field in (
        "vel_x", "vel_y", "vel_z", "rot_x", "rot_y", "rot_z"
    ))
    lines.extend(
        f"monitor_des_usr_var(1,{index}) = .True."
        for index in case.monitor_user_variables
    )
    return lines

def mfix_dat_text(
    case: CaseConfig,
    material: MaterialConfig,
    numerical: NumericalConfig,
    run_name: str,
    run_type: str,
    stop_time: float,
    partition: tuple[int, int, int],
    include_jkr_keywords: bool = True,
) -> str:
    if run_type not in {"new", "restart_1", "restart_2"}:
        raise ValueError(f"unsupported run_type: {run_type}")
    if stop_time <= 0.0:
        raise ValueError("stop_time must be positive")
    if not re.fullmatch(r"[A-Za-z0-9_]+", run_name):
        raise ValueError("invalid run_name")
    contact_time = case_contact_time(case, material)
    domain = numerical.domain_size
    lines = [
        "#! MFIX 26.1.2 JKR validation input",
        "#! Generated from validation/jkr_tools/config.py; do not edit per run.",
        f"#! Classical Hertz contact time: {contact_time:.12e} s",
        f"description = '{case.description}'",
        f"run_name = '{run_name}'",
        "units = 'SI'",
        f"run_type = '{run_type}'",
        f"tstop = {stop_time:.12e}",
        f"dt = {contact_time / numerical.initial_dt_divisor:.12e}",
        f"dt_min = {contact_time / numerical.minimum_dt_divisor:.12e}",
        f"dt_max = {contact_time / numerical.maximum_dt_divisor:.12e}",
        "dt_fac = 0.9999",
        "energy_eq = .False.",
        "turbulence_model = 'NONE'",
        "call_usr = .False.",
        "gravity_x = 0.0",
        "gravity_y = 0.0",
        "gravity_z = 0.0",
        "coordinates = 'CARTESIAN'",
        "cartesian_grid = .False.",
        "use_stl = .False.",
        "x_min = 0.0", f"x_max = {domain:.12e}",
        "y_min = 0.0", f"y_max = {domain:.12e}",
        "z_min = 0.0", f"z_max = {domain:.12e}",
        f"imax = {numerical.mesh_cells}",
        f"jmax = {numerical.mesh_cells}",
        f"kmax = {numerical.mesh_cells}",
        "momentum_x_eq(0) = .False.",
        "momentum_y_eq(0) = .False.",
        "momentum_z_eq(0) = .False.",
        "ro_g0 = 0.0",
        "mu_g0 = 1.8e-5",
        "species_eq(0) = .False.",
        "mmax = 1",
        "solids_model(1) = 'DEM'",
        f"d_p0(1) = {2.0 * material.particle_radius:.12e}",
        f"ro_s0(1) = {material.particle_density:.12e}",
        "species_eq(1) = .False.",
        "nmax_s(1) = 0",
        f"particles = {len(case.particles)}",
        "des_interp_scheme = 'NONE'",
        "des_coll_model = 'HERTZIAN'",
        f"mew = {material.sliding_friction:.12e}",
        f"mew_w = {material.wall_sliding_friction:.12e}",
        f"mew_r = {material.rolling_friction:.12e}",
        f"mew_rw = {material.wall_rolling_friction:.12e}",
        f"des_en_input(1) = {material.normal_restitution:.12e}",
        f"des_en_wall_input(1) = {material.normal_restitution:.12e}",
        f"des_et_input(1) = {material.tangential_restitution:.12e}",
        f"des_et_wall_input(1) = {material.tangential_restitution:.12e}",
        "des_etat_fac = 0.5",
        "des_etat_w_fac = 0.5",
        f"e_young(1) = {material.young_modulus:.12e}",
        f"v_poisson(1) = {material.poisson_ratio:.12e}",
        f"ew_young = {material.wall_young_modulus:.12e}",
        f"vw_poisson = {material.wall_poisson_ratio:.12e}",
        "use_cohesion = .False.",
        "van_der_waals = .False.",
        f"dtsolid_factor = {numerical.dtsolid_factor:.12e}",
        f"des_usr_var_size = {case.des_usr_var_size}",
        f"des_neighbor_search = {numerical.neighbor_search}",
        f"neighbor_search_n = {numerical.neighbor_search_interval}",
        f"factor_rlm = {numerical.neighbor_radius_factor:.12e}",
    ]
    _add_domain(lines, "ic", domain)
    lines.extend((
        "ic_ep_g(1) = 1.0",
        "ic_ep_s(1,1) = 0.0",
        "detect_stall = .False.",
        "max_nit = 50",
        f"res_dt = {case.restart_interval:.12e}",
        "write_part_out = .True.",
        f"write_vtk_files = {boolean_text(case.write_vtk_files)}",
        "time_dependent_filename = .True.",
    ))
    _add_domain(lines, "vtk", domain)
    lines.extend((
        "vtk_filebase(1) = 'PARTICLES'",
        f"vtk_dt(1) = {case.vtk_interval:.12e}",
        "vtk_data(1) = 'P'",
        "vtk_part_vel(1) = .True.",
        "vtk_part_angular_vel(1) = .True.",
        "vtk_part_coll_force(1) = .True.",
        "vtk_part_phase_id(1) = .True.",
    ))
    lines.extend(
        f"vtk_part_usr_var(1,{index}) = .True."
        for index in range(1, case.des_usr_var_size + 1)
    )
    lines.extend(("vtk_part_rank(1) = .True.", "vtk_part_id(1) = .True."))
    lines.extend(_monitor_lines(case, domain))
    lines.extend((
        f"nodesi = {partition[0]}",
        f"nodesj = {partition[1]}",
        f"nodesk = {partition[2]}",
    ))
    if include_jkr_keywords:
        lines.extend((
            f"des_surface_energy(1) = {material.surface_energy:.12e}",
            f"des_surface_energy_wall = {material.wall_surface_energy:.12e}",
            f"jkr_alpha_safety = {material.jkr_safety_factor:.12e}",
            f"jkr_thetar_crit = {material.rolling_angle:.12e}",
            f"use_jkr = {boolean_text(case.use_jkr)}",
            f"wall_jkr_thetar_crit = {material.wall_rolling_angle:.12e}",
        ))

    return "\n".join(lines) + "\n"
def keyword_count(mfix_text: str, keyword: str) -> int:
    pattern = re.compile(
        rf"^\s*{re.escape(keyword)}(?:\([^\n=]*\))?\s*=",
        re.IGNORECASE | re.MULTILINE,
    )
    return len(pattern.findall(mfix_text))


def particle_data_rows(particle_text: str) -> tuple[tuple[float, ...], ...]:
    rows = []
    for line in particle_text.splitlines():
        try:
            values = tuple(float(field) for field in line.split())
        except ValueError:
            continue
        if len(values) >= 9:
            rows.append(values)
    return tuple(rows)


def pair_overlaps(particles: tuple[ParticleConfig, ...]) -> tuple[float, ...]:
    return tuple(
        0.5 * (first.diameter + second.diameter)
        - math.dist(first.position, second.position)
        for index, first in enumerate(particles)
        for second in particles[index + 1 :]
    )


def wall_gap(particle: ParticleConfig) -> float:
    return particle.position[1] - 0.5 * particle.diameter


def validate_generated_case(
    mfix_text: str,
    particle_text: str,
    case: CaseConfig,
    material: MaterialConfig,
    numerical: NumericalConfig,
    include_jkr_keywords: bool = True,
) -> tuple[str, ...]:
    errors = []
    for keyword in ("run_name", "run_type", "dt", "particles"):
        count = keyword_count(mfix_text, keyword)
        if count != 1:
            errors.append(f"{keyword} occurs {count} times; expected once")
    use_jkr_count = keyword_count(mfix_text, "use_jkr")
    if include_jkr_keywords:
        if use_jkr_count != 1:
            errors.append(f"use_jkr occurs {use_jkr_count} times; expected once")
        if f"use_jkr = {boolean_text(case.use_jkr)}".lower() not in mfix_text.lower():
            errors.append("use_jkr does not match case purpose")
    elif use_jkr_count != 0:
        errors.append("official-compatible input contains use_jkr")
    elif case.use_jkr:
        errors.append("cannot omit JKR keywords for a JKR-enabled case")
    rows = particle_data_rows(particle_text)
    if len(rows) != len(case.particles):
        errors.append(f"particle count is {len(rows)}; expected {len(case.particles)}")
    expected_columns = 9 + case.des_usr_var_size
    if any(len(row) != expected_columns for row in rows):
        errors.append(f"particle rows must contain {expected_columns} columns")
    for particle in case.particles:
        radius = particle.diameter / 2.0
        if any(position - radius < -1.0e-15 for position in particle.position):
            errors.append("particle extends below domain minimum")
        if any(
            position + radius > numerical.domain_size + 1.0e-15
            for position in particle.position
        ):
            errors.append("particle extends above domain maximum")
    overlaps = pair_overlaps(case.particles)
    if overlaps and abs(max(overlaps) - case.expected_initial_overlap) > 1.0e-15:
        errors.append("pair geometry differs from expected overlap")
    if case.wall_contact:
        minimum_gap = min(wall_gap(particle) for particle in case.particles)
        if abs(minimum_gap - case.expected_initial_overlap) > 1.0e-15:
            errors.append("wall geometry differs from expected gap")
    contact_time = case_contact_time(case, material)
    if contact_time / numerical.initial_dt_divisor > contact_time / 1.0e4:
        errors.append("time step has fewer than 10000 contact-time steps")
    for keyword in (
        "write_part_out", "monitor_vel_x", "monitor_rot_x", "vtk_part_id"
    ):
        if keyword_count(mfix_text, keyword) != 1:
            errors.append(f"required output {keyword} is missing")
    return tuple(errors)


def write_case_files(
    output_directory: Path,
    case: CaseConfig,
    material: MaterialConfig,
    numerical: NumericalConfig,
    run_name: str,
    run_type: str,
    stop_time: float,
    partition: tuple[int, int, int],
    include_jkr_keywords: bool = True,
) -> None:
    mfix_text = mfix_dat_text(
        case,
        material,
        numerical,
        run_name,
        run_type,
        stop_time,
        partition,
        include_jkr_keywords,
    )
    particle_text = particle_input_text(case.particles)
    errors = validate_generated_case(
        mfix_text, particle_text, case, material, numerical, include_jkr_keywords
    )
    if errors:
        raise ValueError("; ".join(errors))
    output_directory.mkdir(parents=True, exist_ok=False)
    (output_directory / "mfix.dat").write_text(mfix_text, encoding="utf-8")
    (output_directory / "particle_input.dat").write_text(
        particle_text, encoding="utf-8"
    )

def canonical_case_directory(cases_root: Path, case: CaseConfig) -> Path:
    return cases_root / case.category / case.name


def canonical_run_name(case: CaseConfig) -> str:
    return f"R26_{case.name.upper()}"


def canonical_case_texts(
    case: CaseConfig,
    material: MaterialConfig,
    numerical: NumericalConfig,
) -> tuple[str, str]:
    include_jkr_keywords = case.use_jkr
    mfix_text = mfix_dat_text(
        case=case,
        material=material,
        numerical=numerical,
        run_name=canonical_run_name(case),
        run_type="new",
        stop_time=case.stop_time,
        partition=numerical.serial_partition,
        include_jkr_keywords=include_jkr_keywords,
    )
    particle_text = particle_input_text(case.particles)
    errors = validate_generated_case(
        mfix_text,
        particle_text,
        case,
        material,
        numerical,
        include_jkr_keywords,
    )
    if errors:
        raise ValueError("; ".join(errors))
    return mfix_text, particle_text


def write_canonical_case_asset(
    cases_root: Path,
    case: CaseConfig,
    material: MaterialConfig,
    numerical: NumericalConfig,
) -> Path:
    """Write only the two deterministic, reviewable case input assets."""
    output_directory = canonical_case_directory(cases_root, case)
    output_directory.mkdir(parents=True, exist_ok=True)
    mfix_text, particle_text = canonical_case_texts(case, material, numerical)
    (output_directory / "mfix.dat").write_text(mfix_text, encoding="utf-8")
    (output_directory / "particle_input.dat").write_text(
        particle_text, encoding="utf-8"
    )
    return output_directory


def canonical_case_sync_errors(
    cases_root: Path,
    case: CaseConfig,
    material: MaterialConfig,
    numerical: NumericalConfig,
) -> tuple[str, ...]:
    output_directory = canonical_case_directory(cases_root, case)
    expected = dict(zip(
        ("mfix.dat", "particle_input.dat"),
        canonical_case_texts(case, material, numerical),
    ))
    errors = []
    for name, text in expected.items():
        path = output_directory / name
        if not path.is_file():
            errors.append(f"missing canonical asset: {path}")
        elif path.read_bytes() != text.encode("utf-8"):
            errors.append(f"canonical asset differs byte-for-byte: {path}")
    return tuple(errors)


def copy_canonical_case_asset(
    cases_root: Path,
    output_directory: Path,
    case: CaseConfig,
    material: MaterialConfig,
    numerical: NumericalConfig,
) -> None:
    errors = canonical_case_sync_errors(
        cases_root, case, material, numerical
    )
    if errors:
        raise ValueError("; ".join(errors))
    output_directory.mkdir(parents=True, exist_ok=False)
    source_directory = canonical_case_directory(cases_root, case)
    for name in ("mfix.dat", "particle_input.dat"):
        shutil.copy2(source_directory / name, output_directory / name)

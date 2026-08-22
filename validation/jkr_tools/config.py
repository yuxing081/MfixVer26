"""Central MFIX 26.1.2 JKR validation configuration.

All values use SI units. Factory functions return fresh immutable settings so
calculation modules never read experiment parameters from mutable module state.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MaterialConfig:
    particle_radius: float
    particle_density: float
    young_modulus: float
    poisson_ratio: float
    wall_young_modulus: float
    wall_poisson_ratio: float
    normal_restitution: float
    tangential_restitution: float
    sliding_friction: float
    wall_sliding_friction: float
    rolling_friction: float
    wall_rolling_friction: float
    surface_energy: float
    wall_surface_energy: float
    rolling_angle: float
    wall_rolling_angle: float
    jkr_safety_factor: float


@dataclass(frozen=True)
class NumericalConfig:
    domain_size: float
    mesh_cells: int
    dtsolid_factor: float
    initial_dt_divisor: float
    minimum_dt_divisor: float
    maximum_dt_divisor: float
    neighbor_search: int
    neighbor_search_interval: int
    neighbor_radius_factor: float
    serial_partition: tuple[int, int, int]
    dmp_partition: tuple[int, int, int]
    dmp_ranks: int
    coarse_velocity_factors: tuple[float, ...]
    critical_refinement_points: int


@dataclass(frozen=True)
class ParticleConfig:
    position: tuple[float, float, float]
    velocity: tuple[float, float, float]
    diameter: float
    density: float
    surface_energy: float
    diagnostic_values: tuple[float, ...] = ()


@dataclass(frozen=True)
class CaseConfig:
    name: str
    category: str
    description: str
    particles: tuple[ParticleConfig, ...]
    wall_contact: bool
    use_jkr: bool
    stop_time: float
    restart_time: float | None
    monitor_interval: float
    restart_interval: float
    vtk_interval: float
    des_usr_var_size: int
    monitor_user_variables: tuple[int, ...]
    expected_initial_overlap: float
    contact_time_speed: float
    write_vtk_files: bool = False


def default_material_config() -> MaterialConfig:
    return MaterialConfig(
        particle_radius=1.0e-6,
        particle_density=2500.0,
        young_modulus=2.0e8,
        poisson_ratio=0.30,
        wall_young_modulus=1.0e20,
        wall_poisson_ratio=0.30,
        normal_restitution=0.60,
        tangential_restitution=0.60,
        sliding_friction=0.30,
        wall_sliding_friction=0.30,
        rolling_friction=0.01,
        wall_rolling_friction=0.01,
        surface_energy=0.015,
        wall_surface_energy=0.015,
        rolling_angle=0.01,
        wall_rolling_angle=0.01,
        jkr_safety_factor=1.0,
    )


def default_numerical_config() -> NumericalConfig:
    return NumericalConfig(
        domain_size=20.0e-6,
        mesh_cells=6,
        dtsolid_factor=100.0,
        initial_dt_divisor=1.0e4,
        minimum_dt_divisor=1.0e6,
        maximum_dt_divisor=10.0,
        neighbor_search=4,
        neighbor_search_interval=20,
        neighbor_radius_factor=1.2,
        serial_partition=(1, 1, 1),
        dmp_partition=(2, 1, 1),
        dmp_ranks=2,
        coarse_velocity_factors=(
            0.25,
            0.50,
            0.75,
            0.90,
            1.0,
            1.10,
            1.25,
            1.50,
            2.0,
            3.0,
        ),
        critical_refinement_points=7,
    )


def restart_case_configs(
    material: MaterialConfig,
) -> tuple[CaseConfig, ...]:
    pair_particles = (
        ParticleConfig(
            position=(9.99e-6, 4.055e-6, 5.0e-6),
            velocity=(1.1, 0.0, 0.0),
            diameter=2.0e-6,
            density=material.particle_density,
            surface_energy=material.surface_energy,
        ),
        ParticleConfig(
            position=(9.99e-6, 5.945e-6, 5.0e-6),
            velocity=(0.9, 0.0, 0.0),
            diameter=1.8e-6,
            density=material.particle_density,
            surface_energy=material.surface_energy,
        ),
    )
    wall_particles = (
        ParticleConfig(
            position=(10.0e-6, 1.0e-6, 10.0e-6),
            velocity=(0.8, -0.8, 0.0),
            diameter=2.0e-6,
            density=material.particle_density,
            surface_energy=material.surface_energy,
            diagnostic_values=(0.0, 0.0, 0.0),
        ),
    )
    return (
        CaseConfig(
            name="restart_pair_pfr",
            category="restart",
            description="unequal-radius pair with active PFT, PFR, and JKR neck",
            particles=pair_particles,
            wall_contact=False,
            use_jkr=True,
            stop_time=8.0e-8,
            restart_time=5.0e-9,
            monitor_interval=1.0e-9,
            restart_interval=5.0e-9,
            vtk_interval=1.0e-9,
            des_usr_var_size=1,
            monitor_user_variables=(),
            expected_initial_overlap=1.0e-8,
            contact_time_speed=1.0,
        ),
        CaseConfig(
            name="restart_wall_oblique",
            category="restart",
            description="oblique particle-wall contact with wall history",
            particles=wall_particles,
            wall_contact=True,
            use_jkr=True,
            stop_time=6.0e-7,
            restart_time=2.0e-7,
            monitor_interval=1.0e-9,
            restart_interval=2.0e-7,
            vtk_interval=6.0e-7,
            des_usr_var_size=4,
            monitor_user_variables=(2, 3, 4),
            expected_initial_overlap=0.0,
            contact_time_speed=0.8,
        ),
    )


def restart_case_config(
    case_name: str,
    material: MaterialConfig,
) -> CaseConfig:
    matches = tuple(
        case for case in restart_case_configs(material) if case.name == case_name
    )
    if len(matches) != 1:
        raise ValueError(f"Unknown restart case: {case_name}")
    return matches[0]

def hertz_case_configs(
    material: MaterialConfig,
) -> tuple[CaseConfig, ...]:
    diameter = 2.0 * material.particle_radius
    pair_particles = (
        ParticleConfig(
            position=(10.0e-6, 8.9e-6, 10.0e-6),
            velocity=(0.0, 0.5, 0.0),
            diameter=diameter,
            density=material.particle_density,
            surface_energy=material.surface_energy,
        ),
        ParticleConfig(
            position=(10.0e-6, 11.1e-6, 10.0e-6),
            velocity=(0.0, -0.5, 0.0),
            diameter=diameter,
            density=material.particle_density,
            surface_energy=material.surface_energy,
        ),
    )
    wall_normal_particle = (
        ParticleConfig(
            position=(10.0e-6, 1.2e-6, 10.0e-6),
            velocity=(0.0, -1.0, 0.0),
            diameter=diameter,
            density=material.particle_density,
            surface_energy=material.surface_energy,
        ),
    )
    wall_oblique_particle = (
        ParticleConfig(
            position=(10.0e-6, 1.2e-6, 10.0e-6),
            velocity=(0.5, -1.0, 0.0),
            diameter=diameter,
            density=material.particle_density,
            surface_energy=material.surface_energy,
        ),
    )

    def make_case(
        name: str,
        description: str,
        particles: tuple[ParticleConfig, ...],
        wall_contact: bool,
        initial_gap_or_overlap: float,
    ) -> CaseConfig:
        return CaseConfig(
            name=name,
            category="native_hertz",
            description=description,
            particles=particles,
            wall_contact=wall_contact,
            use_jkr=False,
            stop_time=5.0e-7,
            restart_time=None,
            monitor_interval=1.0e-9,
            restart_interval=5.0e-7,
            vtk_interval=1.0e-9,
            des_usr_var_size=1,
            monitor_user_variables=(),
            expected_initial_overlap=initial_gap_or_overlap,
            contact_time_speed=1.0,
            write_vtk_files=True,
        )

    return (
        make_case(
            name="pair_hertz_normal",
            description="native Hertz equal-particle normal collision",
            particles=pair_particles,
            wall_contact=False,
            initial_gap_or_overlap=-0.2e-6,
        ),
        make_case(
            name="wall_hertz_normal",
            description="native Hertz particle-wall normal collision",
            particles=wall_normal_particle,
            wall_contact=True,
            initial_gap_or_overlap=0.2e-6,
        ),
        make_case(
            name="wall_hertz_oblique",
            description="native Hertz particle-wall oblique friction collision",
            particles=wall_oblique_particle,
            wall_contact=True,
            initial_gap_or_overlap=0.2e-6,
        ),
    )


def jkr_normal_case_configs(
    material: MaterialConfig,
    wall_critical_speed: float,
    pair_critical_speed: float,
) -> tuple[CaseConfig, ...]:
    """Return normal cases scaled from the adhesion-only critical speed."""
    if wall_critical_speed <= 0.0 or pair_critical_speed <= 0.0:
        raise ValueError("critical speeds must be positive")
    diameter = 2.0 * material.particle_radius

    def wall_case(name: str, speed: float, behavior: str) -> CaseConfig:
        particle = ParticleConfig(
            position=(10.0e-6, 1.2e-6, 10.0e-6),
            velocity=(0.0, -speed, 0.0),
            diameter=diameter,
            density=material.particle_density,
            surface_energy=material.surface_energy,
        )
        return CaseConfig(
            name=name,
            category="particle_wall",
            description=f"JKR wall normal {behavior} at {speed:.12e} m/s",
            particles=(particle,),
            wall_contact=True,
            use_jkr=True,
            stop_time=2.0e-6,
            restart_time=None,
            monitor_interval=1.0e-9,
            restart_interval=2.0e-6,
            vtk_interval=2.0e-9,
            des_usr_var_size=1,
            monitor_user_variables=(),
            expected_initial_overlap=0.2e-6,
            contact_time_speed=speed,
            write_vtk_files=True,
        )

    def pair_case(name: str, speed: float, behavior: str) -> CaseConfig:
        particles = (
            ParticleConfig(
                position=(10.0e-6, 8.9e-6, 10.0e-6),
                velocity=(0.0, 0.5 * speed, 0.0),
                diameter=diameter,
                density=material.particle_density,
                surface_energy=material.surface_energy,
            ),
            ParticleConfig(
                position=(10.0e-6, 11.1e-6, 10.0e-6),
                velocity=(0.0, -0.5 * speed, 0.0),
                diameter=diameter,
                density=material.particle_density,
                surface_energy=material.surface_energy,
            ),
        )
        return CaseConfig(
            name=name,
            category="particle_pair",
            description=f"JKR pair normal {behavior} at {speed:.12e} m/s",
            particles=particles,
            wall_contact=False,
            use_jkr=True,
            stop_time=2.0e-6,
            restart_time=None,
            monitor_interval=1.0e-9,
            restart_interval=2.0e-6,
            vtk_interval=2.0e-9,
            des_usr_var_size=1,
            monitor_user_variables=(),
            expected_initial_overlap=-0.2e-6,
            contact_time_speed=speed,
            write_vtk_files=True,
        )

    return (
        wall_case(
            "wall_normal_capture",
            0.95 * wall_critical_speed,
            "capture",
        ),
        wall_case(
            "wall_normal_rebound",
            3.0 * wall_critical_speed,
            "rebound",
        ),
        pair_case(
            "pair_normal_capture",
            0.95 * pair_critical_speed,
            "capture",
        ),
        pair_case(
            "pair_normal_rebound",
            3.0 * pair_critical_speed,
            "rebound",
        ),
    )

def comparison_tolerances() -> dict[str, tuple[float, float]]:
    """Return initial absolute/relative tolerances by data source and field."""
    return {
        "double_position": (1.0e-13, 1.0e-10),
        "double_velocity": (1.0e-10, 1.0e-9),
        "double_angular_velocity": (1.0e-4, 1.0e-9),
        "double_history": (1.0e-14, 1.0e-9),
        "monitor_position": (1.0e-11, 1.0e-7),
        "monitor_velocity": (1.0e-7, 1.0e-6),
        "vtp_position": (2.0e-11, 1.0e-7),
        "vtp_velocity": (1.0e-4, 1.0e-5),
        "vtp_angular_velocity": (4.0, 1.0e-5),
        "vtp_force": (1.0e-8, 1.0e-4),
    }

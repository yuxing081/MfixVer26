"""Independent scalar JKR reference functions for MFIX 26.1.2.

The production model stores a surface energy gamma and combines unlike
surfaces as sqrt(gamma_i * gamma_j). For comparison with literature that
uses work of adhesion, this module records W = 2 * gamma_eff explicitly.
Positive scalar normal force is compressive; negative force is tensile.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ContactParameters:
    effective_surface_energy: float
    work_of_adhesion: float
    effective_radius: float
    effective_mass: float
    effective_young_modulus: float
    effective_shear_modulus: float
    equilibrium_contact_radius: float
    critical_separation: float
    pull_off_force_magnitude: float
    normal_stiffness_per_radius: float
    tangential_stiffness_per_radius: float
    damping_alpha: float


def particle_mass(radius: float, density: float) -> float:
    if radius <= 0.0 or density <= 0.0:
        raise ValueError("Particle radius and density must be positive")
    return (4.0 / 3.0) * math.pi * radius**3 * density


def effective_surface_energy(
    first_surface_energy: float,
    second_surface_energy: float,
) -> float:
    if first_surface_energy <= 0.0 or second_surface_energy <= 0.0:
        raise ValueError("JKR surface energies must be positive")
    return math.sqrt(first_surface_energy * second_surface_energy)


def work_of_adhesion(effective_gamma: float) -> float:
    if effective_gamma <= 0.0:
        raise ValueError("Effective surface energy must be positive")
    return 2.0 * effective_gamma


def effective_radius(
    first_radius: float,
    second_radius: float | None,
) -> float:
    if first_radius <= 0.0:
        raise ValueError("Particle radius must be positive")
    if second_radius is None:
        return first_radius
    if second_radius <= 0.0:
        raise ValueError("Second particle radius must be positive")
    return first_radius * second_radius / (first_radius + second_radius)


def shear_modulus(young_modulus: float, poisson_ratio: float) -> float:
    if young_modulus <= 0.0:
        raise ValueError("Young modulus must be positive")
    if not -1.0 < poisson_ratio < 0.5:
        raise ValueError("Poisson ratio must be in (-1, 0.5)")
    return young_modulus / (2.0 * (1.0 + poisson_ratio))


def effective_young_modulus(
    first_young_modulus: float,
    first_poisson_ratio: float,
    second_young_modulus: float,
    second_poisson_ratio: float,
) -> float:
    compliance = (
        (1.0 - first_poisson_ratio**2) / first_young_modulus
        + (1.0 - second_poisson_ratio**2) / second_young_modulus
    )
    return 1.0 / compliance


def effective_shear_modulus(
    first_young_modulus: float,
    first_poisson_ratio: float,
    second_young_modulus: float,
    second_poisson_ratio: float,
) -> float:
    first_shear = shear_modulus(first_young_modulus, first_poisson_ratio)
    second_shear = shear_modulus(second_young_modulus, second_poisson_ratio)
    compliance = (
        (2.0 - first_poisson_ratio) / first_shear
        + (2.0 - second_poisson_ratio) / second_shear
    )
    return 1.0 / compliance


def damping_alpha(restitution: float) -> float:
    if not 0.0 <= restitution <= 1.0:
        raise ValueError("Restitution must be in [0, 1]")
    if restitution == 1.0:
        return 0.0
    return (
        1.2728
        - 4.2783 * restitution
        + 11.087 * restitution**2
        - 22.348 * restitution**3
        + 27.467 * restitution**4
        - 18.022 * restitution**5
        + 4.8218 * restitution**6
    )


def contact_parameters(
    first_radius: float,
    first_density: float,
    first_surface_energy: float,
    first_young_modulus: float,
    first_poisson_ratio: float,
    restitution: float,
    second_radius: float | None,
    second_density: float | None,
    second_surface_energy: float,
    second_young_modulus: float,
    second_poisson_ratio: float,
) -> ContactParameters:
    first_mass = particle_mass(first_radius, first_density)
    if second_radius is None:
        contact_mass = first_mass
    else:
        if second_density is None:
            raise ValueError("Second particle density is required for a pair")
        second_mass = particle_mass(second_radius, second_density)
        contact_mass = first_mass * second_mass / (first_mass + second_mass)
    contact_radius = effective_radius(first_radius, second_radius)
    contact_young = effective_young_modulus(
        first_young_modulus,
        first_poisson_ratio,
        second_young_modulus,
        second_poisson_ratio,
    )
    contact_shear = effective_shear_modulus(
        first_young_modulus,
        first_poisson_ratio,
        second_young_modulus,
        second_poisson_ratio,
    )
    gamma = effective_surface_energy(
        first_surface_energy,
        second_surface_energy,
    )
    equilibrium_radius = (
        9.0 * math.pi * gamma * contact_radius**2 / contact_young
    ) ** (1.0 / 3.0)
    critical_separation = equilibrium_radius**2 / (
        2.0 * 6.0 ** (1.0 / 3.0) * contact_radius
    )
    return ContactParameters(
        effective_surface_energy=gamma,
        work_of_adhesion=work_of_adhesion(gamma),
        effective_radius=contact_radius,
        effective_mass=contact_mass,
        effective_young_modulus=contact_young,
        effective_shear_modulus=contact_shear,
        equilibrium_contact_radius=equilibrium_radius,
        critical_separation=critical_separation,
        pull_off_force_magnitude=3.0 * math.pi * contact_radius * gamma,
        normal_stiffness_per_radius=(4.0 / 3.0) * contact_young,
        tangential_stiffness_per_radius=8.0 * contact_shear,
        damping_alpha=damping_alpha(restitution),
    )


def stable_contact_radius(
    overlap: float,
    parameters: ContactParameters,
) -> float:
    if overlap < -parameters.critical_separation:
        raise ValueError("Overlap is beyond the stable JKR tensile branch")
    minimum_ratio = (1.0 / 6.0) ** (2.0 / 3.0)
    target = overlap / (
        6.0 ** (1.0 / 3.0) * parameters.critical_separation
    )

    def residual(radius_ratio: float) -> float:
        return (
            2.0 * radius_ratio**2
            - (4.0 / 3.0) * math.sqrt(radius_ratio)
            - target
        )

    lower = minimum_ratio
    upper = max(1.0, math.sqrt(max(target, 0.0) / 2.0) + 1.0)
    while residual(upper) < 0.0:
        upper *= 2.0
    for _ in range(80):
        midpoint = 0.5 * (lower + upper)
        if residual(midpoint) < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return parameters.equilibrium_contact_radius * 0.5 * (lower + upper)


def elastic_normal_force(
    contact_radius: float,
    parameters: ContactParameters,
) -> float:
    ratio = contact_radius / parameters.equilibrium_contact_radius
    return 4.0 * parameters.pull_off_force_magnitude * (
        ratio**3 - ratio**1.5
    )


def normal_damping_coefficient(
    contact_radius: float,
    parameters: ContactParameters,
) -> float:
    return parameters.damping_alpha * math.sqrt(
        parameters.effective_mass
        * parameters.normal_stiffness_per_radius
        * contact_radius
    )


def contact_is_active(
    overlap: float,
    was_active: bool,
    critical_separation: float,
) -> bool:
    if overlap >= 0.0:
        return True
    return was_active and overlap >= -critical_separation


def adhesion_only_critical_speed(
    parameters: ContactParameters,
) -> float:
    return math.sqrt(14.18 / parameters.effective_mass) * (
        parameters.work_of_adhesion**5
        * parameters.effective_radius**4
        / parameters.effective_young_modulus**2
    ) ** (1.0 / 6.0)

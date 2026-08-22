import math
import unittest

from validation.jkr_tools.config import default_material_config
from validation.jkr_tools.reference import (
    adhesion_only_critical_speed,
    contact_is_active,
    contact_parameters,
    damping_alpha,
    effective_radius,
    effective_surface_energy,
    elastic_normal_force,
    stable_contact_radius,
)


def wall_parameters():
    material = default_material_config()
    return contact_parameters(
        first_radius=material.particle_radius,
        first_density=material.particle_density,
        first_surface_energy=material.surface_energy,
        first_young_modulus=material.young_modulus,
        first_poisson_ratio=material.poisson_ratio,
        restitution=1.0,
        second_radius=None,
        second_density=None,
        second_surface_energy=material.wall_surface_energy,
        second_young_modulus=material.wall_young_modulus,
        second_poisson_ratio=material.wall_poisson_ratio,
    )


class SurfaceEnergyTests(unittest.TestCase):
    def test_homogeneous_combination_returns_input(self):
        self.assertEqual(effective_surface_energy(0.015, 0.015), 0.015)

    def test_heterogeneous_combination_is_symmetric(self):
        low_high = effective_surface_energy(0.01, 0.04)
        high_low = effective_surface_energy(0.04, 0.01)
        self.assertEqual(low_high, 0.02)
        self.assertEqual(low_high, high_low)

    def test_nonpositive_surface_energy_is_rejected(self):
        with self.assertRaises(ValueError):
            effective_surface_energy(0.0, 0.015)


class EffectiveGeometryTests(unittest.TestCase):
    def test_equal_pair_radius_is_half_particle_radius(self):
        self.assertEqual(effective_radius(1.0e-6, 1.0e-6), 0.5e-6)

    def test_wall_radius_is_particle_radius(self):
        self.assertEqual(effective_radius(1.0e-6, None), 1.0e-6)


class ContactParameterTests(unittest.TestCase):
    def test_wall_parameters_match_mfix_surface_energy_contract(self):
        parameters = wall_parameters()
        self.assertAlmostEqual(parameters.effective_surface_energy, 0.015)
        self.assertAlmostEqual(parameters.work_of_adhesion, 0.030)
        self.assertAlmostEqual(
            parameters.pull_off_force_magnitude,
            3.0 * math.pi * 1.0e-6 * 0.015,
        )

    def test_equilibrium_radius_has_zero_elastic_force(self):
        parameters = wall_parameters()
        force = elastic_normal_force(
            parameters.equilibrium_contact_radius,
            parameters,
        )
        self.assertAlmostEqual(force, 0.0, delta=1.0e-24)

    def test_critical_separation_is_end_of_stable_branch(self):
        parameters = wall_parameters()
        radius = stable_contact_radius(
            -parameters.critical_separation,
            parameters,
        )
        expected_ratio = (1.0 / 6.0) ** (2.0 / 3.0)
        self.assertAlmostEqual(
            radius / parameters.equilibrium_contact_radius,
            expected_ratio,
            delta=1.0e-12,
        )

    def test_beyond_critical_separation_is_rejected(self):
        parameters = wall_parameters()
        with self.assertRaises(ValueError):
            stable_contact_radius(
                -1.01 * parameters.critical_separation,
                parameters,
            )


class ContactStateTests(unittest.TestCase):
    def test_nonnegative_overlap_establishes_contact(self):
        self.assertTrue(contact_is_active(0.0, False, 1.0e-9))

    def test_tensile_overlap_requires_history(self):
        self.assertTrue(contact_is_active(-0.5e-9, True, 1.0e-9))
        self.assertFalse(contact_is_active(-0.5e-9, False, 1.0e-9))

    def test_contact_breaks_at_excessive_gap(self):
        self.assertFalse(contact_is_active(-1.1e-9, True, 1.0e-9))


class DampingAndCriticalSpeedTests(unittest.TestCase):
    def test_unit_restitution_has_zero_damping(self):
        self.assertEqual(damping_alpha(1.0), 0.0)

    def test_marshall_fit_is_preserved(self):
        self.assertAlmostEqual(damping_alpha(0.6), 0.2532703808)

    def test_wall_critical_speed_matches_audited_legacy_value(self):
        speed = adhesion_only_critical_speed(wall_parameters())
        self.assertAlmostEqual(speed, 0.32816731034019575)


if __name__ == "__main__":
    unittest.main()

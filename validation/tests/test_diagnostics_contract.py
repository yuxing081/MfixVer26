"""Static contract tests for the opt-in double-precision JKR diagnostics."""

from pathlib import Path
import re
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
UNPACK_SOURCE = PROJECT_ROOT / "code" / "mpi_unpack_des_mod.f"
NEIGHBOR_SOURCE = PROJECT_ROOT / "code" / "neighbour.f"


class JkrDiagnosticContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.unpack = UNPACK_SOURCE.read_text(encoding="utf-8")
        cls.neighbor = NEIGHBOR_SOURCE.read_text(encoding="utf-8")

    def test_diagnostics_are_default_off_and_environment_opt_in(self):
        self.assertRegex(
            self.unpack,
            r"LOGICAL\s*::\s*JKR_DIAG_CHECKED\s*=\s*\.FALSE\.",
        )
        self.assertRegex(
            self.unpack,
            r"LOGICAL\s*::\s*JKR_DIAG_ENABLED\s*=\s*\.FALSE\.",
        )
        self.assertRegex(
            self.unpack,
            r"INTEGER\s*::\s*JKR_DIAG_UNIT\s*=\s*-1",
        )
        self.assertIn(
            "GET_ENVIRONMENT_VARIABLE('MFIX_JKR_DIAGNOSTICS'",
            self.unpack,
        )
        self.assertRegex(
            self.unpack,
            r"JKR_DIAG_ENABLED\s*=\s*ENV_STATUS\s*==\s*0",
        )
        self.assertIn("IF(.NOT.JKR_DIAG_ENABLED) RETURN", self.unpack)

    def test_header_contains_all_required_double_precision_fields(self):
        required = (
            "time",
            "phase",
            "rank",
            "owner_gid",
            "neigh_gid",
            "wall_id",
            "local_slot",
            "overlap",
            "normal_force",
            "tangential_force",
            "pft",
            "pfr",
            "jkr_state",
        )
        header_match = re.search(
            r"WRITE\(JKR_DIAG_UNIT,'\(A\)'\).*?"
            r"('time phase.*?jkr_state')",
            self.unpack,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(header_match)
        header = header_match.group(1).lower()
        for field in required:
            with self.subTest(field=field):
                self.assertIn(field, header)

    def test_diagnostic_block_does_not_reuse_des_usr_var(self):
        block = self.unpack.split(
            "! Default-off, rank-local double-precision JKR diagnostic stream.",
            1,
        )[1]
        self.assertNotIn("DES_USR_VAR", block.upper())
        self.assertIn("TIME", block)
        self.assertIn("MYPE", block)
        self.assertIn("OWNER_GID", block)
        self.assertIn("NEIGH_GID", block)

    def test_neighbor_stage_labels_cover_rebuild_and_restore_boundaries(self):
        for phase in (
            "neighbor_before_rebuild",
            "pending_restore_direct",
            "pending_restore_reverse",
            "neighbor_after_restore",
        ):
            with self.subTest(phase=phase):
                self.assertIn(phase, self.unpack + self.neighbor)

    def test_diagnostic_calls_have_no_case_specific_fixed_identity(self):
        calls = re.findall(
            r"CALL\s+JKR_DIAG_PAIR_EVENT\((.*?)\)",
            self.unpack + self.neighbor,
            flags=re.IGNORECASE | re.DOTALL,
        )
        self.assertGreaterEqual(len(calls), 4)
        source = (self.unpack + self.neighbor).lower()
        for forbidden in (
            "pair_cross_jkr",
            "pair_new_contact_jkr",
            "restart_pair_pfr",
            "restart_wall_oblique",
        ):
            self.assertNotIn(forbidden, source)
        for call in calls:
            self.assertNotIn("MFIX_JKR_DIAGNOSTICS", call)
        self.assertIn("TIME,TRIM(PHASE_TAG),MYPE,OWNER_GID", self.unpack)


if __name__ == "__main__":
    unittest.main()

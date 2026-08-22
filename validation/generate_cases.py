#!/usr/bin/env python3
"""Regenerate deterministic, Git-tracked MFIX validation case assets."""

from validation.jkr_tools.case_generation import write_canonical_case_asset
from validation.jkr_tools.config import default_numerical_config, material_config
from validation.run_validation import CASES_ROOT, canonical_contact_cases


def main() -> int:
    numerical = default_numerical_config()
    for case in canonical_contact_cases():
        output = write_canonical_case_asset(
            CASES_ROOT,
            case,
            material_config(case.material_profile),
            numerical,
        )
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

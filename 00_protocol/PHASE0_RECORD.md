# Phase 0 Freeze Record

- **Spec file:** 00_protocol/spec.md
- **Git commit hash:** 3c02d34194f45f2f5a8df24492deaee6ead98a80
- **Commit timestamp (UTC):** 2026-08-29 17:08:03 +0000

This record is kept separate from spec.md itself so that freezing the spec
does not require the file to reference its own hash. To verify the spec has
not been altered since freezing, run:

    git show 3c02d34194f45f2f5a8df24492deaee6ead98a80:00_protocol/spec.md

Any change to spec.md after this point must be made as a new, separately
dated commit — never as an edit that overwrites this history — so the full
decision trail (including what was believed before Phase 1 data was seen)
remains auditable.

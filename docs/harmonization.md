# Pipeline harmonization notes

This repository is being aligned with the CGU Hydra pipeline layout described in
the bootstrap/harmonization plan.

## Added in this pass

- `containers/manifest.yaml`
- `profiles/local/config.yaml`
- `workflow/manifest.yaml`

Poirot production overlays currently live in the companion `poirot_config`
repository. The pipeline repository should stay focused on reusable workflow
structure and should not need cluster-expanded reference maps once `cgu/library`
is stable.

## Current workflow target

The main workflow entrypoint is `workflow/Snakefile`. Its local rules, hydra
modules, and script candidates are described in `workflow/manifest.yaml`.

## Remaining refactor candidates

- Decide whether production config should remain in `poirot_config` or be
  packed as a small versioned config bundle beside the workflow.
- Convert container references in `config/config.yaml` to `{{CONTAINER_DIR}}/*.sif`
  after cluster validation of the manifest.
- Keep `workflow/Snakefile` thin and split new local logic into
  `workflow/rules/*.smk`.
- Compare `workflow/scripts/create_peddy_fam.py`,
  `workflow/scripts/create_peddy_mqc_config.py`, and VCF helper scripts with
  Hastings/Marple for possible migration to `hydra-genetics/report`.

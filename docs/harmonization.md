# Pipeline harmonization notes

This repository is being aligned with the CGU Hydra pipeline layout described in
the bootstrap/harmonization plan.

## Added in this pass

- `containers/manifest.yaml`
- `profiles/local/config.yaml`

Poirot production reference/site overlays currently live in the companion
`poirot_config` repository. That repository owns the new `config/site/*.local.yaml`
files for this harmonization pass.

## Current standard invocation target

With the companion config repository staged into the run directory:

```bash
snakemake \
  --profile profiles/slurm \
  --configfile config/config.yaml \
  --configfile config/site/marvin.local.yaml
```

## Remaining refactor candidates

- Decide whether production config should remain in `poirot_config` or be
  vendored into this repo under `config/site/`.
- Convert container references in `config/config.yaml` to `{{CONTAINER_DIR}}/*.sif`
  after cluster validation of the manifest.
- Keep `workflow/Snakefile` thin and split new local logic into
  `workflow/rules/*.smk`.
- Compare `workflow/scripts/create_peddy_fam.py`,
  `workflow/scripts/create_peddy_mqc_config.py`, and VCF helper scripts with
  Hastings/Marple for possible migration to `hydra-genetics/report`.

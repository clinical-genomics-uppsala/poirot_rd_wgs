#!/usr/bin/env python
# coding: utf-8

"""
Convert samples file from trioid/trio_member format to trio/father/mother format
for somalier_trio compatibility. Auto-detects format and only converts if needed.
"""

import pandas as pd
import sys
from pathlib import Path

# Add scripts directory to path to import somalier_utils
sys.path.insert(0, str(Path(snakemake.scriptdir).resolve()))
from somalier_utils import get_samples_for_somalier  # noqa: E402


# Load samples
samples_df = pd.read_table(snakemake.input.samples, dtype=str).set_index("sample", drop=False)

# Convert using the adapter function
samples_converted = get_samples_for_somalier(samples_df)

# Write converted samples
samples_converted.to_csv(snakemake.output.samples, sep="\t", index=False)

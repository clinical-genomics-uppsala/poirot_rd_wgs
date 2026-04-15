#!/usr/bin/env python
# coding: utf-8

"""
Utility functions for somalier trio format conversion.
Shared between workflow scripts and common.smk.
"""

import pandas as pd


def convert_trio_format_to_somalier(samples_dict):
    """
    Adapter function to convert trioid/trio_member format to trio/father/mother format.

    Input format (current pipeline):
        sample  trioid      trio_member
        child1  TRIO001     proband
        dad1    TRIO001     father
        mom1    TRIO001     mother

    Output format (somalier_trio expected):
        sample  trio        father  mother
        child1  TRIO001     dad1    mom1
        dad1    .           .       .
        mom1    .           .       .

    Returns a copy of the DataFrame with added trio/father/mother columns.
    """
    df = samples_dict.copy()

    # Initialize new columns with default values
    df["trio"] = "."
    df["father"] = "."
    df["mother"] = "."

    # Get unique trio IDs (excluding NA, nan, ., 0)
    trio_ids = df[
        df["trioid"].notna() & (df["trioid"] != "NA") &
        (df["trioid"] != ".") & (df["trioid"] != "0")
    ]["trioid"].unique()

    for trio_id in trio_ids:
        trio_samples = df[df["trioid"] == trio_id]

        # Find proband, father, mother samples
        proband_samples = trio_samples[trio_samples["trio_member"] == "proband"]
        father_samples = trio_samples[trio_samples["trio_member"] == "father"]
        mother_samples = trio_samples[trio_samples["trio_member"] == "mother"]

        if len(proband_samples) > 0 and len(father_samples) > 0 and len(mother_samples) > 0:
            proband_name = proband_samples.index[0]
            father_name = father_samples.index[0]
            mother_name = mother_samples.index[0]

            # Set trio information for proband
            df.at[proband_name, "trio"] = trio_id
            df.at[proband_name, "father"] = father_name
            df.at[proband_name, "mother"] = mother_name

    return df


def get_samples_for_somalier(samples_dict):
    """
    Get samples DataFrame in the format expected by somalier_trio.
    Auto-detects format and converts if needed.
    """
    # Check which format we have
    has_somalier_format = all(col in samples_dict.columns for col in ["trio", "father", "mother"])
    has_pipeline_format = all(col in samples_dict.columns for col in ["trioid", "trio_member"])

    if has_somalier_format:
        # Already in correct format
        return samples_dict
    elif has_pipeline_format:
        # Convert from pipeline format to somalier format
        return convert_trio_format_to_somalier(samples_dict)
    else:
        # No trio information available
        return samples_dict


def has_trio_samples(samples_dict):
    """
    Check if samples file has trio information.
    Works with both trioid/trio_member and trio/father/mother formats.
    Uses adapter to convert trioid/trio_member to trio/father/mother if needed.
    """
    samples_for_somalier = get_samples_for_somalier(samples_dict)
    required_cols = ["trio", "father", "mother"]

    # Check if required columns exist in the DataFrame
    if not all(col in samples_for_somalier.columns for col in required_cols):
        return False

    # Check if any sample has non-empty trio information
    return any(
        (samples_for_somalier["trio"].notna()) &
        (samples_for_somalier["trio"] != ".") &
        (samples_for_somalier["trio"] != "0")
    )

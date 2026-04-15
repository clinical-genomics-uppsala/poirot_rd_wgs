#!/usr/bin/env python
# coding: utf-8

"""
Create MultiQC custom content TSV files from somalier output.
Similar to create_peddy_mqc_config.py but for somalier data.
"""

import sys
import pandas as pd
import yaml
import numpy as np
import traceback
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Create Somalier MultiQC config files.')
    parser.add_argument('--pairs', required=True, help='Somalier pairs.tsv input')
    parser.add_argument('--samples', required=True, help='Somalier samples.tsv input')
    parser.add_argument('--ped', required=True, help='PED file input')
    parser.add_argument('--config', required=True, help='MultiQC config YAML input')
    parser.add_argument('--rel-check-mqc', required=True, help='Output relatedness TSV')
    parser.add_argument('--sex-check-mqc', required=True, help='Output sex check TSV')
    parser.add_argument('--general-stats-mqc', required=True, help='Output general stats TSV')
    return parser.parse_args()


def comment_the_config_keys(config_dict):
    """
    Convert config dictionary to commented YAML for embedding in TSV.
    """
    commented_config = '\n'.join(
        ['# ' + line for line in yaml.dump(config_dict).rstrip('\n').split('\n')])

    return commented_config


def get_trio_info(ped_filepath):
    """
    Extract trio membership information from PED file.
    Returns dict mapping sample_id to family_id.
    """
    ped_cols = ['family_id', 'sample_id', 'paternal_id',
                'maternal_id', 'sex', 'phenotype']

    ped_df = pd.read_csv(ped_filepath, sep=None, engine='python', header=None, names=ped_cols, dtype=str)

    # Clean column names
    ped_df.columns = ped_df.columns.str.lstrip('#').str.strip()

    trio_membership_dict = {}
    for row in ped_df.itertuples():
        trio_membership_dict[row.sample_id] = row.family_id

    return trio_membership_dict


def get_relatedness_df(pairs_file_path, samples_file_path, ped_df, trio_dict):
    """
    Process somalier pairs.tsv file and add QC checks similar to Peddy.
    """
    # Read somalier pairs file
    relatedness_df = pd.read_csv(pairs_file_path, sep='\t', dtype=str)

    # Clean column names
    relatedness_df.columns = relatedness_df.columns.str.lstrip('#').str.strip()

    # Read somalier samples file for additional info
    samples_df = pd.read_csv(samples_file_path, sep='\t', dtype=str)
    samples_df.columns = samples_df.columns.str.lstrip('#').str.strip()

    # Convert relatedness to float for comparison
    relatedness_df['relatedness'] = pd.to_numeric(relatedness_df['relatedness'], errors='coerce')

    # Calculate expected relatedness for each pair using PED file
    relatedness_df['expected_relatedness'] = relatedness_df.apply(
        lambda row: get_expected_relatedness(row['sample_a'], row['sample_b'], ped_df),
        axis=1
    )

    # Determine if the relationship passes QC
    # Pass if observed relatedness matches expected (within tolerance)
    def check_relatedness(row):
        obs = row['relatedness']
        exp = row['expected_relatedness']

        if pd.isna(obs):
            return 'Fail'

        # For parent-child or siblings (expected = 0.5), allow 0.35-0.65
        if exp == 0.5:
            return 'Pass' if 0.35 <= obs <= 0.65 else 'Fail'
        # For unrelated (expected = 0.0), allow up to 0.1
        elif exp == 0.0:
            return 'Pass' if obs < 0.1 else 'Fail'
        # For duplicates (expected = 1.0), allow > 0.9
        elif exp == 1.0:
            return 'Pass' if obs > 0.9 else 'Fail'
        else:
            return 'Fail'

    relatedness_df['rel_check_test'] = relatedness_df.apply(check_relatedness, axis=1)

    # Add trio_id column for each sample pair
    # For trio display, we want to group samples that belong to the same trio
    # This includes parent-child pairs and sibling pairs
    def get_trio_id_for_pair(row):
        sample_a = row['sample_a']
        sample_b = row['sample_b']

        # If they have expected relatedness > 0, they're in a trio together
        if row['expected_relatedness'] > 0:
            # Use the child's family_id if one is a parent-child relationship
            info_a = ped_df[ped_df['sample_id'] == sample_a]
            info_b = ped_df[ped_df['sample_id'] == sample_b]

            if not info_a.empty and not info_b.empty:
                info_a = info_a.iloc[0]
                info_b = info_b.iloc[0]

                # If sample_b is child of sample_a, use sample_b's family
                if sample_a == info_b['paternal_id'] or sample_a == info_b['maternal_id']:
                    return info_b['family_id']
                # If sample_a is child of sample_b, use sample_a's family
                elif sample_b == info_a['paternal_id'] or sample_b == info_a['maternal_id']:
                    return info_a['family_id']
                # If they're siblings, use either family (they should be the same)
                else:
                    return info_a['family_id']

        return 'NA'

    relatedness_df['trio_id'] = relatedness_df.apply(get_trio_id_for_pair, axis=1)

    # Filter to show trio pairs (expected_relatedness > 0) or pairs with QC failures
    trio_or_error_df = relatedness_df[
        (relatedness_df['expected_relatedness'] > 0) |
        (relatedness_df['rel_check_test'] == 'Fail')
    ]

    return trio_or_error_df


def get_sex_check_df(samples_file_path):
    """
    Process somalier samples.tsv file and prepare for MultiQC display.
    """
    sex_check_df = pd.read_csv(samples_file_path, sep=None, engine='python', dtype=str)

    # Clean column names
    sex_check_df.columns = sex_check_df.columns.str.lstrip('#').str.strip()

    # Create predicted_sex column from the sex column (1=male, 2=female)
    sex_check_df['predicted_sex'] = sex_check_df['sex'].replace({
        '1': 'male',
        '2': 'female',
        '0': 'unknown'
    })

    # Somalier already has 'sex_check' column (PASS/FAIL) in some outputs
    # If not, we need to create it
    if 'sex_check' not in sex_check_df.columns:
        # Create sex_check based on predicted_sex vs original_pedigree_sex
        def check_sex(row):
            ped_sex = str(row.get('original_pedigree_sex', 'unknown')).lower()
            pred_sex = str(row.get('predicted_sex', 'unknown')).lower()

            if ped_sex == 'unknown' or ped_sex == '0' or pred_sex == 'unknown':
                return 'UNKNOWN'

            return 'PASS' if ped_sex == pred_sex else 'FAIL'

        sex_check_df['sex_check'] = sex_check_df.apply(check_sex, axis=1)

    return sex_check_df


def write_somalier_mqc(somalier_df, somalier_config, outfile):
    """
    Write MultiQC custom content TSV with embedded config.
    """
    with open(outfile, 'w') as outfile_handle:
        print(comment_the_config_keys(somalier_config), file=outfile_handle)
        somalier_df.to_csv(outfile_handle, sep='\t', mode='a', index=False)


def get_trio_id(sample_id, trio_dict):
    """
    Get trio ID for a sample, return 'NA' if sample is the trio itself.
    """
    if sample_id not in trio_dict:
        return 'NA'

    trio_id = trio_dict[sample_id]
    if trio_id == sample_id:
        trio_id = 'NA'
    return trio_id


def get_expected_relatedness(sample_a, sample_b, ped_df):
    """
    Calculate expected relatedness based on pedigree.
    1.0 = identical/duplicate
    0.5 = parent-child or full siblings
    0.25 = half-siblings
    0.0 = unrelated
    """
    # Get pedigree info for both samples
    info_a = ped_df[ped_df['sample_id'] == sample_a]
    info_b = ped_df[ped_df['sample_id'] == sample_b]

    if info_a.empty or info_b.empty:
        return 0.0

    info_a = info_a.iloc[0]
    info_b = info_b.iloc[0]

    # Same sample = duplicate
    if sample_a == sample_b:
        return 1.0

    # Check if parent-child (check this BEFORE family ID comparison)
    if (sample_a == info_b['paternal_id'] or sample_a == info_b['maternal_id'] or
            sample_b == info_a['paternal_id'] or sample_b == info_a['maternal_id']):
        return 0.5

    # Check if siblings (same parents)
    if (info_a['paternal_id'] == info_b['paternal_id'] and
            info_a['maternal_id'] == info_b['maternal_id'] and
            info_a['paternal_id'] != '0' and info_a['maternal_id'] != '0'):
        return 0.5

    # Different families = unrelated (check this LAST)
    if info_a['family_id'] != info_b['family_id']:
        return 0.0

    # Default to unrelated
    return 0.0


try:
    args = parse_args()

    # Load somalier MultiQC config
    with open(args.config, 'r') as report_configs:
        somalier_mqc_configs = yaml.load(report_configs, Loader=yaml.FullLoader)

    # Get trio info from ped file
    trio_dict = get_trio_info(args.ped)

    # Also load PED file as dataframe for relatedness calculations
    ped_cols = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
    ped_df = pd.read_csv(args.ped, sep=None, engine='python', header=None, names=ped_cols, dtype=str)
    ped_df.columns = ped_df.columns.str.lstrip('#').str.strip()

    # Create relatedness check tsv
    try:
        rel_check_df = get_relatedness_df(
            args.pairs,
            args.samples,
            ped_df,
            trio_dict
        )
    except Exception as e:
        print(f"WARNING: Could not process relatedness data: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        rel_check_df = pd.DataFrame(
            columns=['sample_pair', 'sample_a', 'sample_b', 'relatedness',
                     'expected_relatedness', 'relatedness_check']
        )

    # Add trio_id column if not present
    if 'trio_id' not in rel_check_df.columns and 'sample_a' in rel_check_df.columns:
        rel_check_df['trio_id'] = rel_check_df['sample_a'].apply(
            get_trio_id, args=(trio_dict,)
        )

    if 'trio_id' in rel_check_df.columns:
        rel_check_df.sort_values(by=['trio_id'], inplace=True)

    # Create sample_pair column as first column
    if rel_check_df.shape[0] > 0 and 'sample_a' in rel_check_df.columns and 'sample_b' in rel_check_df.columns:
        rel_check_df['sample_pair'] = rel_check_df['sample_a'].astype(str) + '_v_' + rel_check_df['sample_b'].astype(str)
        first_column = rel_check_df.pop('sample_pair')
        rel_check_df.insert(0, 'sample_pair', first_column)

    if somalier_mqc_configs:
        somalier_rel_config = somalier_mqc_configs.get('somalier_rel_check')
    else:
        somalier_rel_config = {}

    write_somalier_mqc(
        rel_check_df,
        somalier_rel_config,
        args.rel_check_mqc
    )

    # Create sex check tsv
    try:
        sex_check_df = get_sex_check_df(args.samples)
    except Exception as e:
        print(f"WARNING: Could not process sex check data: {e}", file=sys.stderr)
        sex_check_df = pd.DataFrame(columns=['sample_id', 'predicted_sex', 'sex_check'])

    if somalier_mqc_configs:
        somalier_sex_config = somalier_mqc_configs.get('somalier_sex_check')
    else:
        somalier_sex_config = {}

    write_somalier_mqc(
        sex_check_df,
        somalier_sex_config,
        args.sex_check_mqc
    )

    # Create general stats tsv (subset of sex check)
    if 'sample_id' in sex_check_df.columns and 'predicted_sex' in sex_check_df.columns and 'sex_check' in sex_check_df.columns:
        general_stats_df = sex_check_df[['sample_id', 'predicted_sex', 'sex_check']].copy()
        general_stats_df.set_index('sample_id', inplace=True)
    else:
        # Create empty general stats df to satisfy Snakemake requirement
        general_stats_df = pd.DataFrame(columns=['sample_id', 'predicted_sex', 'sex_check'])
        general_stats_df.set_index('sample_id', inplace=True)

    if somalier_mqc_configs:
        somalier_general_config = somalier_mqc_configs.get('somalier_general_stats')
    else:
        somalier_general_config = {}

    write_somalier_mqc(
        general_stats_df,
        somalier_general_config,
        args.general_stats_mqc
    )

except FileNotFoundError as e:
    print(f'File not found: {e}', file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f'Error creating somalier MultiQC files: {e}', file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

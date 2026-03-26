#!/usr/bin/env python
# coding: utf-8

"""
Create MultiQC custom content TSV files from somalier output.
Similar to create_peddy_mqc_config.py but for somalier data.
"""

import pandas as pd
import yaml
import numpy as np
import sys


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
    
    ped_df = pd.read_csv(ped_filepath, sep='\t', header=None, names=ped_cols)
    
    trio_membership_dict = {}
    for row in ped_df.itertuples():
        trio_membership_dict[row.sample_id] = row.family_id
    
    return trio_membership_dict


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
    
    # Different families = unrelated
    if info_a['family_id'] != info_b['family_id']:
        return 0.0
    
    # Check if parent-child
    if (sample_a == info_b['paternal_id'] or sample_a == info_b['maternal_id'] or
        sample_b == info_a['paternal_id'] or sample_b == info_a['maternal_id']):
        return 0.5
    
    # Check if siblings (same parents)
    if (info_a['paternal_id'] == info_b['paternal_id'] and 
        info_a['maternal_id'] == info_b['maternal_id'] and
        info_a['paternal_id'] != '0' and info_a['maternal_id'] != '0'):
        return 0.5
    
    # Default to unrelated
    return 0.0


def get_relatedness_df(pairs_file_path, ped_filepath, trio_dict):
    """
    Process somalier pairs.tsv file and add QC checks.
    """
    # Read ped file for expected relatedness calculation
    ped_cols = ['family_id', 'sample_id', 'paternal_id',
                'maternal_id', 'sex', 'phenotype']
    ped_df = pd.read_csv(ped_filepath, sep='\t', header=None, names=ped_cols)
    
    # Read somalier pairs file
    relatedness_df = pd.read_csv(pairs_file_path, sep='\t')
    
    # Rename 'relatedness' column if needed (somalier uses 'relatedness')
    if 'relatedness' in relatedness_df.columns:
        # Column already named correctly
        pass
    elif 'rel' in relatedness_df.columns:
        relatedness_df.rename(columns={'rel': 'relatedness'}, inplace=True)
    
    # Calculate expected relatedness for each pair
    relatedness_df['expected_relatedness'] = relatedness_df.apply(
        lambda row: get_expected_relatedness(row['sample_a'], row['sample_b'], ped_df),
        axis=1
    )
    
    # Add relatedness check (Pass/Fail)
    # Tolerance of 0.15 for relatedness differences
    tolerance = 0.15
    relatedness_df['relatedness_check'] = np.where(
        abs(relatedness_df['relatedness'] - relatedness_df['expected_relatedness']) <= tolerance,
        'Pass',
        'Fail'
    )
    
    # Filter to show only trio pairs or failures
    trio_idx = []
    non_trio_idx = []
    for row in relatedness_df.itertuples():
        if row.sample_a in trio_dict and row.sample_b in trio_dict:
            if trio_dict[row.sample_a] == trio_dict[row.sample_b]:
                trio_idx.append(row.Index)
            else:
                non_trio_idx.append(row.Index)
        else:
            non_trio_idx.append(row.Index)
    
    if len(trio_idx) > 0:
        trio_rel_df = relatedness_df.iloc[trio_idx, :]
        non_trio_rel_df = relatedness_df.iloc[non_trio_idx, :]
        
        # Get failures from non-trio pairs
        error_rel_df = non_trio_rel_df[non_trio_rel_df['relatedness_check'] == 'Fail']
        
        rel_df = pd.concat([trio_rel_df, error_rel_df])
    else:
        rel_df = relatedness_df[relatedness_df['relatedness_check'] == 'Fail']
    
    return rel_df


def get_sex_check_df(samples_file_path):
    """
    Process somalier samples.tsv file and prepare for MultiQC display.
    """
    sex_check_df = pd.read_csv(samples_file_path, sep='\t')
    
    # Rename sample_id column if needed
    if 'sample_id' not in sex_check_df.columns:
        if '#sample_id' in sex_check_df.columns:
            sex_check_df.rename(columns={'#sample_id': 'sample_id'}, inplace=True)
    
    # Somalier already has 'sex_check' column (PASS/FAIL)
    # If not, we need to create it
    if 'sex_check' not in sex_check_df.columns:
        # Create sex_check based on predicted_sex vs original_pedigree_sex
        def check_sex(row):
            ped_sex = str(row.get('original_pedigree_sex', '0'))
            pred_sex = str(row.get('predicted_sex', 'unknown')).lower()
            
            if ped_sex == '0' or pred_sex == 'unknown':
                return 'UNKNOWN'
            
            pred_code = '1' if pred_sex == 'male' else '2' if pred_sex == 'female' else '0'
            
            return 'PASS' if ped_sex == pred_code else 'FAIL'
        
        sex_check_df['sex_check'] = sex_check_df.apply(check_sex, axis=1)
    
    return sex_check_df


def write_somalier_mqc(somalier_df, somalier_config, outfile):
    """
    Write MultiQC custom content TSV with embedded config.
    """
    with open(outfile, 'w') as outfile:
        print(comment_the_config_keys(somalier_config), file=outfile)
        somalier_df.to_csv(outfile, sep='\t', mode='a', index=False)


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


def main():
    try:
        # Get config path from snakemake config
        config = snakemake.config.get("somalier_trio_mqc", {}).get("mqc_config", "")
        
        if not config:
            sys.exit("somalier_trio_mqc.mqc_config not specified in config.yaml")
        
        # Load somalier MultiQC config
        with open(config, 'r') as report_configs:
            somalier_mqc_configs = yaml.load(report_configs, Loader=yaml.FullLoader)
        
        # Get trio information from PED file
        ped_file = snakemake.input.ped
        trio_dict = get_trio_info(ped_file)
        
        # Process relatedness data
        rel_check_df = get_relatedness_df(
            snakemake.input.pairs,
            ped_file,
            trio_dict
        )
        
        # Add trio_id column
        rel_check_df['trio_id'] = rel_check_df['sample_a'].apply(
            get_trio_id, args=(trio_dict,)
        )
        rel_check_df.sort_values(by=['trio_id'], inplace=True)
        
        # Create sample_pair column as first column
        if rel_check_df.shape[0] > 0:
            rel_check_df['sample_pair'] = rel_check_df[['sample_a', 'sample_b']].agg(
                '_v_'.join, axis=1
            )
            first_column = rel_check_df.pop('sample_pair')
            rel_check_df.insert(0, 'sample_pair', first_column)
        
        # Write relatedness check MultiQC file
        somalier_rel_config = somalier_mqc_configs.get('somalier_rel_check')
        write_somalier_mqc(
            rel_check_df,
            somalier_rel_config,
            snakemake.output.rel_check_mqc
        )
        
        # Process sex check data
        sex_check_df = get_sex_check_df(snakemake.input.samples)
        
        # Write sex check MultiQC file
        somalier_sex_config = somalier_mqc_configs.get('somalier_sex_check')
        write_somalier_mqc(
            sex_check_df,
            somalier_sex_config,
            snakemake.output.sex_check_mqc
        )
        
        # Create general stats file with sex check columns
        # Select only the columns we want in general stats
        general_stats_df = sex_check_df[['sample_id', 'predicted_sex', 'sex_check']].copy()
        general_stats_df.set_index('sample_id', inplace=True)
        
        # Write general stats MultiQC file
        somalier_general_config = somalier_mqc_configs.get('somalier_general_stats')
        write_somalier_mqc(
            general_stats_df,
            somalier_general_config,
            snakemake.output.general_stats_mqc
        )
    
    except FileNotFoundError as e:
        sys.exit(f'File not found: {e}')
    except Exception as e:
        sys.exit(f'Error creating somalier MultiQC files: {e}')


if __name__ == '__main__':
    main()

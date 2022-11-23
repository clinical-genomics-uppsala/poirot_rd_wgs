#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import yaml
import numpy as np
import sys

def comment_the_config_keys(config_dict):

    commented_config = '\n'.join(
        ['# ' + line for line in yaml.dump(config_dict).rstrip('\n').split('\n')])

    return commented_config

def get_trio_info(ped_filepath):

    ped_cols = ['family_id', 'sample_id', 'paternal_id',
    'maternal_id', 'sex', 'phenotype']

    ped_df = pd.read_csv(ped_filepath, sep='\t', header=None, names=ped_cols )

    trio_membership_dict = {}
    for row in ped_df.itertuples():
        trio_membership_dict[row.sample_id] = row.family_id

    return trio_membership_dict


def get_relatedness_df(peddy_rel_file_path, trio_dict):
    relatedness_df = pd.read_csv(peddy_rel_file_path)

    relatedness_df["rel_check_test"] = np.where((relatedness_df.parent_error == True) |
    (relatedness_df.sample_duplication_error == True), 'Fail', 'Pass')

    trio_idx = []
    non_trio_idx = []
    for row in relatedness_df.itertuples():
        if trio_dict[row.sample_a] == trio_dict[row.sample_b]:
            trio_idx.append(row.Index) # get the indices of trio pairs
        else:
            non_trio_idx.append(row.Index)

    if len(trio_idx) > 0: # check if there are any trios
        trio_rel_df = relatedness_df.iloc[trio_idx, :]
        non_trio_rel_df = relatedness_df.iloc[non_trio_idx , :]

        error_rel_df = non_trio_rel_df[(non_trio_rel_df['parent_error'] == True) |
            (non_trio_rel_df['sample_duplication_error'] == True)]

        rel_df = pd.concat([trio_rel_df, error_rel_df ])
    else:
        rel_df = relatedness_df[(relatedness_df['parent_error'] == True) | (relatedness_df['sample_duplication_error'] == True)]

    return rel_df


def get_sex_check_df(sex_check_file_path):
    sex_check_df = pd.read_csv(sex_check_file_path)
    sex_check_df["sex_check_test"] = np.where(sex_check_df.error == False,
    'Pass', 'Fail') # report peddy  error check as pass/fail for simplicity
    sex_check_df.rename(columns={'sample_id':'Sample'}, inplace=True)

    return sex_check_df


def write_peddy_mqc(peddy_df, peddy_config, outfile):

    with open(outfile, 'w') as outfile:
        print(comment_the_config_keys(peddy_config),file=outfile)

        peddy_df.to_csv(outfile, sep='\t', mode='a', index=False)

def get_trio_id(sample_id, trio_dict):

    trio_id = trio_dict[sample_id]
    if trio_id == sample_id:
        trio_id = 'NA'
    return trio_id


def main():

    try:

        config = snakemake.config.get("peddy", '').get("config", '')

        with open(config, 'r') as report_configs:
            peddy_mqc_configs = yaml.load(report_configs, Loader=yaml.FullLoader)

            ped_file = snakemake.input.ped
            trio_dict = get_trio_info(ped_file)

            rel_check_df = get_relatedness_df(snakemake.input.peddy_rel_check, trio_dict)
            rel_check_df['trio_id'] = rel_check_df['sample_a'].apply(get_trio_id, args=(trio_dict,))
            rel_check_df.sort_values(by=['trio_id'], inplace=True)

            # create sample pair column as the first column to be used in  multiqc table
            rel_check_df['sample_pair'] = rel_check_df[['sample_a', 'sample_b']].agg('_v_'.join, axis=1)
            first_column = rel_check_df.pop('sample_pair')
            rel_check_df.insert(0, 'sample_pair', first_column)
            peddy_rel_config = peddy_mqc_configs.get('peddy_rel_check')
            write_peddy_mqc(rel_check_df, peddy_rel_config, snakemake.output.rel_check_mqc)

            sex_check_df = get_sex_check_df(snakemake.input.peddy_sex_check)
            peddy_sex_config = peddy_mqc_configs.get('peddy_sex_check')
            write_peddy_mqc(sex_check_df, peddy_sex_config, snakemake.output.sex_check_mqc)

    except FileNotFoundError:
        sys.exit('Path to peddy config file not found/specified in config.yaml')


if __name__ == '__main__':
    main()

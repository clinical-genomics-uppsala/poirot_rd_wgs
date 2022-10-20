#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import yaml
import numpy as np

def comment_the_config_keys(config_dict):

    commented_config = '\n'.join(
        ['# ' + line for line in yaml.dump(config_dict).rstrip('\n').split('\n')])

    return commented_config


def get_trio_info(ped_filepath):
    '''


    '''

    ped_cols = ['family_id', 'sample_id', 'paternal_id',
    'maternal_id', 'sex', 'phenotype']

    ped_df = pd.read_csv(ped_filepath, sep='\t', header=None, names=ped_cols )

    trio_membership_dict = {}
    for row in ped_df.itertuples():
        trio_membership_dict[row.sample_id] = row.family_id

    return (ped_df, trio_membership_dict)



# def get_relatedness_df(peddy_rel_file_path,ped_df, trio_samples):
#     relatedness_df = pd.read_csv(peddy_rel_file_path)
#
#     relatedness_df["rel_check_test"] = np.where(
#     relatedness_df.pedigree_parents ==  relatedness_df.predicted_parents,
#      'Pass', 'Fail') # report peddy  error check as pass/fail for simplicity
#
#     #subset df to only trio samples
#     trio_relatedness_df = relatedness_df[
#     relatedness_df.sample_a.isin(trio_samples)
#     & relatedness_df.sample_b.isin(trio_samples) ]
#
#     trio_relatedness_with_id_df = pd.merge(
#     ped_df[['family_id', 'sample_id']],trio_relatedness_df,how='right',
#     left_on='sample_id', right_on='sample_a')
#
#     trio_relatedness_with_id_df.drop('sample_id', axis=1, inplace=True)
#
#     trio_relatedness_with_id_df['Sample_pair'] = trio_relatedness_with_id_df[
#     ['sample_a', 'sample_b']].agg('_v_'.join, axis=1)
#
#     trio_relatedness_with_id_df = trio_relatedness_with_id_df[
#     sorted(trio_relatedness_with_id_df.columns)]
#
#     return trio_relatedness_with_id_df

def get_relatedness_df(peddy_rel_file_path,ped_df, trio_membership_dict):
    relatedness_df = pd.read_csv(peddy_rel_file_path)

    relatedness_df["rel_check_test"] = np.where(
    relatedness_df.pedigree_parents ==  relatedness_df.predicted_parents,
     'Pass', 'Fail') # report peddy  error check as pass/fail for simplicity

    #subset df to only trio samples
    relatedness_with_sample_a_id = pd.merge(
    ped_df[['family_id', 'sample_id']],relatedness_df,how='right',
    left_on='sample_id', right_on='sample_a')

    # check that sample_a and sample_b are in the same trio
    trio_pair_indices = []
    for sample_pair in relatedness_with_sample_a_id.itertuples():
        sample_a_trio = trio_membership_dict[sample_pair.sample_a]
        sample_b_trio = trio_membership_dict[sample_pair.sample_b]
        if sample_a_trio == sample_b_trio:
            trio_pair_indices.append(sample_pair.Index)

    trio_pairs_df = relatedness_with_sample_a_id.iloc[trio_pair_indices].sort_values(by='family_id')

    trio_pairs_df.drop('sample_id', axis=1, inplace=True)

    trio_pairs_df['Sample_pair'] = trio_pairs_df[
                                ['sample_a', 'sample_b']].agg('_v_'.join, axis=1)

    trio_pairs_sorted_df = trio_pairs_df[sorted(trio_pairs_df.columns)]

    return trio_pairs_sorted_df


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


def main():

    try:
        config = snakemake.config.get("peddy", '').get("config", '')


        with open(config, 'r') as report_configs:
            peddy_mqc_configs = yaml.load(report_configs, Loader=yaml.FullLoader)

            ped_file = snakemake.input.ped
            ped_df, trio_dict = get_trio_info(ped_file)

            rel_check_df = get_relatedness_df(snakemake.input.peddy_rel_check,
                                                ped_df, trio_dict)
            peddy_rel_config = peddy_mqc_configs.get('peddy_rel_check')
            write_peddy_mqc(rel_check_df, peddy_rel_config, snakemake.output.rel_check_mqc)

            sex_check_df = get_sex_check_df(snakemake.input.peddy_sex_check)
            peddy_sex_config = peddy_mqc_configs.get('peddy_sex_check')
            write_peddy_mqc(sex_check_df, peddy_sex_config, snakemake.output.sex_check_mqc)

    except FileNotFoundError:
        sys.exit('Path to peddy config file not found/specified in config.yaml')


if __name__ == '__main__':
    main()

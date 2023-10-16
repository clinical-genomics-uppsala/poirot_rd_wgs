#!/usr/bin/env python
# coding: utf-8

import pandas as pd


def get_ped_sex(sex):

    if sex == 'male':
        sex_code = '1'
    elif sex == 'female':
        sex_code = '2'
    else:
        sex_code = '0'

    return sex_code


samples = pd.read_table(
    snakemake.input[0], dtype=str).set_index("sample", drop=False)

fam_df = samples[['sample', 'sex', 'trioid', 'trio_member']]

child_df = fam_df[samples.trio_member == 'proband']
father_df = fam_df[samples.trio_member == 'father']
mother_df = fam_df[samples.trio_member == 'mother']

with open(snakemake.output[0], 'w') as pedfile:
    phenotype = '-9'

    for sample in fam_df.itertuples():
        sample_id = sample.sample + '_N'  # vcf sample ids have type
        family_id = sample.trioid
        if pd.isna(family_id):
            family_id = sample_id

        if sample.sample in child_df['sample'].tolist():
            child_trio = sample.trioid
            try:
                paternal_id = father_df[father_df.trioid == child_trio].iat[0, 0] + '_N'
            except IndexError:
                paternal_id = '0'
            try:
                maternal_id = mother_df[mother_df.trioid == child_trio].iat[0, 0] + '_N'
            except IndexError:
                maternal_id = '0'
        else:
            paternal_id = '0'
            maternal_id = '0'

        sex = get_ped_sex(sample.sex)
        print(family_id, sample_id, paternal_id, maternal_id,
              sex, phenotype, sep='\t', file=pedfile)

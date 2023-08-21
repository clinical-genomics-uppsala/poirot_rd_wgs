#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
import os


def translate_sex(sex_code):

    if sex_code == "M":
        sex = "male"
    elif sex_code == "K":
        sex = "female"
    elif sex_code == "O":
        sex = "unknown"
    else:
        print('Sex is not specified correctly in the Sample Sheet')
        sys.exit(1)
    return sex


def extract_trio_info(samples, trio_col):  
    trio_member_list = []
    trio_id_list = []
    sex_list = []
    for i in trio_col:
        col_list = i.split('_')
        sex = translate_sex(col_list[1])
        sex_list.append(sex)
        trio_info = col_list[2]
        if trio_info == "NA":
            trio_member_list.append("NA")
            trio_id_list.append("NA")
        else:
            trioid = trio_info.split("-")[0]
            trio_id_list.append(trioid)
            trio_member = trio_info.split("-")[1]
            if trio_member == "Foralder" and sex == "female":
                trio_member = "mother"
            elif trio_member == "Foralder" and sex == "male":
                trio_member = "father"
            else:
                trio_member = "proband"
            trio_member_list.append(trio_member)

    trio_df = pd.DataFrame(data={"sample": samples, "sex": sex_list,
                                 "trioid": trio_id_list,        
                                 "trio_member": trio_member_list})

    return trio_df


def rename_samples(samples, units, sample_id):

    sample_map = {}
    for s1 in sample_id:
        for s2 in samples.itertuples(): 
            try:
                s1_index = s2.sample.index(s1)
                sample_map[s2.sample] = s2.sample[s1_index:]
            except ValueError:
                continue

    samples_renamed = samples.replace({"sample": sample_map})

    units_renamed = units.replace({"sample": sample_map})

    return (samples_renamed, units_renamed)


def main():

    # Read in SampleSheet.csv (illumina sample sheet file)
    # sys.argv[1] is path to SampleSheet.csv file
    sample_sheet_df = pd.read_csv(
        sys.argv[1], header=13, dtype=str).set_index("Sample_ID", drop=False)

    # read in files create by hydra-genetics create-input-files
    samples = pd.read_table("samples.tsv", dtype=str)
    units = pd.read_table("units.tsv", dtype=str)

    # create sample order and replacement files for multiqc
    sample_sheet_df["Sample Order"] = [
        f"sample_{i:03}" for i in range(1, sample_sheet_df.shape[0]+1)]

    sample_replacemen_df = sample_sheet_df[["Sample_ID","Sample Order"]]
    sample_replacemen_df.to_csv(path_or_buf="config/sample_replacement.tsv",
                                sep="\t", index=False, header=False)

    sample_order_df = sample_sheet_df[
        ["Sample Order", "Sample_ID"]].rename(
        columns={"Sample_ID": "Sample Name"})

    sample_order_df.to_csv(
        path_or_buf="config/sample_order.tsv", sep="\t", index=False)

    # add trio info and sex to samples.tsv
    trio_df = extract_trio_info(sample_sheet_df.Sample_ID,
                                sample_sheet_df.Description)

    #  rename samples in units and samples to sample ids in SampleSheet.csv
    samples_renamed = rename_samples(samples, units, sample_sheet_df.Sample_ID)

    samples_renamed_df = samples_renamed[0]
    merged_df = samples_renamed_df.merge(trio_df, on="sample",
                                         validate="one_to_one")
    merged_df.to_csv(path_or_buf="config/samples.tsv", sep='\t', 
                     index=False)

    # fix sample names in units and samples
    units_renamed_df = samples_renamed[1]
    units_renamed_df.to_csv(path_or_buf="config/units.tsv", sep='\t',
                            index=False)


if __name__ == '__main__':
    main()

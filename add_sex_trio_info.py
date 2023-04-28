#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys


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
            # The split on '-' works with current trio namning, needs changing when we see how
            #  how the nex ped ids will look
            trioid = trio_info.split("-")[0]
            trio_id_list.append(trioid)
            trio_member = trio_info.split("-")[1]
            if trio_member == "Foralder" and sex == "female":
                trio_member = "mother"
            elif trio_member == "Foralder" and sex == "male":
                trio_member = "father"
            elif trio_member == "NA":
                trio_member = "NA"
            else:
                trio_member = "proband"
            trio_member_list.append(trio_member)

    trio_df = pd.DataFrame(data={"sample": samples, "sex": sex_list,
                                 "trioid": trio_id_list,        
                                 "trio_member": trio_member_list})

    return trio_df


def main():
    df = pd.read_csv(sys.argv[1], header=13,
                     dtype=str).set_index("Sample_ID", drop=False)
    trio_df = extract_trio_info(df.Sample_ID, df.Description)
    samples = pd.read_table("samples.tsv", dtype=str)
    merged_df = samples.merge(trio_df, on="sample", validate="one_to_one")
    merged_df.to_csv(path_or_buf="config/samples.tsv", sep='\t', index=False)


if __name__ == '__main__':
    main()


import pandas as pd
import sys
import os
import argparse


def translate_sex(sex_code):
    """
    Translates a sex code ("M", "K", or "O") to its English equivalent ("male", "female", or "unknown").

    Args:
        sex_code (str): The sex code to translate.

    Returns:
        str: The English translation of the sex code.

    Raises:
        SystemExit: If the sex code is invalid.
    """

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


def translate_trio_member(trio_member, sex):
    """
    Translates a trio member descriptor ("Foralder", "Index", etc.) and sex ("male" or "female")
    to its English equivalent ("mother", "father", "proband", or "NA").

    Args:
        trio_member (str): The trio member descriptor to translate.
        sex (str): The sex of the individual ("male" or "female").

    Returns:
        str: The English translation of the trio member descriptor.
    """

    if trio_member == "Foralder" and sex == "female":
        trio_member = "mother"
    elif trio_member == "Foralder" and sex == "male":
        trio_member = "father"
    elif trio_member == "Index":
        trio_member = "proband"
    else:
        trio_member = "NA"

    return trio_member


def get_sample_sheet_order(fastq_path):
    """
    Extract the sample sheet numeric order from a
    fastq file path

    Args:
        fastq_path (str): path to fastq1 from the units file

    Returns:
        int: the numeric order of sample in the sample sheet
    """

    fq1_name = os.path.basename(fastq_path)
    # get the 's#' part of the illumina fastq file name
    sample_order = fq1_name.split('_')[2]
    numeric_order = int(sample_order[1:])

    return numeric_order


def format_sample_order(numeric_order):

    return f"sample_{numeric_order:03}"


def main(samples_file, units_file, order_file, replacement_file):

    try:
        # Attempt to read the samples file
        samples = pd.read_csv(samples_file, sep='\t', dtype=str)
    except FileNotFoundError as e:
        print(f"Error: Samples file not found! {e}")
        sys.exit(1)

    # translate the sex code to english
    samples["sex"] = samples["sex"].apply(translate_sex)

    # split trio col
    try:
        samples[["trioid", "trio_member"]] = samples.trio.str.split(
            "-", expand=True)
    except ValueError:  # manually create cols with NAs when no trio present
        samples["trioid"] = ["NA"] * samples.shape[0]
        samples["trio_member"] = ["NA"] * samples.shape[0]
    
    # get the trio member in english
    samples["trio_member"] = samples.apply(
        lambda x: translate_trio_member(x.trio_member, x.sex), axis=1)

    # select the columns we need
    samples = samples[["sample", "sex", "trioid", "trio_member"]]

    samples.to_csv(path_or_buf="samples_with_info.tsv",  na_rep='NA',
                   sep='\t', index=False)

    # read in units.tsv to get the sample order info
    try:
        units = pd.read_csv(units_file, sep='\t', dtype=str)
    except FileNotFoundError as e:
        print(f"Error: Units file not found! {e}")
        sys.exit(1)

    units["sample_order"] = units.fastq1.apply(get_sample_sheet_order)

    sample_order_df = units[
        ["sample_order", "sample"]].drop_duplicates().sort_values(
            by="sample_order")

    sample_order_df["sample_order"] = sample_order_df["sample_order"].apply(
        format_sample_order)

    sample_order_df = sample_order_df.rename(
        columns={"sample_order": "Sample Order", "sample": "Sample Name"})

    sample_order_df.to_csv(
        path_or_buf=order_file, sep="\t", index=False)

    sample_replacement_df = sample_order_df[["Sample Name", "Sample Order"]]
    sample_replacement_df.to_csv(replacement_file, sep="\t",
                                 index=False, header=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract trio member info, translate sex."
        "Generate sample order and replacement files for multiqc.")
    parser.add_argument('-s', '--samples_file', type=str, 
                        help='Path to the samples file')
    parser.add_argument('-u', '--units_file', type=str, 
                        help='Path to the units file')
    parser.add_argument('-o', '--sample_order', type=str, 
                        help='Path to the units file')
    parser.add_argument('-r', '--sample_replacement', type=str, 
                        help='Path to the units file')
    args = parser.parse_args()
    main(args.samples_file, args.units_file, args.sample_order,
         args.sample_replacement)

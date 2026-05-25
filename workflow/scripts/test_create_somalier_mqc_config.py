#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Julia Höglund"
__copyright__ = "Copyright 2026, Julia Höglund"
__email__ = "julia.hoglund@scilifelab.uu.se"
__license__ = "GPL-3"

import sys
import os
import unittest
from unittest.mock import Mock, patch
import tempfile
import pandas as pd
import yaml
import numpy as np

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = TEST_DIR  # Scripts are in the same directory

# Add script directory to path for importing
sys.path.insert(0, SCRIPT_DIR)

# Import the module directly - main block won't execute thanks to if __name__ guard
from create_somalier_mqc_config import (
    comment_the_config_keys,
    get_trio_info,
    get_relatedness_df,
    get_sex_check_df,
    get_expected_relatedness,
)


class TestCommentTheConfigKeys(unittest.TestCase):
    def test_comment_config(self):
        """Test that config dictionary is properly commented"""
        config_dict = {
            "id": "somalier_sex_check",
            "section_name": "Somalier Sex Check"
        }
        result = comment_the_config_keys(config_dict)

        for line in result.split('\n'):
            self.assertTrue(line.startswith('#'), f"Line not commented: {line}")

        self.assertIn("id:", result)
        self.assertIn("somalier_sex_check", result)
        self.assertIn("section_name:", result)


class TestGetTrioInfo(unittest.TestCase):
    def test_get_trio_info(self):
        """Test extraction of trio information from PED file"""
        ped_file = os.path.join(TEST_DIR, ".tests", "create_somalier_mqc.trio.ped")
        trio_dict = get_trio_info(ped_file)

        # Check trio dict structure (maps sample_id to family_id)
        self.assertIn('Child_N', trio_dict)
        self.assertEqual(trio_dict['Child_N'], 'Trio1')
        self.assertEqual(trio_dict['Mother_N'], 'Trio1')
        self.assertEqual(trio_dict['Father_N'], 'Trio1')
        self.assertEqual(trio_dict['Unrelated_N'], 'Single1')


class TestGetExpectedRelatedness(unittest.TestCase):
    def setUp(self):
        """Set up test PED dataframe"""
        ped_file = os.path.join(TEST_DIR, ".tests", "create_somalier_mqc.trio.ped")
        # Load PED directly as dataframe
        ped_cols = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
        self.ped_df = pd.read_csv(ped_file, sep=None, engine='python', header=None, names=ped_cols, dtype=str)

    def test_parent_child_relatedness(self):
        """Test expected relatedness for parent-child (0.5)"""
        rel = get_expected_relatedness('Mother_N', 'Child_N', self.ped_df)
        self.assertEqual(rel, 0.5)

        rel = get_expected_relatedness('Father_N', 'Child_N', self.ped_df)
        self.assertEqual(rel, 0.5)

    def test_parent_parent_relatedness(self):
        """Test expected relatedness for parent-parent (0.0, not related)"""
        rel = get_expected_relatedness('Mother_N', 'Father_N', self.ped_df)
        self.assertEqual(rel, 0.0)

    def test_unrelated_samples(self):
        """Test expected relatedness for unrelated samples (0.0)"""
        rel = get_expected_relatedness('Unrelated_N', 'Mother_N', self.ped_df)
        self.assertEqual(rel, 0.0)


class TestGetRelatednessDF(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = os.path.join(TEST_DIR, ".tests")
        self.pairs_file = os.path.join(self.test_dir, "create_somalier_mqc.trio.pairs.tsv")
        self.samples_file = os.path.join(self.test_dir, "create_somalier_mqc.trio.samples.tsv")
        self.ped_file = os.path.join(self.test_dir, "create_somalier_mqc.trio.ped")
        self.trio_dict = get_trio_info(self.ped_file)
        # Load PED as dataframe
        ped_cols = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
        self.ped_df = pd.read_csv(self.ped_file, sep=None, engine='python', header=None, names=ped_cols, dtype=str)

    def test_relatedness_df_structure(self):
        """Test that relatedness dataframe has required columns"""
        df = get_relatedness_df(self.pairs_file, self.samples_file, self.ped_df, self.trio_dict)

        # Check for key columns (use actual column names from function)
        required_cols = ['sample_a', 'sample_b', 'trio_id', 'relatedness', 'expected_relatedness',
                         'relationship_status', 'rel_check_test']  # rel_check_test not relationship_status_check
        for col in required_cols:
            self.assertIn(col, df.columns, f"Missing column: {col}")

    def test_parent_child_status(self):
        """Test that parent-child pairs are marked as 'Expected (Family)' with Pass"""
        df = get_relatedness_df(self.pairs_file, self.samples_file, self.ped_df, self.trio_dict)

        # Mother-Child pair (use sample_a and sample_b columns)
        mc_pair = df[((df['sample_a'] == 'Mother_N') & (df['sample_b'] == 'Child_N')) |
                     ((df['sample_a'] == 'Child_N') & (df['sample_b'] == 'Mother_N'))]
        self.assertEqual(len(mc_pair), 1)
        self.assertEqual(mc_pair.iloc[0]['relationship_status'], 'Expected (Family)')
        self.assertEqual(mc_pair.iloc[0]['rel_check_test'], 'Pass')  # 'Pass' not 'PASS'
        self.assertEqual(mc_pair.iloc[0]['trio_id'], 'Trio1')

        # Father-Child pair
        fc_pair = df[((df['sample_a'] == 'Father_N') & (df['sample_b'] == 'Child_N')) |
                     ((df['sample_a'] == 'Child_N') & (df['sample_b'] == 'Father_N'))]
        self.assertEqual(len(fc_pair), 1)
        self.assertEqual(fc_pair.iloc[0]['relationship_status'], 'Expected (Family)')
        self.assertEqual(fc_pair.iloc[0]['rel_check_test'], 'Pass')  # 'Pass' not 'PASS'
        self.assertEqual(fc_pair.iloc[0]['trio_id'], 'Trio1')

    def test_parent_parent_trio_id(self):
        """Test that parent-parent pair gets trio_id via shared children"""
        df = get_relatedness_df(self.pairs_file, self.samples_file, self.ped_df, self.trio_dict)

        # Parent-Parent pair should have trio_id = Trio1 (via shared child)
        pp_pair = df[((df['sample_a'] == 'Mother_N') & (df['sample_b'] == 'Father_N')) |
                     ((df['sample_a'] == 'Father_N') & (df['sample_b'] == 'Mother_N'))]
        self.assertEqual(len(pp_pair), 1)
        self.assertEqual(pp_pair.iloc[0]['trio_id'], 'Trio1')

    def test_unrelated_samples_status(self):
        """Test that unrelated samples are marked appropriately"""
        df = get_relatedness_df(self.pairs_file, self.samples_file, self.ped_df, self.trio_dict)

        # Check that we have at least 3 pairs (the trio)
        self.assertGreaterEqual(len(df), 3)

        # Verify a known unrelated fixture pair is present and classified correctly
        unrelated_pair = df[
            ((df['sample_a'] == 'Unrelated_N') & (df['sample_b'] == 'Mother_N')) |
            ((df['sample_a'] == 'Mother_N') & (df['sample_b'] == 'Unrelated_N'))
        ]
        self.assertEqual(len(unrelated_pair), 1)
        self.assertEqual(unrelated_pair.iloc[0]['relationship_status'], 'No Relatedness Expected')


class TestGetSexCheckDF(unittest.TestCase):
    def test_sex_check_df_structure(self):
        """Test that sex check dataframe has required columns"""
        test_dir = os.path.join(TEST_DIR, ".tests")
        samples_file = os.path.join(test_dir, "create_somalier_mqc.trio.samples.tsv")

        df = get_sex_check_df(samples_file)

        required_cols = ['sample_id', 'predicted_sex', 'original_pedigree_sex', 'sex_check']
        for col in required_cols:
            self.assertIn(col, df.columns, f"Missing column: {col}")

    def test_sex_check_pass(self):
        """Test that matching sex is marked as Pass"""
        test_dir = os.path.join(TEST_DIR, ".tests")
        samples_file = os.path.join(test_dir, "create_somalier_mqc.trio.samples.tsv")

        df = get_sex_check_df(samples_file)

        # All samples in test data should pass (sex matches original_pedigree_sex)
        mother = df[df['sample_id'] == 'Mother_N']
        self.assertEqual(mother.iloc[0]['sex_check'], 'PASS')
        self.assertEqual(mother.iloc[0]['predicted_sex'], 'female')

        father = df[df['sample_id'] == 'Father_N']
        self.assertEqual(father.iloc[0]['sex_check'], 'PASS')
        self.assertEqual(father.iloc[0]['predicted_sex'], 'male')


class TestFullWorkflow(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = os.path.join(TEST_DIR, ".tests")

    def test_full_relatedness_workflow(self):
        """Test the complete relatedness check workflow"""
        pairs_file = os.path.join(self.test_dir, "create_somalier_mqc.trio.pairs.tsv")
        samples_file = os.path.join(self.test_dir, "create_somalier_mqc.trio.samples.tsv")
        ped_file = os.path.join(self.test_dir, "create_somalier_mqc.trio.ped")

        # Process files
        trio_dict = get_trio_info(ped_file)
        ped_cols = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
        ped_df = pd.read_csv(ped_file, sep=None, engine='python', header=None, names=ped_cols, dtype=str)
        rel_df = get_relatedness_df(pairs_file, samples_file, ped_df, trio_dict)

        # Verify all pairs were processed (note: function may filter some pairs)
        self.assertGreaterEqual(len(rel_df), 3)  # At least the trio pairs

        # Verify trio pairs have correct trio_id
        trio_pairs = rel_df[rel_df['trio_id'] == 'Trio1']
        self.assertEqual(len(trio_pairs), 3)  # Mother-Father, Mother-Child, Father-Child

        # Verify status checks are working
        pass_checks = rel_df[rel_df['rel_check_test'] == 'Pass']  # 'Pass' not 'PASS'
        self.assertGreaterEqual(len(pass_checks), 3)  # At least trio pairs should pass

    def test_full_sex_check_workflow(self):
        """Test the complete sex check workflow"""
        samples_file = os.path.join(self.test_dir, "create_somalier_mqc.trio.samples.tsv")

        sex_df = get_sex_check_df(samples_file)

        # Verify all samples were processed
        self.assertEqual(len(sex_df), 4)

        # Verify all pass (test data has matching sex)
        pass_checks = sex_df[sex_df['sex_check'] == 'PASS']
        self.assertEqual(len(pass_checks), 4)


if __name__ == "__main__":
    unittest.main()

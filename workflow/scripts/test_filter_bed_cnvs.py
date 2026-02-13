import pytest
import tempfile
import os
import pysam
from filter_bed_cnvs import filter_vcf_to_gz


class TestFilterBedCnvs:

    @pytest.fixture
    def sample_bed_content(self):
        """Sample BED file content with exclusion regions"""
        return """chr1\t1000\t2000
chr1\t5000\t6000
chr2\t3000\t4000"""

    @pytest.fixture
    def sample_vcf_header(self):
        """Sample VCF header for testing"""
        return """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=INS,Description="Insertion">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

    def create_temp_vcf(self, vcf_content):
        """Helper to create a temporary VCF file"""
        temp_vcf = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
        temp_vcf.write(vcf_content)
        temp_vcf.close()
        return temp_vcf.name

    def create_temp_bed(self, bed_content):
        """Helper to create a temporary BED file"""
        temp_bed = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        temp_bed.write(bed_content)
        temp_bed.close()
        return temp_bed.name

    def test_filter_dup_variants_overlapping_bed_regions(self, sample_bed_content, sample_vcf_header):
        """Test that DUP variants overlapping with BED regions are filtered out"""
        vcf_content = sample_vcf_header + "\n" + """chr1\t1500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=1800
chr1\t5500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=5800
chr2\t3500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=3800"""

        input_vcf = self.create_temp_vcf(vcf_content)
        bed_file = self.create_temp_bed(sample_bed_content)
        output_vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False).name

        try:
            filter_vcf_to_gz(input_vcf, bed_file, output_vcf)

            # Check that output file has no records (all DUPs were filtered out)
            vcf_out = pysam.VariantFile(output_vcf)
            records = list(vcf_out.fetch())
            vcf_out.close()

            assert len(records) == 0, "All DUP variants overlapping BED regions should be filtered out"

        finally:
            os.unlink(input_vcf)
            os.unlink(bed_file)
            os.unlink(output_vcf)

    def test_keep_dup_variants_not_overlapping_bed_regions(self, sample_bed_content, sample_vcf_header):
        """Test that DUP variants not overlapping with BED regions are kept"""
        vcf_content = sample_vcf_header + "\n" + """chr1\t500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=800
chr1\t7000\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=8000
chr3\t1000\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=2000"""

        input_vcf = self.create_temp_vcf(vcf_content)
        bed_file = self.create_temp_bed(sample_bed_content)
        output_vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False).name

        try:
            filter_vcf_to_gz(input_vcf, bed_file, output_vcf)

            # Check that all records are kept
            vcf_out = pysam.VariantFile(output_vcf)
            records = list(vcf_out.fetch())
            vcf_out.close()

            assert len(records) == 3, "DUP variants not overlapping BED regions should be kept"

        finally:
            os.unlink(input_vcf)
            os.unlink(bed_file)
            os.unlink(output_vcf)

    def test_keep_non_dup_variants(self, sample_bed_content, sample_vcf_header):
        """Test that non-DUP variants are always kept regardless of overlap"""
        vcf_content = sample_vcf_header + "\n" + """chr1\t1500\t.\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1800
chr1\t5500\t.\tA\t<INV>\t60\tPASS\tSVTYPE=INV;END=5800
chr2\t3500\t.\tA\t<INS>\t60\tPASS\tSVTYPE=INS;END=3600"""

        input_vcf = self.create_temp_vcf(vcf_content)
        bed_file = self.create_temp_bed(sample_bed_content)
        output_vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False).name

        try:
            filter_vcf_to_gz(input_vcf, bed_file, output_vcf)

            # Check that all non-DUP records are kept
            vcf_out = pysam.VariantFile(output_vcf)
            records = list(vcf_out.fetch())
            vcf_out.close()

            assert len(records) == 3, "Non-DUP variants should always be kept"

        finally:
            os.unlink(input_vcf)
            os.unlink(bed_file)
            os.unlink(output_vcf)

    def test_mixed_variants_filtering(self, sample_bed_content, sample_vcf_header):
        """Test filtering behavior with mixed variant types"""
        vcf_content = sample_vcf_header + "\n" + """chr1\t1500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=1800
chr1\t1600\t.\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1700
chr1\t5500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=5800
chr1\t7000\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=8000
chr2\t3500\t.\tA\t<INS>\t60\tPASS\tSVTYPE=INS;END=3600"""

        input_vcf = self.create_temp_vcf(vcf_content)
        bed_file = self.create_temp_bed(sample_bed_content)
        output_vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False).name

        try:
            filter_vcf_to_gz(input_vcf, bed_file, output_vcf)

            vcf_out = pysam.VariantFile(output_vcf)
            records = list(vcf_out.fetch())
            vcf_out.close()

            # Should keep: DEL at 1600 (non-DUP), DUP at 7000 (no overlap), INS at 3500 (non-DUP)
            # Should filter: DUP at 1500 and 5500 (overlap with BED regions)
            assert len(records) == 3, "Should keep non-DUP variants and non-overlapping DUP variants"

            # Check specific variants that should be kept
            kept_positions = [record.pos for record in records]
            assert 1600 in kept_positions, "DEL variant should be kept"
            assert 7000 in kept_positions, "Non-overlapping DUP should be kept"
            assert 3500 in kept_positions, "INS variant should be kept"

        finally:
            os.unlink(input_vcf)
            os.unlink(bed_file)
            os.unlink(output_vcf)

    def test_empty_bed_file(self, sample_vcf_header):
        """Test behavior with empty BED file (no filtering should occur)"""
        vcf_content = sample_vcf_header + "\n" + """chr1\t1500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=1800
chr1\t5500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=5800"""

        input_vcf = self.create_temp_vcf(vcf_content)
        bed_file = self.create_temp_bed("")  # Empty BED file
        output_vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False).name

        try:
            filter_vcf_to_gz(input_vcf, bed_file, output_vcf)

            vcf_out = pysam.VariantFile(output_vcf)
            records = list(vcf_out.fetch())
            vcf_out.close()

            assert len(records) == 2, "All variants should be kept with empty BED file"

        finally:
            os.unlink(input_vcf)
            os.unlink(bed_file)
            os.unlink(output_vcf)

    def test_bed_file_with_comments(self, sample_vcf_header):
        """Test that BED file comments and track lines are properly ignored"""
        bed_content = """# This is a comment
track name="test" description="test track"
chr1\t1000\t2000
# Another comment
browser position chr1:1000-2000
chr2\t3000\t4000"""

        vcf_content = sample_vcf_header + "\n" + """chr1\t1500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=1800
chr2\t3500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=3800
chr3\t1000\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=2000"""

        input_vcf = self.create_temp_vcf(vcf_content)
        bed_file = self.create_temp_bed(bed_content)
        output_vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False).name

        try:
            filter_vcf_to_gz(input_vcf, bed_file, output_vcf)

            vcf_out = pysam.VariantFile(output_vcf)
            records = list(vcf_out.fetch())
            vcf_out.close()

            # Only chr3 DUP should be kept (no overlap with BED regions)
            assert len(records) == 1, "Only non-overlapping DUP should be kept"
            assert records[0].chrom == "chr3", "chr3 DUP should be kept"

        finally:
            os.unlink(input_vcf)
            os.unlink(bed_file)
            os.unlink(output_vcf)

    def test_edge_case_overlap_boundaries(self, sample_vcf_header):
        """Test edge cases for overlap detection at boundaries"""
        bed_content = "chr1\t1000\t2000"

        # Test variants at exact boundaries
        vcf_content = sample_vcf_header + "\n" + """chr1\t999\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=1000
chr1\t2000\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=2001
chr1\t1000\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=2000
chr1\t1500\t.\tA\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=1800"""

        input_vcf = self.create_temp_vcf(vcf_content)
        bed_file = self.create_temp_bed(bed_content)
        output_vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False).name

        try:
            filter_vcf_to_gz(input_vcf, bed_file, output_vcf)

            vcf_out = pysam.VariantFile(output_vcf)
            records = list(vcf_out.fetch())
            vcf_out.close()

            # Check which variants remain
            kept_positions = [record.pos for record in records]

            # Only the variant ending exactly at boundary should be kept (no overlap)
            # Variant 999-1000: touches boundary but doesn't overlap (END=1000 means up to position 1000)
            assert 999 in kept_positions, "Variant ending at boundary should be kept"

            # These should be filtered (overlap with BED region)
            assert 2000 not in kept_positions, "Variant starting at boundary overlaps and should be filtered"
            assert 1000 not in kept_positions, "Variant starting at boundary overlaps and should be filtered"
            assert 1500 not in kept_positions, "Variant inside region should be filtered"

            # Only one variant should remain
            assert len(kept_positions) == 1, f"Expected 1 variant, got {len(kept_positions)}: {kept_positions}"

        finally:
            os.unlink(input_vcf)
            os.unlink(bed_file)
            os.unlink(output_vcf)


if __name__ == "__main__":
    pytest.main([__file__])

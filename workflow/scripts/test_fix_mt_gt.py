import pytest
from fix_mt_gt import fix_gt_field
import pysam


@pytest.fixture
def test_vcf(tmp_path):
    # Create a temporary VCF file with test data
    vcf_path = tmp_path / "test.vcf"
    header = (
        '##fileformat=VCFv4.2\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##contig=<ID=chrM>\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n'
    )
    content = (
        'chrM\t100\t.\tA\tATTT\t.\tPASS\t.\tGT\t0/./1/.\n'
        'chrM\t100\t.\tATT\tATT\t.\tPASS\t.\tGT\t0/1/.\n'
        'chrM\t200\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\n'
    )

    with open(vcf_path, 'w') as f:
        f.write(header + content)

    return str(vcf_path)


def test_fix_gt_field(test_vcf, tmp_path):
    output_vcf = str(tmp_path / "output.vcf")

    fix_gt_field(test_vcf, output_vcf)

    with pysam.VariantFile(output_vcf) as vcf:
        variants = list(vcf)

        # Test first variant - should be fixed
        assert variants[0].samples['SAMPLE1']['GT'] == (0, 1)

        # Test second variant - should be fixed
        assert variants[1].samples['SAMPLE1']['GT'] == (0, 1)

        # Test third variant - should remain unchanged
        assert variants[2].samples['SAMPLE1']['GT'] == (0, 1)

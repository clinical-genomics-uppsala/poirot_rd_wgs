#!/bin/python3.6
import sys
from pysam import VariantFile

vcf_in = VariantFile(sys.argv[1])  # dosen't matter if bgziped or not. Automatically recognizes

# Add reference_description descriptions to new header
new_header = vcf_in.header
#new_header.add_line("reference="+ sys.argv[2])
new_header.add_line("##reference=" + sys.argv[2])

# start new vcf with the new_header
vcf_out = VariantFile(sys.argv[3], 'w', header=new_header)


for record in vcf_in.fetch():
    vcf_out.write(record)


#import pdb; pdb.set_trace()

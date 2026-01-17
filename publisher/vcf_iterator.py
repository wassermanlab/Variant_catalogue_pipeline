#!/usr/bin/env python
# -*- coding: utf-8 -*-
import vcfpy

# Open file, this will read in the header

file = "../test_case/SNV/SNV_filtered_frequ_only_SNV_annotation_table_merged_22_truncated.vcf"
#file = "vcfs/vcf1.vcf.gz"

reader = vcfpy.Reader.from_path(file)

# Build and print header
header = ["#CHROM", "POS", "REF", "ALT", "INFO"] + reader.header.samples.names
print("\t".join(header))

for record in reader:
    if not record.is_snv():
        continue
    line = [record.CHROM, record.POS, record.REF, record.INFO.get("AF_tot_XX_XY")]
    line += [alt.value for alt in record.ALT]
    line += [call.data.get("GT") or "./." for call in record.calls]
    print("\t".join(map(str, line)))
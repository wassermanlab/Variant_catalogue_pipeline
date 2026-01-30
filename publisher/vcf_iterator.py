#!/usr/bin/env python
# -*- coding: utf-8 -*-
import vcfpy
import gzip
import io


# Open file, this will read in the header

#file = "../test_case/SNV/SNV_filtered_frequ_only_SNV_annotation_table_merged_22_truncated.vcf"
file = '../test_case/crop3.vcf'
file = '../test_case/ben-big.vcf.gz'

vcf_reader = None

def crop_region(vcf_reader, regions, output_path="cropped_region.vcf"):
    """
    Write a new VCF file containing only records from the specified region(s).
    regions: list of (chrom, start, end) tuples
    output_path: output VCF file path (string)
    """
    import collections
    # Accept either a single region or a list of regions
    if not isinstance(regions, list):
        regions = [regions]
    # Build a lookup for fast region matching
    region_dict = collections.defaultdict(list)
    for c, s, e in regions:
        region_dict[str(c)].append((s, e))

    # Defensive: re-open the vcf_reader to avoid exhausted iterator
    # vcf_reader may be exhausted if used before; so we require a fresh reader
    # Instead, pass the file path and open a new reader here
    # If vcf_reader is a vcfpy.Reader, get its source file path
    # We'll try to get the file path from vcf_reader, else error
    vcf_source = getattr(vcf_reader, 'stream', None)
    if hasattr(vcf_source, 'name'):
        vcf_path = vcf_source.name
        # Try to detect gzip
        is_gz = vcf_path.endswith('.gz')
        import gzip, io
        if is_gz:
            with gzip.open(vcf_path, 'rb') as gz:
                with io.TextIOWrapper(gz, encoding='utf-8', errors='replace') as f:
                    reader = vcfpy.Reader(stream=f)
                    _write_cropped(reader, region_dict, output_path)
        else:
            with io.TextIOWrapper(open(vcf_path, 'rb'), encoding='utf-8', errors='replace') as f:
                reader = vcfpy.Reader(stream=f)
                _write_cropped(reader, region_dict, output_path)
    else:
        # fallback: try to use the passed vcf_reader (may be exhausted)
        _write_cropped(vcf_reader, region_dict, output_path)

def _write_cropped(reader, region_dict, output_path):
    writer = vcfpy.Writer.from_path(output_path, reader.header)
    count = 0
    for record in reader:
        # Try both '22' and 'chr22' style matching
        chrom_keys = [str(record.CHROM), str(record.CHROM).replace('chr','') if str(record.CHROM).startswith('chr') else 'chr'+str(record.CHROM)]
        found = False
        for chrom in chrom_keys:
            chrom_regions = region_dict.get(chrom, [])
            for s, e in chrom_regions:
                if s <= record.POS <= e:
                    writer.write_record(record)
                    count += 1
                    found = True
                    break
            if found:
                break
    writer.close()
    if count == 0:
        print("Warning: No records written to", output_path)

def read_vcf(reader):
    header = ["#CHROM", "POS", "REF", "ALT", "INFO"] + vcf_reader.header.samples.names
    print("\t".join(header))

    for record in vcf_reader:
        if not record.is_snv():
            print("Skipping non-SNV record")
            print(record.INFO.get("type"))
    #        exit()
            continue
        line = [record.CHROM, record.POS, record.REF, record.INFO.get("AF_tot_XX_XY")]
        line += [alt.value for alt in record.ALT]
        line += [call.data.get("GT") or "./." for call in record.calls]
        print("\t".join(map(str, line)))
    
vcf_reader = None
if file.endswith('.gz'):
    with gzip.open(file, 'rb') as gz:
        with io.TextIOWrapper(gz, encoding='utf-8', errors='replace') as f:  # or errors='ignore'
            vcf_reader = vcfpy.Reader(stream=f)
#            read_vcf(vcf_reader)
            regions = [("22", 27000075, 29999884), ("X", 2701436, 2996505), ("Y", 2781761, 2974173)]

            crop_region(vcf_reader, regions, output_path="cropped_regions.vcf")
                
else:
    with io.TextIOWrapper(open(file, 'rb'), encoding='utf-8', errors='replace') as f:
        vcf_reader = vcfpy.Reader(stream=f)
 #       read_vcf(vcf_reader)



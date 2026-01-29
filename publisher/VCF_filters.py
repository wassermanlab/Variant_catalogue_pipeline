"""
CallFilter classes for extracting table data from VCF files.

Each CallFilter class corresponds to a table in the database. The class reads
VCF files in its constructor and provides a getTableRows() method that returns
a list of dictionaries representing rows in the table.

VCF files contain "calls" (variant calls), and these filters extract and transform
the relevant data for each specific table (genes, frequencies, transcripts, etc.).
"""

import logging
import gzip
import io
from typing import List, Dict, Any, Optional
from abc import ABC, abstractmethod
import urllib.parse

import vcfpy
import os

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

NA = "." # fallback None value filler
    
class CallFilter(ABC):
    vcf_records: List[vcfpy.Record]
    
    def __init__(self, vcf_file_path: str):
        self.vcf_records = []
        self.csq_fields = []
        self.csq_index_map = {}
        self.severity_map = {}
        
        child_class_name = self.__class__.__name__
        logger.info("booting up %s", child_class_name)
        
        #read severity table file
        severity_table_path = os.path.join(os.path.dirname(__file__), "severities.tsv")
        try:
            with open(severity_table_path, "r") as f:
                for line in f.readlines()[1:]:
                    parts = line.strip().split("\t")
                    if len(parts) == 3:
                        severity, consequence = parts
                        self.severity_map[consequence] = int(severity)
        except FileNotFoundError:
            logger.warning("Severity table file not found: %s", severity_table_path)
        
        self.load_vcf_file(vcf_file_path)
        
    
    def load_vcf_file(self, file, type = "SNV"):
        
        def read_vcf(reader: vcfpy.Reader):
                
            if type == "SNV":
            
                csq = reader.header.get_info_field_info("CSQ")
                csq_elements = csq.description.split("Format: ")[1]
                self.csq_fields = csq_elements.split("|")
                self.csq_index_map = {field: index for index, field in enumerate(self.csq_fields)}
                
                for record in reader:
                    try:
                        self.vcf_records.append(record)
                    except Exception as e:
                        logger.error(f"Error processing record {l} in file {file}: {e}")
                        
            elif type == "MT":
                # ????? TBA
                
                self.csq_index_map = {field: index for index, field in enumerate(self.csq_fields)}

                for record in reader:    
                    self.vcf_records.append(record)
                
        if file.endswith('.gz'):
            with gzip.open(file, 'rb') as gz:
                with io.TextIOWrapper(gz, encoding='utf-8', errors='replace') as f:  # or errors='ignore'
                    vcf_reader = vcfpy.Reader(stream=f)
                    read_vcf(vcf_reader)
                        
        else:
            with io.TextIOWrapper(open(file, 'rb'), encoding='utf-8', errors='replace') as f:
                vcf_reader = vcfpy.Reader(stream=f)
                read_vcf(vcf_reader)
        
    
    def get_csq_values(self, record: vcfpy.Record, field_name: str) -> List[str]:
        """
        Helper method to extract a specific CSQ field value from a VCF record.
        
        Args:
            record: VCF record object
            field_name: Name of the CSQ field to extract
        """
        index = self.csq_index_map.get(field_name)
        values = []
        if index is None:
            return []
        csq_list = record.INFO.get("CSQ", [])
        if not csq_list:
            return []
        for list in csq_list:
            csq_parts = list.split("|")
            if index >= len(csq_parts):
                return []
            else:
                values.append(csq_parts[index])
        return values
    
    def get_info_value(self, record: vcfpy.Record, field_name: str, fallback = None) -> str:
        """
        Helper method to extract a specific INFO field value from a VCF record.
        
        Args:
            record: VCF record object
            field_name: Name of the INFO field to extract
        """
        return record.INFO.get(field_name, fallback)

    @abstractmethod
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract table rows from the loaded VCF data.
        
        Returns:
            List of dictionaries where each dict represents a row in the table.
        """
        pass
    

class GenesCallFilter(CallFilter):
    """
    Generates the 'genes' table.
    
    Data Source: VCF INFO field, VEP annotation CSQ subfield SYMBOL
    
    Processing:
    - Parse the CSQ annotation from VCF INFO field
    - Split by pipe delimiter to extract SYMBOL field
    - Filter out NA/empty values
    - Output unique gene short names
    
    Output Fields: short_name
    """
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract unique gene symbols from VEP CSQ annotations.
        
        Returns:
            List of dicts with structure: {'short_name': str}
        """

        short_names = set()
        for record in self.vcf_records:
            
            for gene_symbol in self.get_csq_values(record, "SYMBOL"):
                if gene_symbol and gene_symbol != "NA":
                    short_names.add(gene_symbol)
        return [{'short_name': name} for name in sorted(short_names)]


class TranscriptsCallFilter(CallFilter):
    """
    Generates the 'transcripts' table.
    
    Data Source: VCF INFO field, VEP annotation CSQ subfields: Feature, SYMBOL, SOURCE, TSL
    
    Processing:
    - Parse CSQ annotation from VCF INFO field
    - Extract: transcript ID (Feature), gene (SYMBOL), source (SOURCE), TSL
    - Recode SOURCE: "Ensembl" → "E", "RefSeq" → "R"
    - Filter out entries without transcript IDs
    - Output unique transcript records
    
    Output Fields: transcript_id, gene, transcript_type, tsl
    
    Note: Only applies to SNV and MT variants (not SV, MEI, STR)
    """
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Associate transcripts with genes from VEP annotations.
        
        Returns:
            List of dicts with structure: 
            {'transcript_id': str, 'gene': str, 'transcript_type': str, 'tsl': str}
        """
        # TODO: Parse CSQ field from loaded VCF records and extract Feature, SYMBOL, SOURCE, TSL
        # TODO: Recode SOURCE: Ensembl->E, RefSeq->R
        # TODO: Filter entries without transcript IDs
        # TODO: Return unique transcript records
        
        transcripts = {};
        
        for record in self.vcf_records:
            feature = self.get_csq_values(record, "Feature")
            
            if feature == "" or feature == "NA":
                logger.warning("skipping transcript with no feature: %s", feature)
                continue
            symbol = self.get_csq_values(record, "SYMBOL")
            source = self.get_csq_values(record, "SOURCE")
            tsl = self.get_csq_values(record, "TSL")
            
            l = len(feature)
            if not (l == len(symbol) == len(source) == len(tsl)):
                logger.warning("mismatched lengths for transcript feature, symbol, source, tsl: %d vs %d vs %d vs %d", 
                               l, len(symbol), len(source), len(tsl))
                continue
            for i in range(l):
                transcript_id = feature[i]
                if transcript_id not in transcripts:
                    transcript_type = source[i]
                    if transcript_type == "Ensembl":
                        transcript_type = "E"
                    elif transcript_type == "RefSeq" or transcript_type == "Refseq":
                        transcript_type = "R"
                    else:
                        transcript_type = NA
                        
                    transcripts[transcript_id] = {
                        'transcript_id': transcript_id,
                        'gene': symbol[i],
                        'transcript_type': transcript_type,
                        'tsl': tsl[i]
                    }
            
#                logger.info("seen this transcript before: %s", feature)
        return list(transcripts.values())

class VariantsCallFilter(CallFilter):
    """
    Generates the 'variants' table (master variant list).
    
    Data Source: VCF ID field (variant identifier: chr_pos_ref_alt format)
    
    Processing:
    - Extract variant ID from VCF ID field
    - Assign variant type: "SNV", "MT", or "SV"
    - Output unique variant identifiers with type
    
    Output Fields: variant_id, var_type
    """
    
    def __init__(self, vcf_file_path: str):
        super().__init__(vcf_file_path)
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        # TODO: Extract variant ID from VCF ID field in loaded records
        # TODO: Assign variant type
        # TODO: Return unique variant records
    
        variants = {}
        for record in self.vcf_records:
            variant_id = record.ID[0]
            filter = ";".join(record.FILTER)
            for csq in record.INFO.get("CSQ", []):
                csq_parts = csq.split("|")
                var_type = self.get_csq_values(record, "VARIANT_CLASS")
                if var_type == []:
                    var_type = self.get_info_value(record, "TYPE")
                filter =  ";".join(record.FILTER)
                if variant_id not in variants:
                    variants[variant_id] = {
                        'variant_id': variant_id, 
                        'var_type': var_type[0] if var_type and var_type != [] else NA,
                        'filter': filter
                    }
                else:
                    continue
        
        return list(variants.values())


class VariantsTranscriptsCallFilter(CallFilter):
    """
    Generates the 'variants_transcripts' table.
    
    Data Source: VCF INFO field, VEP annotation CSQ subfields: Feature, variant ID, HGVSc
    
    Processing:
    - Parse CSQ annotation from VCF INFO field
    - Extract: transcript ID (Feature), variant ID, HGVSc (coding sequence change)
    - Filter out intergenic variants (where transcript is NA)
    - Output unique variant-transcript associations
    
    Output Fields: transcript, variant, hgvsc
    """
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Link variants to transcripts with HGVS coding notation.
        
        Args:
            vcf_records: List of VCF records with INFO field containing CSQ annotation
        
        Returns:
            List of dicts with structure: {'transcript': str, 'variant': str, 'hgvsc': str}
        """
        # TODO: Parse CSQ field and extract Feature, variant ID, HGVSc
        # TODO: Filter intergenic variants
        # TODO: Return unique variant-transcript associations
        
        variantsTranscripts = []
        
        for record in self.vcf_records:
            transcript = self.get_csq_values(record, "Feature")
            variant = record.ID[0]
            hgvsc = self.get_csq_values(record, "HGVSc")
            
            l = len(transcript)
            lh = len(hgvsc)
            
            if l != lh:
                logger.warning("mismatched lengths for transcript and hgvsc: %d vs %d", l, lh)
                exit()
            else:
                for i in range(l):
                    t = transcript[i]
                    h = hgvsc[i]
                    if t == "NA" or t == "":
                        continue
                    variantsTranscripts.append({
                        'transcript': t,
                        'variant': variant,
                        'hgvsc': h
                    })
        return variantsTranscripts

class VariantsAnnotationsCallFilter(CallFilter):
    """
    Generates the 'variants_annotations' table.
    
    Data Source: VCF INFO field, VEP annotation CSQ subfields: HGVSp, SIFT, PolyPhen, Feature, variant ID
    
    Processing:
    - Parse CSQ annotation from VCF INFO field
    - Extract: HGVSp (protein change), SIFT, PolyPhen, transcript ID, variant ID
    - Decode URL-encoded characters in HGVSp (%3D → =)
    - Filter to entries with valid HGVSp (protein-coding variants only)
    - Output unique variant-transcript annotation records
    
    Output Fields: hgvsp, sift, polyphen, transcript, variant
    """
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract protein-level annotations with pathogenicity predictions.
        
        Args:
            vcf_records: List of VCF records with INFO field containing CSQ annotation
        
        Returns:
            List of dicts with structure: 
            {'hgvsp': str, 'sift': str, 'polyphen': str, 'transcript': str, 'variant': str}
        """
        # TODO: Parse CSQ field and extract HGVSp, SIFT, PolyPhen, Feature, variant ID
        # TODO: Decode URL-encoded characters in HGVSp
        # TODO: Filter entries without valid HGVSp
        # TODO: Return unique annotation records
        
        annotations = []
        for record in self.vcf_records:
            variant = record.ID[0]
            hgvsp_list = self.get_csq_values(record, "HGVSp")
            sift_list = self.get_csq_values(record, "SIFT")
            polyphen_list = self.get_csq_values(record, "PolyPhen")
            transcript_list = self.get_csq_values(record, "Feature")
            
            if hgvsp_list and hgvsp_list[0] == "":
                continue    
            l = len(hgvsp_list)
            if not (l == len(sift_list) == len(polyphen_list) == len(transcript_list)):
                # pad sift and polyphen lists with NAs to match transcript list length
                # Pad lists with "NA" to match the length of hgvsp_list
                def pad_list(lst, target_len):
                    return lst + [NA] * (target_len - len(lst))
                sift_list = pad_list(sift_list, l)
                polyphen_list = pad_list(polyphen_list, l)
            for i in range(l):
                hgvsp = urllib.parse.unquote(hgvsp_list[i])
                sift = sift_list[i]
                polyphen = polyphen_list[i]
                transcript = transcript_list[i]
                
                if hgvsp == "" or hgvsp == "NA":
                    continue
                
                annotations.append({
                    'hgvsp': hgvsp,
                    'sift': sift,
                    'polyphen': polyphen,
                    'transcript': transcript,
                    'variant': variant
                })
        return annotations

class VariantsConsequencesCallFilter(CallFilter):
    """
    Generates the 'variants_consequences' table.
    
    Data Source: VCF INFO field, VEP annotation CSQ Consequence field, severity_table.tsv
    
    Processing:
    - Parse CSQ annotation from VCF INFO field
    - Extract Consequence field (may contain multiple "&"-separated consequences)
    - Split compound consequences into separate rows
    - Join with severity_table.tsv to convert consequence terms to numeric severity
    - Filter out intergenic variants
    - Output variant-transcript-severity associations
    
    Output Fields: severity, variant, transcript
    """
    
    def __init__(self, vcf_file_path: str):
        super().__init__(vcf_file_path)
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Map variants to transcripts with numeric severity scores.
        
        Returns:
            List of dicts with structure: {'severity': int, 'variant': str, 'transcript': str}
        """
        # TODO: Parse CSQ Consequence field from loaded VCF records
        # TODO: Split compound consequences (separated by &)
        # TODO: Map consequence terms to severity numbers
        # TODO: Filter intergenic variants
        # TODO: Return variant-transcript-severity associations 
        variant_consequences = []
        for record in self.vcf_records:
            variant = record.ID[0]
            transcript_list = self.get_csq_values(record, "Feature")
            consequence_list = self.get_csq_values(record, "Consequence")
            l = len(transcript_list)
            if l != len(consequence_list):
                logger.warning("mismatched lengths for transcript and consequence: %d vs %d", l, len(consequence_list))
                continue
            for i in range(l):
                transcript = transcript_list[i]
                consequences = consequence_list[i].split("&")
                if transcript == "NA" or transcript == "":
                    continue
                for consequence in consequences:
                    severity = self.severity_map.get(consequence)
                    if severity is not None:
                        variant_consequences.append({
                            'severity': severity,
                            'variant': variant,
                            'transcript': transcript
                        })
        return variant_consequences

class SnvsCallFilter(CallFilter):
    """
    Generates the 'snvs' table (SNV-specific annotations).
    
    Data Source: VCF (CHROM, POS, ID, REF, ALT, QUAL, INFO with VEP CSQ)
    
    Processing:
    - Parse VCF fixed fields and VEP CSQ annotation
    - Extract VARIANT_CLASS, CADD_PHRED, Existing_variation, VAR_SYNONYMS
    - Calculate variant length and CADD interpretation (≤15=Tolerable, >15=Damaging)
    - Calculate max SpliceAI score from DS_AG/AL/DG/DL fields
    - Extract dbSNP ID, ClinVar VCV number
    - Filter to unique variants
    
    Output Fields: variant, type, length, chr, pos, ref, alt, cadd_score, cadd_intr,
                   dbsnp_id, dbsnp_url, ucsc_url, ensembl_url, clinvar_url, gnomad_url,
                   clinvar_vcv, splice_ai
    """
    
    def __init__(self, vcf_file_path: str, assembly: Optional[str] = None):
        """
        Initialize with VCF files and optional assembly version.
        
        Args:
            vcf_file_path: List of VCF file paths
            assembly: Genome assembly ("GRCh37" or "GRCh38"), auto-detected if None
        """
        super().__init__(vcf_file_path)
        self.assembly = assembly
        # TODO: Auto-detect assembly from VCF ##contig headers if not provided
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Generate SNV-specific annotations with scores and browser URLs.
        
        Returns:
            List of dicts with structure: {'variant': str, 'type': str, 'length': int,
            'chr': str, 'pos': int, 'ref': str, 'alt': str, 'cadd_score': float,
            'cadd_intr': str, 'dbsnp_id': str, 'dbsnp_url': str, 'ucsc_url': str,
            'ensembl_url': str, 'clinvar_url': str, 'gnomad_url': str, 
            'clinvar_vcv': str, 'splice_ai': float}
        """
        # TODO: Parse VCF fields and CSQ annotation from loaded records
        # TODO: Calculate variant length based on VARIANT_CLASS
        # TODO: Derive CADD interpretation
        # TODO: Calculate max SpliceAI score
        # TODO: Extract dbSNP and ClinVar IDs
        # TODO: Generate browser URLs based on assembly
        # TODO: Return unique SNV annotations
        
        snvs = {}
        for record in self.vcf_records:
            variant = record.ID[0]
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = record.ALT[0].value  # assuming single ALT allele
            qual = record.QUAL
            info = record.INFO
            
            # Extract CSQ values
            csq_list = self.get_csq_values(record, "VARIANT_CLASS")
            cadd_phred_list = self.get_csq_values(record, "CADD_PHRED")
            existing_variation_list = self.get_csq_values(record, "Existing_variation")
            var_synonyms_list = self.get_csq_values(record, "VAR_SYNONYMS")
            ds_ag_list = self.get_csq_values(record, "DS_AG")
            ds_al_list = self.get_csq_values(record, "DS_AL")
            ds_dg_list = self.get_csq_values(record, "DS_DG")
            ds_dl_list = self.get_csq_values(record, "DS_DL")
            dbsnp_ids = []
            clinvar_vcvs = []
            
            for ev in existing_variation_list:
                if ev.startswith("rs"):
                    dbsnp_ids.append(ev)
                if ev.startswith("VCV"):
                    clinvar_vcvs.append(ev)
            
            variant_class = csq_list[0] if csq_list else NA
            cadd_phred = float(cadd_phred_list[0]) if cadd_phred_list and cadd_phred_list[0] != "" else None
            
            # Determine variant length
            if variant_class == "SNV":
                var_length = 1
            elif variant_class in ["INS", "DEL"]:
                var_length = abs(len(alt) - len(ref))
            else:
                var_length = None
            
            # CADD interpretation
            if cadd_phred is not None:
                cadd_intr = "Damaging" if cadd_phred > 15 else "Tolerable"
            else:
                cadd_intr = NA
            #clinvar_vcv
            
            
            # Max SpliceAI score
            splice_ai_scores = []
            for score_list in [ds_ag_list, ds_al_list, ds_dg_list, ds_dl_list]:
                for score in score_list:
                    try:
                        splice_ai_scores.append(float(score))
                    except ValueError:
                        continue
            max_splice_ai = max(splice_ai_scores) if splice_ai_scores else None
            if snvs.get(variant) is None:
                snvs[variant] = {
                    'variant': variant,
                    'type': variant_class,
                    'length': var_length,
                    'chr': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'cadd_score': cadd_phred,
                    'cadd_intr': cadd_intr,
                    'dbsnp_id': dbsnp_ids[0] if dbsnp_ids else NA,
                    "clinvar_vcv": clinvar_vcvs[0] if clinvar_vcvs else NA,
                    "splice_ai": max_splice_ai if max_splice_ai else NA
                }
            else:
#                logger.info("seen this snv before: %s", variant)
                continue
        return list(snvs.values())  

class GenomicIbvlFrequenciesCallFilter(CallFilter):

    """
    Generates the 'genomic_ibvl_frequencies' table.
    
    Data Source: VCF INFO field, Hail-calculated frequencies: AF_tot_XX_XY, AC_tot_XX_XY,
                 AN_tot_XX_XY, hom_tot_XX_XY
    
    Processing:
    - Parse Hail-added INFO fields (comma-separated: total, XX female, XY male)
    - Split into af_tot/xx/xy, ac_tot/xx/xy, an_tot/xx/xy, hom_tot/xx/xy
    - Include QUAL field from VCF
    - Output unique variant frequency records
    
    Output Fields: variant, af_tot, af_xx, af_xy, ac_tot, ac_xx, ac_xy, an_tot,
                   an_xx, an_xy, hom_tot, hom_xx, hom_xy, quality
    """
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract internal cohort allele frequencies stratified by sex.
        
        Args:
            vcf_records: List of VCF records with Hail-added INFO fields
                        Format: {'INFO': {'AF_tot_XX_XY': '0.001,0.002,0.0005', ...}}
        
        Returns:
            List of dicts with structure: {'variant': str, 'af_tot': float, 'af_xx': float,
            'af_xy': float, 'ac_tot': int, 'ac_xx': int, 'ac_xy': int, 'an_tot': int,
            'an_xx': int, 'an_xy': int, 'hom_tot': int, 'hom_xx': int, 'hom_xy': int,
            'quality': float}
        """
        # TODO: Parse Hail INFO fields
        # TODO: Split comma-separated values (tot, XX, XY)
        # TODO: Extract QUAL field
        # TODO: Return unique frequency records
        
        rows = []
        for record in self.vcf_records:
            variant = record.ID[0]
            qual = record.QUAL
            info = record.INFO

            # Each field is a comma-separated string: total, XX, XY
            def parse_info_field(field):
                values = info.get(field, None)
                if values is None:
                    return [None, None, None]
                return values[0:3]

            if "AF_tot_XX_XY" in info:
                # ibvl format
                af_tot, af_xx, af_xy = parse_info_field("AF_tot_XX_XY")
                ac_tot, ac_xx, ac_xy = parse_info_field("AC_tot_XX_XY")
                an_tot, an_xx, an_xy = parse_info_field("AN_tot_XX_XY")
                hom_tot, hom_xx, hom_xy = parse_info_field("hom_tot_XX_XY")
                
                row = {
                            'variant': variant,
                            'af_tot': validate_get(af_tot, i),
                            'ac_tot': validate_get(ac_tot, i),
                            'an_tot': validate_get(an_tot, i),
                            'hom_tot': validate_get(hom_tot, i),
                            'hemi_tot': validate_get(hemi_tot, i),
                            'af_xx': validate_get(af_xx, i),
                            'af_xy': validate_get(af_xy, i),
                            'ac_xy': validate_get(ac_xy, i),
                            'an_xx': validate_get(an_xx, i),
                            'ac_xx': validate_get(ac_xx, i),
                            'an_xy': validate_get(an_xy, i),
                            'hom_xx': validate_get(hom_xx, i),
                            'hom_xy': validate_get(hom_xy, i),
                            'hemi_xx': validate_get(hemi_xx, i),
                            'hemi_xy': validate_get(hemi_xy, i),
                            'quality': qual
                        }
                
                rows.append(row)
                
            else:
                # variome format
                af_tot = info.get("AF", None)
                af_xx = info.get("AF_XX", None)
                af_xy = info.get("AF_XY", None)
                ac_tot = info.get("AC", None)
                ac_xx = info.get("AC_XX", None)
                ac_xy = info.get("AC_XY", None)
                an_tot = info.get("AN", None)
                an_xx = info.get("AN_XX", None)
                an_xy = info.get("AN_XY", None)
                hom_tot = info.get("AC_Hom", None)
                hom_xx = info.get("AC_Hom_XX", None)
                hom_xy = info.get("AC_Hom_XY", None)
                hemi_tot = info.get("AC_Hemi", None)
                hemi_xx = info.get("AC_Hemi_XX", None)
                hemi_xy = info.get("AC_Hemi_XY", None)

                i = 0
                l = len(af_tot)
                
                for i in range(l):
                    
                    def validate_get(v, index):
                        if v is [] or v is None:
                            return NA
                        # if type of v is not list, return v
                        if not isinstance(v, list):
                            return v
                        if len(v) <= index:
                            return NA
                        val = v[index]
                        if v in [None, ""]:
                            return NA
                        return val
                
                    try:
                        row = {
                            'variant': variant,
                            'af_tot': validate_get(af_tot, i),
                            'ac_tot': validate_get(ac_tot, i),
                            'an_tot': validate_get(an_tot, i),
                            'hom_tot': validate_get(hom_tot, i),
                            'hemi_tot': validate_get(hemi_tot, i),
                            'af_xx': validate_get(af_xx, i),
                            'af_xy': validate_get(af_xy, i),
                            'ac_xy': validate_get(ac_xy, i),
                            'an_xx': validate_get(an_xx, i),
                            'ac_xx': validate_get(ac_xx, i),
                            'an_xy': validate_get(an_xy, i),
                            'hom_xx': validate_get(hom_xx, i),
                            'hom_xy': validate_get(hom_xy, i),
                            'hemi_xx': validate_get(hemi_xx, i),
                            'hemi_xy': validate_get(hemi_xy, i),
                            'quality': qual
                        }                        
                        rows.append(row)
                        
                    except Exception as e:
                        logger.warning("Error parsing frequency fields for variant %s: %s", variant, e)
        return rows

class MtsCallFilter(CallFilter):
    """
    Generates the 'mts' table (mitochondrial variant annotations).
    
    Data Source: VCF (CHROM, POS, ID, REF, ALT, INFO with VEP CSQ)
    
    Processing:
    - Parse VCF fixed fields and VEP CSQ annotation
    - Adjust variant IDs for indels to match gnomAD format (position +1, trimmed)
    - Extract dbSNP ID and ClinVar VCV from VEP annotations
    - Generate URLs: UCSC, MitoMap, gnomAD, dbSNP, ClinVar
    - Filter to unique variants
    
    Output Fields: variant, pos, ref, alt, ucsc_url, mitomap_url, gnomad_url,
                   dbsnp_id, dbsnp_url, clinvar_url, clinvar_vcv
    """
    
    def __init__(self, vcf_file_path: str, 
                 assembly: Optional[str] = None):
        """
        Initialize with VCF files, gnomAD data, and optional assembly.
        
        Args:
            vcf_file_path: List of VCF file paths
            assembly: Genome assembly, auto-detected if None
        """
        super().__init__(vcf_file_path)
        self.assembly = assembly
        self.gnomad_variants = set()
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Generate MT-specific annotations with MT database URLs.
        
        Returns:
            List of dicts with structure: {'variant': str, 'pos': int, 'ref': str,
            'alt': str, 'ucsc_url': str, 'mitomap_url': str, 'gnomad_url': str,
            'dbsnp_id': str, 'dbsnp_url': str, 'clinvar_url': str, 'clinvar_vcv': str}
        """
        # TODO: Parse VCF fields and CSQ annotation from loaded records
        # TODO: Adjust variant IDs for indels (gnomAD format)
        # TODO: Extract dbSNP and ClinVar IDs
        # TODO: Generate MT-specific URLs (MitoMap, etc.)
        # TODO: Return unique MT annotations
        
        mts = {}
        for record in self.vcf_records:
            variant = record.ID[0]
            pos = record.POS
            ref = record.REF
            alt = record.ALT[0].value  # assuming single ALT allele

            # Extract CSQ values
            existing_variation_list = self.get_csq_values(record, "Existing_variation")
            dbsnp_id = NA
            clinvar_vcv = NA
            for ev in existing_variation_list:
                if ev.startswith("rs"):
                    dbsnp_id = ev
                if ev.startswith("VCV"):
                    clinvar_vcv = ev

            # URLs (placeholders, adjust as needed)
            ucsc_url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={self.assembly or 'hg38'}&position=chrM%3A{pos}-{pos}"
            mitomap_url = f"https://www.mitomap.org/foswiki/bin/view/MITOMAP/MutationsCodingControl#{pos}"
            gnomad_url = f"https://gnomad.broadinstitute.org/variant/M-{pos}-{ref}-{alt}?dataset=gnomad_r3"
            dbsnp_url = f"https://www.ncbi.nlm.nih.gov/snp/{dbsnp_id}" if dbsnp_id != NA else ""
            clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_vcv[3:]}" if clinvar_vcv != NA else ""

            if variant not in mts:
                mts[variant] = {
                    'variant': variant,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'ucsc_url': ucsc_url,
                    'mitomap_url': mitomap_url,
                    'gnomad_url': gnomad_url,
                    'dbsnp_id': dbsnp_id,
                    'dbsnp_url': dbsnp_url,
                    'clinvar_url': clinvar_url,
                    'clinvar_vcv': clinvar_vcv
                }
        return list(mts.values())

class MtIbvlFrequenciesCallFilter(CallFilter):
    """
    Generates the 'mt_ibvl_frequencies' table.
    
    Data Source: VCF GT fields with Hail-calculated MT-specific metrics
    
    Processing:
    - Parse Hail-added GT fields (not INFO): AC_hom, AC_het, AF_hom, AF_het, AN,
      max_observed_heteroplasmy, heteroplasmy_histogram
    - Extract histogram values and format as comma-separated string
    - Filter out variants where AN=0
    - Adjust variant IDs for indels (gnomAD format)
    - Output unique MT variant frequencies
    
    Output Fields: variant, an, ac_hom, ac_het, af_hom, af_het, hl_hist, max_hl
    """
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract mitochondrial frequencies with heteroplasmy information.
        
        Args:
            vcf_records: List of VCF records with Hail MT-specific GT fields
                        Format: {'GT_fields': {'AC_hom': ..., 'heteroplasmy_histogram': ...}}
        
        Returns:
            List of dicts with structure: {'variant': str, 'an': int, 'ac_hom': int,
            'ac_het': int, 'af_hom': float, 'af_het': float, 'hl_hist': str, 'max_hl': float}
        """
        # TODO: Parse GT fields (not INFO)
        # TODO: Extract and format heteroplasmy histogram
        # TODO: Filter AN=0 variants
        # TODO: Adjust indel variant IDs
        # TODO: Return unique MT frequency records
        
        rows = []
        for record in self.vcf_records:
            variant = record.ID[0]
            gt_fields = record.INFO  # Placeholder; actual GT fields extraction may differ

            an = gt_fields.get("AN", 0)
            if an == 0:
                continue

            ac_hom = gt_fields.get("AC_hom", 0)
            ac_het = gt_fields.get("AC_het", 0)
            af_hom = gt_fields.get("AF_hom", 0.0)
            af_het = gt_fields.get("AF_het", 0.0)
            max_hl = gt_fields.get("max_observed_heteroplasmy", 0.0)
            hl_histogram = gt_fields.get("heteroplasmy_histogram", [])
            hl_hist = ",".join(map(str, hl_histogram)) if hl_histogram else NA

            row = {
                'variant': variant,
                'an': an,
                'ac_hom': ac_hom,
                'ac_het': ac_het,
                'af_hom': af_hom,
                'af_het': af_het,
                'hl_hist': hl_hist,
                'max_hl': max_hl
            }
            rows.append(row)
        return rows


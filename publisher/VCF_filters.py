"""
CallFilter classes for extracting table data from VCF files.

Each CallFilter class corresponds to a table in the database. The class reads
VCF files in its constructor and provides a getTableRows() method that returns
a list of dictionaries representing rows in the table.

VCF files contain "calls" (variant calls), and these filters extract and transform
the relevant data for each specific table (genes, frequencies, transcripts, etc.).
"""

import logging
from typing import List, Dict, Any, Optional
from abc import ABC, abstractmethod

import vcfpy


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
    
class CallFilter(ABC):
    vcf_records: List[vcfpy.Record]
    
    def __init__(self, vcf_file_paths: List[str]):
        self.vcf_records = []
        self.csq_fields = []
        self.csq_index_map = {}
        
        self.load_vcf_files(vcf_file_paths)
    
    def load_vcf_files(self, vcf_file_paths: List[str]):
        """
        Load and parse VCF files into internal records structure.
        Subclasses can override this if they need custom parsing.
        """
        
        for file in vcf_file_paths:
            reader = vcfpy.Reader.from_path(file)
            csq = reader.header.get_info_field_info("CSQ")
            csq_elements = csq.description.split("Format: ")[1]
            self.csq_fields = csq_elements.split("|")
            self.csq_index_map = {field: index for index, field in enumerate(self.csq_fields)}

            for record in reader:
                if not record.is_snv():
                    continue
                self.vcf_records.append(record)
    
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
        
        print("get genees from", self.csq_fields)
        gene_index = self.csq_fields.index("SYMBOL")

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
        feature_index = self.csq_fields.index("Feature")
        symbol_index = self.csq_fields.index("SYMBOL")
        source_index = self.csq_fields.index("SOURCE")
        tsl_index = self.csq_fields.index("TSL")
        for record in self.vcf_records:
            csq = record.INFO.get("CSQ", [])[0]
            csq_parts = csq.split("|")
            feature = csq_parts[feature_index]
            if feature == "" or feature == "NA":
                continue
            symbol = csq_parts[symbol_index]
            source = csq_parts[source_index]
            tsl = csq_parts[tsl_index]
            transcript = {'transcript_id': feature,
                          'gene': symbol,
                          'transcript_type': 'E' if source == 'Ensembl' else ('R' if source == 'RefSeq' else source),
                          'tsl': tsl}
            if feature not in transcripts:
                transcripts[feature] = transcript
            else:
                continue
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
    
    def __init__(self, vcf_file_paths: List[str]):
        super().__init__(vcf_file_paths)
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        # TODO: Extract variant ID from VCF ID field in loaded records
        # TODO: Assign variant type
        # TODO: Return unique variant records
        
        typeIndex = self.csq_fields.index("VARIANT_CLASS")
    
    
        variants = {}
        for record in self.vcf_records:
            variant_id = record.ID[0]
            filter = ";".join(record.FILTER)
            csq = record.INFO.get("CSQ", [])[0]
            csq_parts = csq.split("|")
            variant_class = csq_parts[typeIndex]
#            filter = record.INFO.split("|")[filterIndex]
            if variant_id not in variants:
                variants[variant_id] = {
                    'variant_id': variant_id, 
                    'var_type': variant_class,
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
        raise NotImplementedError("VariantsAnnotationsCallFilter.getTableRows() not yet implemented")


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
    
    def __init__(self, vcf_file_paths: List[str]):
        super().__init__(vcf_file_paths)
        # TODO: Load severity table in _load_vcf_files or separate method
    
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
        raise NotImplementedError("VariantsConsequencesCallFilter.getTableRows() not yet implemented")


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
    - Generate URLs: dbSNP, UCSC, Ensembl, ClinVar, gnomAD
    - Filter to unique variants
    
    Output Fields: variant, type, length, chr, pos, ref, alt, cadd_score, cadd_intr,
                   dbsnp_id, dbsnp_url, ucsc_url, ensembl_url, clinvar_url, gnomad_url,
                   clinvar_vcv, splice_ai
    """
    
    def __init__(self, vcf_file_paths: List[str], assembly: Optional[str] = None):
        """
        Initialize with VCF files and optional assembly version.
        
        Args:
            vcf_file_paths: List of VCF file paths
            assembly: Genome assembly ("GRCh37" or "GRCh38"), auto-detected if None
        """
        super().__init__(vcf_file_paths)
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
        raise NotImplementedError("SnvsCallFilter.getTableRows() not yet implemented")


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
    
    def __init__(self, vcf_file_paths: List[str], gnomad_file_path: str, 
                 assembly: Optional[str] = None):
        """
        Initialize with VCF files, gnomAD data, and optional assembly.
        
        Args:
            vcf_file_paths: List of VCF file paths
            gnomad_file_path: Path to gnomAD MT TSV file
            assembly: Genome assembly, auto-detected if None
        """
        super().__init__(vcf_file_paths)
        self.gnomad_file_path = gnomad_file_path
        self.assembly = assembly
        self.gnomad_variants = set()
        # TODO: Load gnomAD variant IDs in initialization
    
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
        raise NotImplementedError("MtsCallFilter.getTableRows() not yet implemented")


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
        raise NotImplementedError("GenomicIbvlFrequenciesCallFilter.getTableRows() not yet implemented")


class GenomicGnomadFrequenciesCallFilter(CallFilter):
    """
    Generates the 'genomic_gnomad_frequencies' table.
    
    Data Source: External gnomAD VCF files (pre-processed via gnomad_frequency_table.nf)
    
    Processing:
    - Read pre-processed gnomAD TSV (CHROM, POS, REF, ALT, FILTER, AF, AC, AN, nhomalt)
    - Create variant ID: chr_pos_ref_alt
    - Adjust chromosome labels to match pipeline format
    - Intersect with IBVL variants (keep only cohort variants)
    - Output gnomAD frequencies for matching variants
    
    Output Fields (GRCh37): variant, af_tot, ac_tot, an_tot, hom_tot, FILTER
    Output Fields (GRCh38): +exomes_filters, genomes_filters
    """
    
    def __init__(self, vcf_file_paths: List[str], gnomad_file_path: str,
                 assembly: Optional[str] = None):
        """
        Initialize with VCF files and gnomAD data.
        
        Args:
            vcf_file_paths: List of VCF file paths (to get IBVL variant IDs)
            gnomad_file_path: Path to pre-processed gnomAD TSV
            assembly: Genome assembly, auto-detected if None
        """
        super().__init__(vcf_file_paths)
        self.gnomad_file_path = gnomad_file_path
        self.assembly = assembly
        self.ibvl_variants = set()
        self.gnomad_records = []
        # TODO: Extract IBVL variant IDs and load gnomAD records
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract gnomAD population frequencies for cohort variants.
        
        Returns:
            List of dicts with structure: {'variant': str, 'af_tot': float, 'ac_tot': int,
            'an_tot': int, 'hom_tot': int, 'FILTER': str, 'exomes_filters': str (GRCh38),
            'genomes_filters': str (GRCh38)}
        """
        # TODO: Create variant IDs from gnomAD records
        # TODO: Adjust chromosome labels
        # TODO: Intersect with IBVL variants
        # TODO: Include assembly-specific fields
        # TODO: Return gnomAD frequencies for matching variants
        raise NotImplementedError("GenomicGnomadFrequenciesCallFilter.getTableRows() not yet implemented")


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
        raise NotImplementedError("MtIbvlFrequenciesCallFilter.getTableRows() not yet implemented")


class MtGnomadFrequenciesCallFilter(CallFilter):
    """
    Generates the 'mt_gnomad_frequencies' table.
    
    Data Source: External gnomAD mitochondrial TSV file
    
    Processing:
    - Read gnomAD MT table (chromosome, position, ref, alt, AN, AC_hom, AC_het,
      AF_hom, AF_het, max_observed_heteroplasmy)
    - Create variant ID: chr_pos_ref_alt
    - Adjust IDs for indels
    - Intersect with IBVL MT variants
    - Output gnomAD frequencies for matching variants
    
    Output Fields: variant, an, ac_hom, ac_het, af_hom, af_het, max_hl
    """
    
    def __init__(self, vcf_file_paths: List[str], gnomad_mt_file_path: str):
        """
        Initialize with VCF files and gnomAD MT data.
        
        Args:
            vcf_file_paths: List of VCF file paths (to get IBVL MT variant IDs)
            gnomad_mt_file_path: Path to gnomAD MT TSV file
        """
        super().__init__(vcf_file_paths)
        self.gnomad_mt_file_path = gnomad_mt_file_path
        self.ibvl_variants = set()
        self.gnomad_mt_records = []
        # TODO: Extract IBVL MT variant IDs and load gnomAD MT records
    
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Extract gnomAD MT population frequencies for cohort variants.
        
        Returns:
            List of dicts with structure: {'variant': str, 'an': int, 'ac_hom': int,
            'ac_het': int, 'af_hom': float, 'af_het': float, 'max_hl': float}
        """
        # TODO: Create variant IDs from gnomAD MT records
        # TODO: Adjust indel variant IDs
        # TODO: Intersect with IBVL MT variants
        # TODO: Return gnomAD MT frequencies for matching variants
        raise NotImplementedError("MtGnomadFrequenciesCallFilter.getTableRows() not yet implemented")

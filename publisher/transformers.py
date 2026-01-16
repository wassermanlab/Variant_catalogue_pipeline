"""
Transformer classes for converting VCF data to database-ready structures.

Each transformer class corresponds to a table in the database and implements
a transform() method that takes raw data sources and returns a list of dictionaries
representing rows in the table.
"""

from typing import List, Dict, Any, Optional
from abc import ABC, abstractmethod


class BaseTransformer(ABC):
    """Base class for all transformers."""
    
    @abstractmethod
    def transform(self, *args, **kwargs) -> List[Dict[str, Any]]:
        """
        Transform raw data into a list of dictionaries.
        
        Returns:
            List of dictionaries where each dict represents a row in the table.
        """
        pass


class GenesTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Extract unique gene symbols from VEP CSQ annotations.
        
        Args:
            vcf_records: List of VCF records with INFO field containing CSQ annotation
                        Format: {'INFO': {'CSQ': 'allele|consequence|...|SYMBOL|...'}}
        
        Returns:
            List of dicts with structure: {'short_name': str}
        """
        # TODO: Implement parsing of CSQ field
        # TODO: Extract SYMBOL from pipe-delimited CSQ
        # TODO: Filter NA/empty values
        # TODO: Return unique gene names
        raise NotImplementedError("GenesTransformer.transform() not yet implemented")


class TranscriptsTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Associate transcripts with genes from VEP annotations.
        
        Args:
            vcf_records: List of VCF records with INFO field containing CSQ annotation
        
        Returns:
            List of dicts with structure: 
            {'transcript_id': str, 'gene': str, 'transcript_type': str, 'tsl': str}
        """
        # TODO: Parse CSQ field and extract Feature, SYMBOL, SOURCE, TSL
        # TODO: Recode SOURCE: Ensembl->E, RefSeq->R
        # TODO: Filter entries without transcript IDs
        # TODO: Return unique transcript records
        raise NotImplementedError("TranscriptsTransformer.transform() not yet implemented")


class VariantsTransformer(BaseTransformer):
    """
    Generates the 'variants' table (master variant list).
    
    Data Source: VCF ID field (variant identifier: chr_pos_ref_alt format)
    
    Processing:
    - Extract variant ID from VCF ID field
    - Assign variant type: "SNV", "MT", or "SV"
    - Output unique variant identifiers with type
    
    Output Fields: variant_id, var_type
    """
    
    def transform(self, vcf_records: List[Dict[str, Any]], variant_type: str) -> List[Dict[str, Any]]:
        """
        Create master variant list with type classification.
        
        Args:
            vcf_records: List of VCF records with ID field
                        Format: {'ID': 'chr_pos_ref_alt'}
            variant_type: One of "SNV", "MT", or "SV"
        
        Returns:
            List of dicts with structure: {'variant_id': str, 'var_type': str}
        """
        # TODO: Extract variant ID from VCF ID field
        # TODO: Assign variant type
        # TODO: Return unique variant records
        raise NotImplementedError("VariantsTransformer.transform() not yet implemented")


class VariantsTranscriptsTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
        raise NotImplementedError("VariantsTranscriptsTransformer.transform() not yet implemented")


class VariantsAnnotationsTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
        raise NotImplementedError("VariantsAnnotationsTransformer.transform() not yet implemented")


class VariantsConsequencesTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]], 
                  severity_table: Dict[str, int]) -> List[Dict[str, Any]]:
        """
        Map variants to transcripts with numeric severity scores.
        
        Args:
            vcf_records: List of VCF records with INFO field containing CSQ annotation
            severity_table: Dict mapping consequence terms to numeric severity scores
        
        Returns:
            List of dicts with structure: {'severity': int, 'variant': str, 'transcript': str}
        """
        # TODO: Parse CSQ Consequence field
        # TODO: Split compound consequences (separated by &)
        # TODO: Map consequence terms to severity numbers
        # TODO: Filter intergenic variants
        # TODO: Return variant-transcript-severity associations
        raise NotImplementedError("VariantsConsequencesTransformer.transform() not yet implemented")


class SnvsTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]], 
                  assembly: str) -> List[Dict[str, Any]]:
        """
        Generate SNV-specific annotations with scores and browser URLs.
        
        Args:
            vcf_records: List of VCF records with full annotation
            assembly: Genome assembly ("GRCh37" or "GRCh38")
        
        Returns:
            List of dicts with structure: {'variant': str, 'type': str, 'length': int,
            'chr': str, 'pos': int, 'ref': str, 'alt': str, 'cadd_score': float,
            'cadd_intr': str, 'dbsnp_id': str, 'dbsnp_url': str, 'ucsc_url': str,
            'ensembl_url': str, 'clinvar_url': str, 'gnomad_url': str, 
            'clinvar_vcv': str, 'splice_ai': float}
        """
        # TODO: Parse VCF fields and CSQ annotation
        # TODO: Calculate variant length based on VARIANT_CLASS
        # TODO: Derive CADD interpretation
        # TODO: Calculate max SpliceAI score
        # TODO: Extract dbSNP and ClinVar IDs
        # TODO: Generate browser URLs based on assembly
        # TODO: Return unique SNV annotations
        raise NotImplementedError("SnvsTransformer.transform() not yet implemented")


class MtsTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]], 
                  assembly: str,
                  gnomad_variants: List[str]) -> List[Dict[str, Any]]:
        """
        Generate MT-specific annotations with MT database URLs.
        
        Args:
            vcf_records: List of VCF records with full annotation
            assembly: Genome assembly
            gnomad_variants: List of variant IDs present in gnomAD
        
        Returns:
            List of dicts with structure: {'variant': str, 'pos': int, 'ref': str,
            'alt': str, 'ucsc_url': str, 'mitomap_url': str, 'gnomad_url': str,
            'dbsnp_id': str, 'dbsnp_url': str, 'clinvar_url': str, 'clinvar_vcv': str}
        """
        # TODO: Parse VCF fields and CSQ annotation
        # TODO: Adjust variant IDs for indels (gnomAD format)
        # TODO: Extract dbSNP and ClinVar IDs
        # TODO: Generate MT-specific URLs (MitoMap, etc.)
        # TODO: Return unique MT annotations
        raise NotImplementedError("MtsTransformer.transform() not yet implemented")


class GenomicIbvlFrequenciesTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
        raise NotImplementedError("GenomicIbvlFrequenciesTransformer.transform() not yet implemented")


class GenomicGnomadFrequenciesTransformer(BaseTransformer):
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
    
    def transform(self, gnomad_records: List[Dict[str, Any]], 
                  ibvl_variants: List[str],
                  assembly: str) -> List[Dict[str, Any]]:
        """
        Extract gnomAD population frequencies for cohort variants.
        
        Args:
            gnomad_records: Pre-processed gnomAD records from TSV
            ibvl_variants: List of variant IDs present in the cohort
            assembly: Genome assembly ("GRCh37" or "GRCh38")
        
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
        raise NotImplementedError("GenomicGnomadFrequenciesTransformer.transform() not yet implemented")


class MtIbvlFrequenciesTransformer(BaseTransformer):
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
    
    def transform(self, vcf_records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
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
        raise NotImplementedError("MtIbvlFrequenciesTransformer.transform() not yet implemented")


class MtGnomadFrequenciesTransformer(BaseTransformer):
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
    
    def transform(self, gnomad_mt_records: List[Dict[str, Any]], 
                  ibvl_variants: List[str]) -> List[Dict[str, Any]]:
        """
        Extract gnomAD MT population frequencies for cohort variants.
        
        Args:
            gnomad_mt_records: gnomAD MT records from TSV
            ibvl_variants: List of MT variant IDs in the cohort
        
        Returns:
            List of dicts with structure: {'variant': str, 'an': int, 'ac_hom': int,
            'ac_het': int, 'af_hom': float, 'af_het': float, 'max_hl': float}
        """
        # TODO: Create variant IDs from gnomAD MT records
        # TODO: Adjust indel variant IDs
        # TODO: Intersect with IBVL MT variants
        # TODO: Return gnomAD MT frequencies for matching variants
        raise NotImplementedError("MtGnomadFrequenciesTransformer.transform() not yet implemented")

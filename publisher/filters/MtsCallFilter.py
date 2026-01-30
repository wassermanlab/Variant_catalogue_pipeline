from .CallFilter import CallFilter
from typing import List, Dict, Any, Optional
from constants import NA

class MtsCallFilter(CallFilter):
    """
    Generates the 'mts' table (mitochondrial variant annotations).
    """
    def __init__(self, vcf_file_path: str, assembly: Optional[str] = None):
        super().__init__(vcf_file_path)
        self.assembly = assembly
        self.gnomad_variants = set()
    def getTableRows(self) -> List[Dict[str, Any]]:
        mts = {}
        for record in self.vcf_records:
            variant = self.make_variant_id(record)
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

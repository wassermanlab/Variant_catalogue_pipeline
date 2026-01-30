from .CallFilter import CallFilter
from typing import List, Dict, Any, Optional
from constants import NA, CHR_NOTATION

class SnvsCallFilter(CallFilter):
    """
    Generates the 'snvs' table (SNV-specific annotations).
    """
    def __init__(self, vcf_file_path: str, assembly: Optional[str] = None):
        super().__init__(vcf_file_path)
        self.assembly = assembly
    def getTableRows(self) -> List[Dict[str, Any]]:
        snvs = {}
        for record in self.vcf_records:
            variant = self.make_variant_id(record)
            chrom = record.CHROM.replace("chr", "") if not CHR_NOTATION else record.CHROM
            pos = record.POS
            ref = record.REF
            alt = record.ALT[0].value  # assuming single ALT allele

            class_list = self.get_csq_values(record, "VARIANT_CLASS")
            cadd_phred_list = record.INFO.get("CADD_PHREDscore", self.get_csq_values(record, "CADD_PHRED")) # variome: CADD_PHREDscore in info; ibvl: CADD_PHRED in csq!
            existing_variation_list = self.get_csq_values(record, "Existing_variation")
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
            variant_class = record.INFO.get("TYPE")[0] if record.INFO.get("TYPE") else (class_list[0] if class_list else NA)
            cadd_phred = float(cadd_phred_list[0]) if cadd_phred_list and cadd_phred_list[0] != "" else None
            
            # Determine variant length
            if variant_class in ["SNV", "SNP"]:
                var_length = 1
            elif variant_class in ["INS", "DEL"]:
                var_length = abs(len(alt) - len(ref))
            elif variant_class in ["INDEL"]:
                var_length = len(alt)
            else:
                var_length = None

            # CADD interpretation
            if cadd_phred is not None:
                cadd_intr = "Damaging" if cadd_phred > 15 else "Tolerable"
            else:
                cadd_intr = NA

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
                continue
        return list(snvs.values())

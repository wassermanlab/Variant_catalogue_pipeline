from .CallFilter import CallFilter
from typing import List, Dict, Any
from constants import NA

class MtBvlFrequenciesCallFilter(CallFilter):
    """
    Generates the 'mt_bvl_frequencies' table.
    """
    def getTableRows(self) -> List[Dict[str, Any]]:
        rows = []
        for record in self.vcf_records:
            variant = self.make_variant_id(record)
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
            rows.append({
                'variant': variant,
                'an': an,
                'ac_hom': ac_hom,
                'ac_het': ac_het,
                'af_hom': af_hom,
                'af_het': af_het,
                'hl_hist': hl_hist,
                'max_hl': max_hl
            })
        return rows

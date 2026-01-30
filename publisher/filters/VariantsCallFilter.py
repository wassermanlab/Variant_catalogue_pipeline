from .CallFilter import CallFilter
from typing import List, Dict, Any

from constants import NA

class VariantsCallFilter(CallFilter):
    """
    Generates the 'variants' table (master variant list).
    """
    def __init__(self, vcf_file_path: str):
        super().__init__(vcf_file_path)
    def getTableRows(self) -> List[Dict[str, Any]]:
        variants = {}
        for record in self.vcf_records:
            variant_id = self.make_variant_id(record)
            filter = ";".join(record.FILTER)
            for csq in record.INFO.get("CSQ", []):
                var_type = self.get_csq_values(record, "VARIANT_CLASS")
                if var_type == []:
                    var_type = self.get_info_value(record, "TYPE")
                filter =  ";".join(record.FILTER)
                if variant_id not in variants:
                    variants[variant_id] = {
                        'variant_id': variant_id, 
                        'var_type': var_type[0] if var_type and var_type != [] else NA ,
                        'filter': filter if filter != "" else NA
                    }
                else:
                    continue
        return list(variants.values())

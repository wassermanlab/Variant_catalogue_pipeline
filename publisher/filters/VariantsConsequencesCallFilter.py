from .CallFilter import CallFilter
from typing import List, Dict, Any
import logging
logger = logging.getLogger(__name__)

class VariantsConsequencesCallFilter(CallFilter):
    """
    Generates the 'variants_consequences' table.
    """
    def getTableRows(self) -> List[Dict[str, Any]]:
        variant_consequences = []
        for record in self.vcf_records:
            variant = self.make_variant_id(record)
            transcript_list = self.get_csq_values(record, "Feature")
            consequence_list = self.get_csq_values(record, "Consequence")
            l = len(transcript_list)
            if l != len(consequence_list):
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

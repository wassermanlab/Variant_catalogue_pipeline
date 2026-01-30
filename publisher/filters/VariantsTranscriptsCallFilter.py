from .CallFilter import CallFilter
from typing import List, Dict, Any
import logging
logger = logging.getLogger(__name__)

class VariantsTranscriptsCallFilter(CallFilter):
    """
    Generates the 'variants_transcripts' table.
    """
    def getTableRows(self) -> List[Dict[str, Any]]:
        variantsTranscripts = []
        for record in self.vcf_records:
            transcript = self.get_csq_values(record, "Feature")
            variant = self.make_variant_id(record)
            hgvsc = self.get_csq_values(record, "HGVSc")
            l = len(transcript)
            lh = len(hgvsc)
            if l != lh:
                import logging
                logging.warning(f"mismatched lengths for transcript and hgvsc: {l} vs {lh}")
                continue
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

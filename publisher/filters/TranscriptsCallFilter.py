from .CallFilter import CallFilter
from typing import List, Dict, Any
from constants import NA
import logging
logger = logging.getLogger(__name__)

class TranscriptsCallFilter(CallFilter):
    """
    Generates the 'transcripts' table.
    """
    def getTableRows(self) -> List[Dict[str, Any]]:
        """
        Associate transcripts with genes from VEP annotations.
        Returns:
            List of dicts with structure: 
            {'transcript_id': str, 'gene': str, 'transcript_type': str, 'tsl': str}
        """
        transcripts = {}
        for record in self.vcf_records:
            feature = self.get_csq_values(record, "Feature")
            if feature == "" or feature == "NA":
                logger.warning("skipping transcript with no feature: %s", feature)
                continue
            symbol = self.get_csq_values(record, "SYMBOL")
            source = self.get_csq_values(record, "SOURCE")
            tsl = self.get_csq_values(record, "TSL")
            biotype = self.get_csq_values(record, "BIOTYPE")
            l = len(feature)
            if not (l == len(symbol) == len(source) == len(tsl)):
                logger.warning("mismatched lengths for transcript feature, symbol, source, tsl: %d vs %d vs %d vs %d", l, len(symbol), len(source), len(tsl))
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
                        'tsl': tsl[i],
                        'biotype': biotype[i] if biotype and len(biotype) > i else NA
                    }
        return list(transcripts.values())

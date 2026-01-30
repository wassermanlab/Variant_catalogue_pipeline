from .CallFilter import CallFilter
from typing import List, Dict, Any

class GenesCallFilter(CallFilter):
    """
    Generates the 'genes' table.
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

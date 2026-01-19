"""
VCF Publisher
This module orchestrates the transformation of VCF files and related data sources
into database tables
"""

from typing import Dict, List, Any, Optional
from pathlib import Path
import logging
import os

from VCF_filters import (
    GenesCallFilter,
    TranscriptsCallFilter,
    VariantsCallFilter,
    VariantsTranscriptsCallFilter,
    VariantsAnnotationsCallFilter,
    VariantsConsequencesCallFilter,
    SnvsCallFilter,
    MtsCallFilter,
    GenomicIbvlFrequenciesCallFilter,
    GenomicGnomadFrequenciesCallFilter,
    MtIbvlFrequenciesCallFilter,
    MtGnomadFrequenciesCallFilter,
)

vcfs = [
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_case/SNV/SNV_filtered_frequ_only_SNV_annotation_table_merged_22_truncated.vcf')
]

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

class VariantPublisher:
    
    def __init__(self):
        pass
    
    def start(self, config: Dict[str, Any]) -> None:
        
        results = {}
        
        # Transform data using appropriate transformers
        results['genes'] = GenesCallFilter(vcfs).getTableRows()
        results['transcripts'] = TranscriptsCallFilter(vcfs).getTableRows()
        results['variants'] = VariantsCallFilter(vcfs).getTableRows()
        results['variants_transcripts'] = VariantsTranscriptsCallFilter(vcfs).getTableRows()
#        results['variants_annotations'] = VariantsAnnotationsCallFilter(vcfs).getTableRows()
#        results['variants_consequences'] = VariantsConsequencesCallFilter(vcfs).getTableRows()
 #       results['snvs'] = SnvsCallFilter(vcfs).getTableRows( )
#        results['genomic_ibvl_frequencies'] = GenomicIbvlFrequenciesCallFilter(vcfs).getTableRows()
        
        for table_name, records in results.items():
            for r in records[:5]:
                logger.info(f"SNV Table {table_name} record: {r}")
            #send it to the database

def main():
    """Command-line entry point."""
    import argparse
    import pprint
    
    parser = argparse.ArgumentParser(
        description='filter (extract) BVL data out from VCF files and publish to database'
    )
    parser.add_argument('--assembly', choices=['GRCh37', 'GRCh38'],
                       help='Genome assembly', default='GRCh38')
    parser.add_argument('--severity-table', default='severities.tsv', 
                       help='Path to severity_table.tsv')
    parser.add_argument('--gnomad-snv-tsv', help='Path to gnomAD SNV TSV')
    parser.add_argument('--gnomad-mt-tsv', help='Path to gnomAD MT TSV')
    
    args = parser.parse_args()
    
    config = {
        'assembly': args.assembly,
        'severity_table_path': args.severity_table,
    }
    logger.info("configuration:\n%s", pprint.pformat(config))
    
    publish_job = VariantPublisher()
    publish_job.start(config)
    logger.info("completed")

if __name__ == '__main__':
    main()
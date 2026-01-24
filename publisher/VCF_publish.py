"""
VCF Publisher
This module orchestrates the transformation of VCF files and related data sources
into database tables
"""

from typing import Dict, List, Any, Optional
from pathlib import Path
import logging
import os
from datetime import datetime

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
#    GenomicGnomadFrequenciesCallFilter,
    MtIbvlFrequenciesCallFilter,
#    MtGnomadFrequenciesCallFilter,
)

snv_vcfs = [
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_case/HG002-4_chr21_SNV_v7.vcf')
#    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures/vcf/mock_snv.vcf')
]

mt_vcfs = [
#    os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_case/HG002-4_MT_v3.vcf')
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
        now = datetime.now()
        
        results['genes'] = GenesCallFilter(snv_vcfs).getTableRows()
        results['transcripts'] = TranscriptsCallFilter(snv_vcfs).getTableRows()
        results['variants'] = VariantsCallFilter(snv_vcfs).getTableRows()
        results['variants_transcripts'] = VariantsTranscriptsCallFilter(snv_vcfs).getTableRows()
        results['variants_annotations'] = VariantsAnnotationsCallFilter(snv_vcfs).getTableRows()
        results['variants_consequences'] = VariantsConsequencesCallFilter(snv_vcfs).getTableRows()
        results['snvs'] = SnvsCallFilter(snv_vcfs).getTableRows( )
        results['genomic_ibvl_frequencies'] = GenomicIbvlFrequenciesCallFilter(snv_vcfs).getTableRows()
        results['mts'] = MtsCallFilter(mt_vcfs).getTableRows()
        results['mt_ibvl_frequencies'] = MtIbvlFrequenciesCallFilter(mt_vcfs).getTableRows()
        
        for table_name, records in results.items():
            for r in records[:5]:
                logger.info(f"SNV Table {table_name} record: {r}")
            #send it to the database
            
        output_dir = Path(os.path.dirname(os.path.abspath(__file__))) / "output"
        output_dir.mkdir(parents=True, exist_ok=True)

        for table_name, records in results.items():
            if not records:
                continue
            output_path = output_dir / f"{table_name}.tsv"
            with open(output_path, "w", encoding="utf-8") as f:
                # Write header
                header = records[0].keys() if isinstance(records[0], dict) else []
                if header:
                    f.write("\t".join(header) + "\n")
                # Write rows
                for row in records:
                    if isinstance(row, dict):
                        f.write("\t".join(str(row.get(col, "")) for col in header) + "\n")
                    
        logger.info(f"TSV files written to {output_dir}")
        logger.info(f"Processing completed at {datetime.now()}")
        duration = datetime.now() - now
        logger.info(f"Total duration: {duration}")
        duration_seconds = duration.total_seconds()
        logger.info(f"Total duration in seconds: {duration_seconds}")

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
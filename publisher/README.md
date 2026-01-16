# Variant Catalogue Publisher

This module provides Python classes to transform VCF files and related data sources directly into database-ready structures, replacing the R script + TSV workflow.

## Overview

The publisher module consists of:

- **transformers.py**: Stub classes for each database table transformation
- **pipeline.py**: Main orchestrator that coordinates reading data, transforming, and loading to database
- **test_transformers.py**: Unit tests with mock data defining expected input/output structures

## Architecture

Each transformer class:
- Corresponds to one database table
- Implements a `transform()` method
- Takes raw data sources as input (VCF records, external files)
- Returns a list of dictionaries (analogous to TSV rows)

## Usage

### Running Tests

```bash
python -m pytest publisher/test_transformers.py -v
```

or

```bash
python publisher/test_transformers.py
```

### Using the Pipeline

```python
from publisher.pipeline import VariantCataloguePipeline

pipeline = VariantCataloguePipeline(db_connection=your_db_conn)

# Process SNV data
snv_results = pipeline.process_snv_vcf(
    vcf_path='path/to/snv.vcf',
    assembly='GRCh38',
    severity_table_path='severity_table.tsv',
    gnomad_tsv_path='gnomad.tsv'
)

# Results is a dict mapping table names to record lists
for table_name, records in snv_results.items():
    pipeline.load_to_database(table_name, records)
```

### Command Line

```bash
python -m publisher.pipeline \
    --snv-vcf path/to/snv.vcf \
    --mt-vcf path/to/mt.vcf \
    --assembly GRCh38 \
    --severity-table severity_table.tsv \
    --gnomad-snv-tsv gnomad_snv.tsv \
    --gnomad-mt-tsv gnomad_mt.tsv
```

## Transformer Classes

### Core Tables

1. **GenesTransformer** - Extracts unique gene symbols
2. **TranscriptsTransformer** - Associates transcripts with genes
3. **VariantsTransformer** - Master variant list with type classification
4. **VariantsTranscriptsTransformer** - Links variants to transcripts with HGVS
5. **VariantsAnnotationsTransformer** - Protein-level annotations (SIFT, PolyPhen)
6. **VariantsConsequencesTransformer** - Consequence severity scores

### Variant-Specific Tables

7. **SnvsTransformer** - SNV annotations with CADD scores and browser URLs
8. **MtsTransformer** - Mitochondrial variant annotations

### Frequency Tables

9. **GenomicIbvlFrequenciesTransformer** - Internal cohort frequencies (stratified by sex)
10. **GenomicGnomadFrequenciesTransformer** - gnomAD population frequencies
11. **MtIbvlFrequenciesTransformer** - MT frequencies with heteroplasmy data
12. **MtGnomadFrequenciesTransformer** - gnomAD MT population frequencies

## Data Sources

### VCF Files
- Annotated with VEP (Variant Effect Predictor)
- Include Hail-calculated frequency fields in INFO
- MT VCFs include heteroplasmy data in GT fields

### External Files
- **severity_table.tsv**: Maps VEP consequence terms to numeric severity
- **gnomAD TSV files**: Pre-processed population frequency data
- **gnomAD MT TSV**: Mitochondrial population frequencies

## Implementation Status

All transformers are currently **stubs** (raise `NotImplementedError`). The framework provides:

- Type hints and docstrings defining inputs/outputs
- Comments describing transformation logic
- Unit tests with mock data nailing down data structures

## Next Steps

1. Implement VCF parsing (consider pysam or cyvcf2)
2. Implement CSQ field parsing and VEP annotation extraction
3. Implement each transformer's `transform()` method
4. Implement database loading logic
5. Add integration tests with real VCF files
6. Add error handling and validation

## Testing

Tests validate:
- Each transformer raises `NotImplementedError` (current state)
- Expected output structures are documented
- Mock data represents realistic VCF structures

Run tests to verify the framework:
```bash
python publisher/test_transformers.py
```

All tests should pass (verifying stubs are in place).


## Common Data Sources

### Pipeline VCF Files
All variant types are first processed through:
1. **Variant calling tools** (DeepVariant for SNVs, Mutect2 for MT, Manta/Smoove for SVs, MELT for MEIs, ExpansionHunter for STRs)
2. **Hail** for quality control and frequency calculations (AF, AC, AN, homozygote counts stratified by sex: tot/XX/XY)
3. **VEP (Variant Effect Predictor)** for functional annotation

### External Reference Files
- **gnomAD VCF files**: Population frequency databases for SNVs and MT variants
- **severity_table.tsv**: Maps VEP consequence terms to numeric severity scores
- **STR variant catalogue JSON**: Reference for STR loci with gene annotations

---

## Table Descriptions

### 1. genes
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L320), [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L232), [`SV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SV_data_organization.R#L289), [`MEI_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MEI_data_organization.R#L224), [`STR_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/STR_data_organization.R#L150)

**Data Source**: VCF INFO field, VEP annotation CSQ subfield `SYMBOL`

**Description**: Extracts unique gene symbols from VEP annotations.

**Processing**:
- Parse the CSQ annotation from VCF INFO field
- Split by pipe delimiter to extract SYMBOL field
- Filter out NA/empty values
- Output unique gene short names

**Output Fields**: `short_name`

---

### 2. transcripts
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L335), [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L248)

**Data Source**: VCF INFO field, VEP annotation CSQ subfields: `Feature`, `SYMBOL`, `SOURCE`, `TSL`

**Description**: Associates transcripts with genes and indicates transcript source (Ensembl or RefSeq).

**Processing**:
- Parse CSQ annotation from VCF INFO field
- Extract: transcript ID (Feature), gene (SYMBOL), source (SOURCE), transcript support level (TSL)
- Recode SOURCE: "Ensembl" → "E", "RefSeq" → "R"
- Filter out entries without transcript IDs
- Output unique transcript records

**Output Fields**: `transcript_id`, `gene`, `transcript_type`, `tsl`

**Note**: SV, MEI, and STR variants are not associated with specific transcripts, so these tables are only generated for SNV and MT variants.

---

### 3. variants
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L291), [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L163), [`SV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SV_data_organization.R#L276), [`MEI_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MEI_data_organization.R#L211), [`STR_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/STR_data_organization.R#L142)

**Data Source**: VCF ID field (variant identifier: chr_pos_ref_alt format)

**Description**: Master variant list with variant type classification.

**Processing**:
- Extract variant ID from VCF ID field
- Assign variant type: "SNV" for single nucleotide variants, "MT" for mitochondrial variants, "SV" for structural variants
- Output unique variant identifiers with type

**Output Fields**: `variant_id`, `var_type`

---

### 4. variants_transcripts
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L246), [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L174)

**Data Source**: VCF INFO field, VEP annotation CSQ subfields: `Feature`, variant ID, `HGVSc`

**Description**: Links variants to transcripts with HGVS coding sequence notation.

**Processing**:
- Parse CSQ annotation from VCF INFO field
- Extract: transcript ID (Feature), variant ID, HGVSc (coding sequence change)
- Filter out intergenic variants (where transcript is NA)
- Output unique variant-transcript associations

**Output Fields**: `transcript`, `variant`, `hgvsc`

---

### 5. variants_annotations
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L283), [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L205)

**Data Source**: VCF INFO field, VEP annotation CSQ subfields: `HGVSp`, `SIFT`, `PolyPhen`, `Feature`, variant ID

**Description**: Protein-level annotations for variants including pathogenicity predictions.

**Processing**:
- Parse CSQ annotation from VCF INFO field
- Extract: HGVSp (protein change), SIFT score/prediction, PolyPhen score/prediction, transcript ID, variant ID
- Decode URL-encoded characters in HGVSp (%3D → =)
- Filter to entries with valid HGVSp values (protein-coding variants only)
- Output unique variant-transcript annotation records

**Output Fields**: `hgvsp`, `sift`, `polyphen`, `transcript`, `variant`

---

### 6. variants_consequences
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L265), [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L194)

**Data Source**: VCF INFO field, VEP annotation CSQ subfield `Consequence`, severity_table.tsv

**Description**: Maps variants to transcripts with numeric severity scores based on consequence type.

**Processing**:
- Parse CSQ annotation from VCF INFO field
- Extract Consequence field (may contain multiple consequences separated by "&")
- Split compound consequences into separate rows
- Join with severity_table.tsv to convert consequence terms to numeric severity
- Filter out intergenic variants
- Output variant-transcript-severity associations

**Output Fields**: `severity`, `variant`, `transcript`

---

### 7. snvs
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L234)

**Data Source**: VCF (CHROM, POS, ID, REF, ALT, QUAL, INFO field with VEP CSQ)

**Description**: SNV-specific annotation table with variant details, CADD scores, database IDs, and browser URLs.

**Processing**:
- Parse VCF fixed fields and VEP CSQ annotation
- Extract VARIANT_CLASS, CADD_PHRED, Existing_variation, VAR_SYNONYMS
- Calculate variant length (1 for SNVs, or ref/alt length difference for indels)
- Derive CADD interpretation: ≤15 = "Tolerable", >15 = "Damaging"
- Calculate maximum SpliceAI score from DS_AG, DS_AL, DS_DG, DS_DL fields
- Extract dbSNP ID from Existing_variation (if contains "rs")
- Extract ClinVar VCV number from VAR_SYNONYMS
- Generate URLs:
  - dbSNP: `https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rsID>`
  - UCSC: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chr<chr>:<pos>-<pos>&position=chr<chr>:<pos-25>-<pos+25>`
  - Ensembl: `https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>:<pos-25>-<pos+25>` (GRCh38) or grch37.ensembl.org (GRCh37)
  - ClinVar: `https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/`
  - gnomAD: `https://gnomad.broadinstitute.org/variant/<chr>-<pos>-<ref>-<alt>?dataset=gnomad_r3` (GRCh38) or gnomad_r2_1 (GRCh37)
- Filter to unique variants

**Output Fields**: `variant`, `type`, `length`, `chr`, `pos`, `ref`, `alt`, `cadd_score`, `cadd_intr`, `dbsnp_id`, `dbsnp_url`, `ucsc_url`, `ensembl_url`, `clinvar_url`, `gnomad_url`, `clinvar_vcv`, `splice_ai`

---

### 8. mts
**Script**: [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L154)

**Data Source**: VCF (CHROM, POS, ID, REF, ALT, INFO field with VEP CSQ)

**Description**: Mitochondrial variant annotation table with MT-specific details and URLs.

**Processing**:
- Parse VCF fixed fields and VEP CSQ annotation
- Adjust variant IDs for indels to match gnomAD format (position +1, trimmed sequence)
- Extract dbSNP ID and ClinVar VCV from VEP annotations
- Generate URLs:
  - UCSC: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM:<pos>-<pos>&position=chrM:<pos-25>-<pos+25>`
  - MitoMap: `https://mitomap.org/cgi-bin/search_allele?variant=<pos><ref>><alt>`
  - gnomAD: `https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3` (if in gnomAD)
  - dbSNP: `https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rsID>`
  - ClinVar: `https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/`
- Filter to unique variants

**Output Fields**: `variant`, `pos`, `ref`, `alt`, `ucsc_url`, `mitomap_url`, `gnomad_url`, `dbsnp_id`, `dbsnp_url`, `clinvar_url`, `clinvar_vcv`

**Note**: The "svs" naming convention is also used for structural variants (see below).

---

### 9. genomic_ibvl_frequencies
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L225), [`SV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SV_data_organization.R#L237), [`MEI_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MEI_data_organization.R#L179)

**Data Source**: VCF INFO field, Hail-calculated frequencies: `AF_tot_XX_XY`, `AC_tot_XX_XY`, `AN_tot_XX_XY`, `hom_tot_XX_XY`

**Description**: Internal cohort (IBVL) allele frequencies stratified by sex.

**Processing**:
- Parse Hail-added INFO fields
- Each field contains comma-separated values: total, XX (female), XY (male)
- Split values: 
  - af_tot, af_xx, af_xy (allele frequencies)
  - ac_tot, ac_xx, ac_xy (allele counts)
  - an_tot, an_xx, an_xy (allele numbers - total chromosomes tested)
  - hom_tot, hom_xx, hom_xy (homozygote counts)
- Include QUAL field from VCF
- Output unique variant frequency records

**Output Fields**: `variant`, `af_tot`, `af_xx`, `af_xy`, `ac_tot`, `ac_xx`, `ac_xy`, `an_tot`, `an_xx`, `an_xy`, `hom_tot`, `hom_xx`, `hom_xy`, `quality`

---

### 10. genomic_gnomad_frequencies
**Script**: [`SNV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SNV_data_organization.R#L311) (uses [`gnomad_frequency_table.nf`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/gnomad_frequency_table.nf))

**Data Source**: External gnomAD VCF files, processed through `gnomad_frequency_table.nf`

**Description**: gnomAD population frequencies for variants present in the cohort.

**Processing**:
1. **Pre-processing** (gnomad_frequency_table.nf): Extract subset of fields from gnomAD VCF using GATK VariantsToTable:
   - Fields: CHROM, POS, REF, ALT, QUAL, FILTER, AF, AC, AN, nhomalt
   - Additional GRCh38 fields: exomes_filters, genomes_filters
   - Split by chromosome for efficiency
2. **R script processing**:
   - Create variant ID: chr_pos_ref_alt
   - Adjust chromosome labels to match pipeline format (remove "chr" prefix if present)
   - Intersect with IBVL variants (keep only variants present in cohort)
   - Output gnomAD frequencies for matching variants

**Output Fields**: 
- GRCh37: `variant`, `af_tot`, `ac_tot`, `an_tot`, `hom_tot`, `FILTER`
- GRCh38: `variant`, `af_tot`, `ac_tot`, `an_tot`, `hom_tot`, `FILTER`, `exomes_filters`, `genomes_filters`

---

### 11. mt_ibvl_frequencies
**Script**: [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L145)

**Data Source**: VCF GT (genotype) fields with Hail-calculated MT-specific metrics

**Description**: Mitochondrial variant frequencies including heteroplasmy information.

**Processing**:
- Parse Hail-added GT fields (not INFO):
  - AC_hom, AC_het (homoplasmic and heteroplasmic allele counts)
  - AF_hom, AF_het (homoplasmic and heteroplasmic allele frequencies)
  - AN (total allele number)
  - max_observed_heteroplasmy (maximum heteroplasmy level observed)
  - heteroplasmy_histogram (distribution of heteroplasmy levels)
- Extract histogram values and format as comma-separated string
- Filter out variants where AN=0 (no observations)
- Adjust variant IDs for indels to match gnomAD format
- Output unique MT variant frequencies

**Output Fields**: `variant`, `an`, `ac_hom`, `ac_het`, `af_hom`, `af_het`, `hl_hist`, `max_hl`

---

### 12. mt_gnomad_frequencies
**Script**: [`MT_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MT_data_organization.R#L219)

**Data Source**: External gnomAD mitochondrial TSV file (pre-processed, not VCF format)

**Description**: gnomAD mitochondrial population frequencies for variants in the cohort.

**Processing**:
- Read pre-processed gnomAD MT table (TSV format with fields: chromosome, position, ref, alt, AN, AC_hom, AC_het, AF_hom, AF_het, max_observed_heteroplasmy)
- Create variant ID: chr_pos_ref_alt
- Adjust IDs for indels to match format
- Intersect with IBVL MT variants
- Output gnomAD frequencies for matching variants

**Output Fields**: `variant`, `an`, `ac_hom`, `ac_het`, `af_hom`, `af_het`, `max_hl`

---

## Additional SV-Specific Tables

The pipeline also generates SV-specific consequence and annotation tables that differ from the SNV/MT pattern:

### svs (Structural Variants)
**Scripts**: [`SV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SV_data_organization.R#L246), [`MEI_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MEI_data_organization.R#L188)

**Data Source**: VCF with SV-specific INFO fields (SVLEN, SVTYPE, AVG_START, AVG_END, etc.)

**Description**: Structural variant annotations including SV type, size, and genomic coordinates.

**Processing**:
- Parse SV-specific INFO fields:
  - For SVs: SVLEN, SVTYPE, AVG_START, AVG_END, AVG_LEN, IDLIST (to determine algorithm: Manta or Smoove)
  - For MEIs: SVLEN, SVTYPE, TSD, ASSESS, INTERNAL (algorithm is always MELT)
- Extract VEP CSQ annotations
- Calculate coordinates and construct UCSC browser URLs
- Generate gnomAD SV region URLs (GRCh37 only): `https://gnomad.broadinstitute.org/region/<chr>-<start>-<end>?dataset=gnomad_sv_r2_1`

**Output Fields**: `variant`, `chr1`, `chr1_pos1`, `chr1_pos2`, `sv_type`, `sv_length`, `algorithm`, `ucsc_url`, `gnomad_id`, `gnomad_url`

### sv_consequences
**Scripts**: [`SV_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/SV_data_organization.R), [`MEI_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/MEI_data_organization.R), [`STR_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/STR_data_organization.R#L146)

**Data Source**: VCF INFO field, VEP CSQ Consequence field

**Description**: Gene-level consequences for structural variants (no transcript associations).

**Processing**:
- Parse CSQ Consequence field
- Split compound consequences (separated by "&") into separate rows
- Associate with gene symbols
- No transcript-level associations (SVs affect genes broadly)
- For STR: consequence is "NA", gene comes from JSON variant catalogue

**Output Fields**: `gene`, `variant`, `consequence`

### str
**Script**: [`STR_data_organization.R`](https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/main/modules/STR_data_organization.R#L138)

**Data Source**: VCF (REF, RU, REPID, END from INFO), ExpansionHunter variant catalogue JSON, genotype REPCN field

**Description**: Short tandem repeat allele distributions and repeat unit information.

**Processing**:
- Read VCF and extract STR-specific fields:
  - RU: repeat unit sequence
  - REF: reference repeat count
  - END: end position
  - REPID: repeat ID
- Read JSON variant catalogue to map genomic regions to gene symbols (LocusId)
- Parse REPCN genotype field (repeat count) for all samples
- Calculate allele size distribution:
  - Extract allele 1 and allele 2 from each sample
  - Missing genotypes = homozygous reference
  - Count frequency of each allele size
  - Format as comma-separated: "size1:count1,size2:count2,..."
- Calculate min and max repeat counts observed
- Generate gnomAD STR URL: `https://gnomad.broadinstitute.org/short-tandem-repeat/<gene>?dataset=gnomad_r3`

**Output Fields**: `variant`, `repeat_unit`, `min_n_repeat`, `max_n_repeat`, `allele_distrib`

---
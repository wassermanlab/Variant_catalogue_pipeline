# helper for local environment
# TODO add position, remove uniqueness on variant

from sqlalchemy import (
    create_engine,
    MetaData,
    Table,
    Column,
    Integer,
    String,
    Float,
    UniqueConstraint,
    ForeignKey,
    text
)

from sqlalchemy.exc import DataError, IntegrityError, ProgrammingError
import pandas as pd
import sys
import os
from dotenv import load_dotenv
from datetime import datetime

load_dotenv()

# get command line arguments
rootDir = os.environ.get("PIPELINE_OUTPUT_PATH")
chunk_size = int(os.environ.get("CHUNK_SIZE"))
verbose = os.environ.get("VERBOSE") == "true"
dbConnectionString = os.environ.get("DB")
# Replace 'sqlite:///:memory:' with your actual database connection string

if (False and verbose):
    engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})
else:
    engine = create_engine(dbConnectionString, pool_pre_ping=True, pool_recycle=3600, connect_args={'autocommit': True})

metadata = MetaData()

# Define the "genes" table
genes_table = Table(
    "genes",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("short_name", String(30), nullable=False),
    UniqueConstraint("short_name", name="unique"),
)


transcripts_table = Table(
    "transcripts",
    metadata,
    Column("transcript_id", String(255)),
    Column("id", Integer, primary_key=True),
    Column("gene", Integer, ForeignKey("genes.id", ondelete="CASCADE")),
    Column("transcript_type", String(1)),
    Column("tsl", String(255)),
    UniqueConstraint("transcript_id", name="transcripts_unique"),
#    ForeignKey("gene", name="transcripts_genes_fk", ondelete="CASCADE", table="genes"),
)

# Define the "genomic_gnomad_frequencies" table
genomic_gnomad_frequencies_table = Table(
    "genomic_gnomad_frequencies",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("af_tot", Float),
    Column("ac_tot", Integer),
    Column("an_tot", Integer),
    Column("hom_tot", Integer),
    UniqueConstraint("variant", name="geno_gnomad_freq_unique"),
)

# Define the "genomic_ibvl_frequencies" table
genomic_ibvl_frequencies_table = Table(
    "genomic_ibvl_frequencies",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("af_tot", Float),
    Column("af_xx", Float),
    Column("af_xy", Float),
    Column("ac_tot", Integer),
    Column("ac_xx", Integer),
    Column("ac_xy", Integer),
    Column("an_tot", Integer),
    Column("an_xx", Integer),
    Column("an_xy", Integer),
    Column("hom_tot", Integer),
    Column("hom_xx", Integer),
    Column("hom_xy", Integer),
    Column("quality", Integer),
    UniqueConstraint("variant", name="geno_ibvl_freq_unique")
)

mt_gnomad_frequencies_table = Table(
    "mt_gnomad_frequencies",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("an", Integer),
    Column("ac_hom", Integer),
    Column("ac_het", Integer),
    Column("af_hom", Float),
    Column("af_het", Float),
    Column("max_hl", Integer),
    UniqueConstraint("variant", name="mt_gnomad_freq_unique")
)

mt_ibvl_frequencies_table = Table(
    "mt_ibvl_frequencies",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("an", Integer),
    Column("ac_hom", Integer),
    Column("ac_het", Integer),
    Column("af_hom", Float),
    Column("af_het", Float),
    Column("hl_hist", String(30)),
    Column("max_hl", Integer),
    UniqueConstraint("variant", name="mt_ibvl_freq_unique"),
)


mts_table = Table(
    "mts",
    metadata,
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("id", Integer, primary_key=True),
    Column("pos", Integer),
    Column("ref", String(60)),
    Column("alt", String(30)),
    Column("ucsc_url", String(511)),
    Column("mitomap_url", String(511)),
    Column("gnomad_url", String(511)),
    Column("dbsnp_id", String(30)),
    Column("dbsnp_url", String(511)),
    Column("clinvar_url", String(511)),
    Column("clinvar_vcv", Integer),
    UniqueConstraint("variant", name="mts_unique"),
)

snvs_table = Table(
    "snvs",
    metadata,
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("id", Integer, primary_key=True),
    Column("type", String(30)),
    Column("length", Integer),
    Column("chr", String(2)),
    Column("pos", Integer),
    Column("ref", String(255)),
    Column("alt", String(255)),
    Column("cadd_intr", String(255)),
    Column("cadd_score", Integer),
    Column("dbsnp_id", String(30)),
    Column("dbsnp_url", String(511)),
    Column("ucsc_url", String(511)),
    Column("ensembl_url", String(511)),
    Column("clinvar_vcv", Integer),
    Column("clinvar_url", String(511)),
    Column("gnomad_url", String(511)),
    Column("splice_ai", Integer),
    UniqueConstraint("variant", name="snvs_unique")
)

# Define the "str" table
str_table = Table(
    "str",
    metadata,
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("id", Integer, primary_key=True),
    Column("repeat_unit", String(20)),
    Column("min_n_repeat", Integer),
    Column("max_n_repeat", Integer),
    Column("allele_distrib", String(255)),
    Column("reference_region", String(40)),
    UniqueConstraint("variant", name="str_unique")
)

# Define the "sv_consequences" table
sv_consequences_table = Table(
    "sv_consequences",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("gene", Integer, ForeignKey("genes.id", ondelete="CASCADE")),
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("consequence", String(255)),
    UniqueConstraint("variant", name="sv_consequences_unique")
)

svs_table = Table(
    "svs",
    metadata,
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("id", Integer, primary_key=True),
    Column("chr1", String(2)),
    Column("chr1_pos1", Integer),
    Column("chr1_pos2", Integer),
    Column("sv_type", String(30)),
    Column("sv_length", Integer),
    Column("algorithm", String(30)),
    Column("ucsc_url", String(511)),
    Column("gnomad_id", String(30)),
    Column("gnomad_url", String(511)),
    UniqueConstraint("variant", name="svs_unique"),
#    ForeignKey("variant", name="svs_fk", ondelete="CASCADE", table="variants"),
)

svs_ctx_table = Table(
    "svs_ctx",
    metadata,
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("id", Integer, primary_key=True),
    Column("chr2", String(2)),
    Column("chr2_pos1", Integer),
    Column("ucsc_url2", String(511)),
    Column("gnomad_id2", String(30)),
    Column("gnomad_url2", String(511)),
    UniqueConstraint("variant", name="svs_ctx_unique"),
#    ForeignKey("variant", name="svs_ctx_fk", ondelete="CASCADE", table="variants"),
)


variants_table = Table(
    "variants",
    metadata,
    Column("variant_id", String(255)),
    Column("id", Integer, primary_key=True),
    Column("var_type", String(30)),
    UniqueConstraint("variant_id", name="variants_unique"),
)

variants_annotations_table = Table(
    "variants_annotations",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("hgvsp", String(255)),
    Column("polyphen", String(255)),
    Column("sift", String(255)),
    Column("variant_transcript", Integer, ForeignKey("variants_transcripts.id", ondelete="CASCADE")),
    UniqueConstraint("variant_transcript", name="variants_annotations_unique"),
)
# ... (Previous code)

# Define the "variants_consequences" table
variants_consequences_table = Table(
    "variants_consequences",
    metadata,
    Column("id", Integer, primary_key=True),
    Column("severity", Integer, ForeignKey("severities.id", ondelete="CASCADE")),
    Column("variant_transcript", Integer, ForeignKey("variants_transcripts.id", ondelete="CASCADE")),
)


variants_transcripts_table = Table(
    "variants_transcripts",
    metadata,
    Column("transcript", Integer, ForeignKey("transcripts.id", ondelete="CASCADE")),
    Column("id", Integer, primary_key=True),
    Column("variant", Integer, ForeignKey("variants.id", ondelete="CASCADE")),
    Column("hgvsc", String(255)),
    UniqueConstraint("transcript", "variant", name="variants_transcripts_unique"),
)

severities_table = Table(
    "severities",
    metadata,
    Column("severity_number", Integer),
    Column("id", Integer, primary_key=True),
    Column("consequence", String(255)),
    UniqueConstraint("severity_number", name="severities_unique"),
)

metadata.create_all(engine)

severities_sql = "INSERT INTO `severities` (`severity_number`, `id`, `consequence`) VALUES(1, 1, 'transcript_ablation'),(2, 2, 'splice_acceptor_variant'),(3, 3, 'splice_donor_variant'),(4, 4, 'stop_gained'),(5, 5, 'frameshift_variant'),(6, 6, 'stop_lost'),(7, 7, 'start_lost'),(8, 8, 'transcript_amplification'),(9, 9, 'inframe_insertion'),(10, 10, 'inframe_deletion'),(11, 11, 'missense_variant'),(12, 12, 'protein_altering_variant'),(13, 13, 'regulatory_region_ablation'),(14, 14, 'splice_region_variant'),(15, 15, 'incomplete_terminal_codon_variant'),(16, 16, 'start_retained_variant'),(17, 17, 'stop_retained_variant'),(18, 18, 'synonymous_variant'),(19, 19, 'coding_sequence_variant'),(20, 20, 'mature_miRNA_variant'),(21, 21, '5_prime_UTR_variant'),(22, 22, '3_prime_UTR_variant'),(23, 23, 'non_coding_transcript_exon_variant'),(24, 24, 'intron_variant'),(25, 25, 'NMD_transcript_variant'),(26, 26, 'non_coding_transcript_variant'),(27, 27, 'upstream_gene_variant'),(28, 28, 'downstream_gene_variant'),(29, 29, 'TFBS_ablation'),(30, 30, 'TFBS_amplification'),(31, 31, 'TF_binding_site_variant'),(32, 32, 'regulatory_region_amplification'),(33, 33, 'feature_elongation'),(34, 34, 'regulatory_region_variant'),(35, 35, 'feature_truncation'),(36, 36, 'intergenic_variant');"
severities_sql_oracle = "INSERT INTO "
# insert into engine
with engine.connect() as connection:
    connection.execute(text(severities_sql))


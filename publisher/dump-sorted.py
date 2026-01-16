import pandas as pd

from dotenv import load_dotenv
import os
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import (
    create_engine,
    MetaData,
    Table,
    select,
)
from publish_test import *

load_dotenv()

dbConnectionString = os.environ.get("DB")
schema = os.environ.get("SCHEMA")


print("connecting...")
engine = create_engine(dbConnectionString, echo=False, pool_pre_ping=True, pool_recycle=3600)

metadata = MetaData()
metadata.reflect(bind=engine, schema=schema)

Variants = Table("variants", metadata, schema=schema)
Snvs = Table("snvs", metadata, schema=schema)
VariantsTranscripts = Table("variants_transcripts", metadata, schema=schema)
Transcripts = Table("transcripts", metadata, schema=schema)
Genes = Table("genes", metadata, schema=schema)
VariantsAnnotations = Table("variants_annotations", metadata, schema=schema)
VariantsConsequences = Table("variants_consequences", metadata, schema=schema)
Severities = Table("severities", metadata, schema=schema)
IbvlFrequencies = Table("genomic_ibvl_frequencies",metadata, schema=schema)
GnomadFrequencies = Table("genomic_gnomad_frequencies",metadata, schema=schema)
MtIbvlFrequencies = Table("mt_ibvl_frequencies",metadata, schema=schema)
MtGnomadFrequencies = Table("mt_gnomad_frequencies",metadata, schema=schema)
Mts = Table("mts",metadata, schema=schema)



with engine.connect() as connection:
    statement = select(
            Variants.c.variant_id, 
            Snvs.c.type,
            Snvs.c.pos,
            VariantsTranscripts.c.hgvsc, 
            Transcripts.c.transcript_id, 
            Genes.c.short_name,
            VariantsAnnotations.c.hgvsp,
            VariantsConsequences.c.severity,
            IbvlFrequencies.c.af_tot,
            GnomadFrequencies.c.af_tot,
            MtGnomadFrequencies.c.af_hom,
            MtIbvlFrequencies.c.an,
            MtIbvlFrequencies.c.af_hom,
            Mts.c.pos
        ).join(VariantsTranscripts, Variants.c.id == VariantsTranscripts.c.variant, isouter=True
        ).join(Transcripts, VariantsTranscripts.c.transcript == Transcripts.c.id, isouter=True
        ).join(VariantsAnnotations, VariantsTranscripts.c.id == VariantsAnnotations.c.variant_transcript, isouter=True
        ).join(VariantsConsequences, VariantsTranscripts.c.id == VariantsConsequences.c.variant_transcript, isouter=True
        ).join(Genes, Transcripts.c.gene == Genes.c.id, isouter=True
        ).join(Snvs, Snvs.c.variant == Variants.c.id, isouter=True
        ).join(Mts, Mts.c.variant == Variants.c.id, isouter=True
        ).join(IbvlFrequencies, IbvlFrequencies.c.variant == Variants.c.id, isouter=True
        ).join(GnomadFrequencies, GnomadFrequencies.c.variant == Variants.c.id, isouter=True
        ).join(MtGnomadFrequencies, MtGnomadFrequencies.c.variant == Variants.c.id, isouter=True
        ).join(MtIbvlFrequencies, MtIbvlFrequencies.c.variant == Variants.c.id, isouter=True
        
        ).order_by(
            Variants.c.variant_id, 
            Transcripts.c.transcript_id, 
            VariantsTranscripts.c.hgvsc,
            VariantsAnnotations.c.hgvsp,
            VariantsConsequences.c.severity,
            
        )
    print("executing statement...")
    result = connection.execute(statement)
    print("building dataframe...")
    df = pd.DataFrame(result.fetchall())
    tsv_output = df.to_csv(index=False, sep="\t")
    print("saving to tsv...")
    with open("variants-6chr1.tsv", "w") as f:
        f.write(tsv_output + "\n schema:"+schema)
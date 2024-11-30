

#todo add key expression
#todo add depends_on models 
#todo add cache_by_chromosome boolean


model_import_actions = {
    
    "genes": {
        "name": "genes",
        "table": "genes",
        "map_key_expression": lambda row: row.short_name.upper(),
        "tsv_map_key_expression": lambda row: row["short_name"].upper(),
        "pk_lookup_col": ["short_name"],
        "cache_by_chromosome": False,
        "fk_map": {},
        "filters": {"short_name": lambda x: x.upper() if x is not None else None},
    },
    "transcripts": {
        "name": "transcripts",
        "table": "transcripts",
        "map_key_expression": lambda row: row.transcript_id,
        "tsv_map_key_expression": lambda row: row["transcript_id"],
        "pk_lookup_col": ["transcript_id"],
        "cache_by_chromosome": False,
        "fk_map": {"gene": "genes"},
        "filters": {
            "transcript_type": lambda x: x.replace("RefSeq", "R") if x is not None else None,
        }
    },
    "variants": {
        "name": "variants",
        "table": "variants",
        "map_key_expression": lambda row: row.variant_id,
        "tsv_map_key_expression": lambda row: row["variant_id"],
        "pk_lookup_col": ["variant_id"],
        "fk_map": {},
    },
    "variants_transcripts": {
        "name": "variants_transcripts",
        "table": "variants_transcripts",
        "map_key_expression": lambda row: (row.variant_id, row.transcript_id),
        "tsv_map_key_expression": lambda row: (row["variant"], row["transcript"]),
        "pk_lookup_col": [],
        "fk_map": {"transcript": "transcripts", "variant": "variants"},
    },
    "variants_annotations": {
        "name": "variants_annotations",
        "table": "variants_annotations",
        "map_key_expression": lambda row: (row.variant_transcript),
        "tsv_map_key_expression": lambda row: (row["variant"], row["transcript"]),
        "pk_lookup_col": ["variant_transcript"],
        "fk_map": {"variant_transcript": "variants_transcripts"},
        "filters": {
            "hgvsp": lambda x: x.replace("%3D", "=") if x is not None else None
        },
    },
    "variants_consequences": {
        "name": "variants_consequences",
        "table": "variants_consequences",
        "map_key_expression": lambda row: (row.variant_transcript, row.severity), 
        "tsv_map_key_expression": lambda row: (row["variant"], row["transcript"], row["severity"]),
        "pk_lookup_col": ["variant_transcript", "severity"],
        "fk_map": {"variant_transcript": "variants_transcripts"},
    },
    "snvs": {
        "name": "snvs",
        "table": "snvs",
        "map_key_expression": lambda row: row.variant_id,
        "tsv_map_key_expression": lambda row: row["variant"],
        "pk_lookup_col": [],
        "fk_map": {"variant": "variants"},
        "filters": {
            "dbsnp_id": lambda x: x.split("&")[0] if x is not None else None,
            "chr": lambda x: str(x).replace("chr", "") if x is not None else None,
            },
    },
    "mts": {
        "name": "mts",
        "table": "mts",
        "map_key_expression": lambda row: row.variant_id,
        "tsv_map_key_expression": lambda row: row["variant"],
        "pk_lookup_col": [],
        "fk_map": {"variant": "variants"},
    },
    "genomic_ibvl_frequencies": {
        "name": "genomic_ibvl_frequencies",
        "table": "genomic_ibvl_frequencies",
        "map_key_expression": lambda row: row.variant_id,
        "tsv_map_key_expression": lambda row: row["variant"],
        "pk_lookup_col": [],
        "fk_map": {"variant": "variants"},
    },
    "genomic_gnomad_frequencies": {
        "name": "genomic_gnomad_frequencies",
        "table": "genomic_gnomad_frequencies",
        "map_key_expression": lambda row: row.variant_id,
        "tsv_map_key_expression": lambda row: row["variant"],
        "pk_lookup_col": [],
        "fk_map": {"variant": "variants"},
    },
    "mt_ibvl_frequencies": {
        "name": "mt_ibvl_frequencies",
        "table": "mt_ibvl_frequencies",
        "map_key_expression": lambda row: row.variant_id,
        "tsv_map_key_expression": lambda row: row["variant"],
        "pk_lookup_col": [],
        "fk_map": {"variant": "variants"},
    },
    "mt_gnomad_frequencies": {
        "name": "mt_gnomad_frequencies",
        "table": "mt_gnomad_frequencies",
        "map_key_expression": lambda row: row.variant_id,
        "tsv_map_key_expression": lambda row: row["variant"],
        "pk_lookup_col": [],
        "fk_map": {"variant": "variants"},
    },
    
    
    
    
    
    
    
#   "svs": {
#       "name": "svs",
#       "table": "svs",
#       "pk_lookup_col": None,
#       "fk_map": {"variant": "variants"},
#   },
#    "sv_consequences": {
#        "name": "sv_consequences",
#        "table": "sv_consequences",
#        "pk_lookup_col": None,
#        "fk_map": {"gene": "genes", "variant": "variants"},
#    },
#   "svs_ctx": {
#       "name": "svs_ctx",
#       "table": "svs_ctx",
#       "pk_lookup_col": None,
#       "fk_map": {"variant": "variants"},
#   },
#   "str": {
#       "name": "str",
#       "table": "str",
#       "pk_lookup_col": None,
#       "fk_map": {"variant": "variants"},
#   },
    
    
    
    
    
}

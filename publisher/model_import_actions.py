
model_import_actions = {
    "genes": {
        "name": "genes",
        "table": "genes",
        "pk_lookup_col": "short_name",
        "fk_map": {},
        "filters": {"short_name": lambda x: x.upper() if x is not None else None},
    },
    "transcripts": {
        "name": "transcripts",
        "table": "transcripts",
        "pk_lookup_col": "transcript_id",
        "fk_map": {"gene": "genes"},
        "filters": {
            "transcript_type": lambda x: x.replace("RefSeq", "R") if x is not None else None,
        }
    },
    "variants": {
        "name": "variants",
        "table": "variants",
        "pk_lookup_col": "variant_id",
        "fk_map": {},
    },
    "variants_transcripts": {
        "name": "variants_transcripts",
        "table": "variants_transcripts",
        "pk_lookup_col": ["transcript", "variant"],
        "fk_map": {"transcript": "transcripts", "variant": "variants"},
    },
    "variants_annotations": {
        "name": "variants_annotations",
        "table": "variants_annotations",
        "pk_lookup_col": None,
        "fk_map": {"DO_COMPOUND_FK": "for variants_transcripts"},
        "filters": {
            "hgvsp": lambda x: x.replace("%3D", "=") if x is not None else None
        },
    },
    #    "severities":{
    #        "name":"severities",
    #        "pk_lookup_col": None,
    #        "fk_map": {},
    #        },
    "variants_consequences": {
        "name": "variants_consequences",
        "table": "variants_consequences",
        "pk_lookup_col": None,
        "fk_map": {"DO_COMPOUND_FK": "for variants_transcripts"},
    },
    "sv_consequences": {
        "name": "sv_consequences",
        "table": "sv_consequences",
        "pk_lookup_col": None,
        "fk_map": {"gene": "genes", "variant": "variants"},
    },
    "snvs": {
        "name": "snvs",
        "table": "snvs",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
        "filters": {"dbsnp_id": lambda x: x.split("&")[0] if x is not None else None},
    },
    "svs": {
        "name": "svs",
        "table": "svs",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
    "svs_ctx": {
        "name": "svs_ctx",
        "table": "svs_ctx",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
    "str": {
        "name": "str",
        "table": "str",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
    "mts": {
        "name": "mts",
        "table": "mts",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
    "genomic_ibvl_frequencies": {
        "name": "genomic_ibvl_frequencies",
        "table": "genomic_ibvl_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
    "genomic_gnomad_frequencies": {
        "name": "genomic_gnomad_frequencies",
        "table": "genomic_gnomad_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
    "mt_ibvl_frequencies": {
        "name": "mt_ibvl_frequencies",
        "table": "mt_ibvl_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
    "mt_gnomad_frequencies": {
        "name": "mt_gnomad_frequencies",
        "table": "mt_gnomad_frequencies",
        "pk_lookup_col": None,
        "fk_map": {"variant": "variants"},
    },
}

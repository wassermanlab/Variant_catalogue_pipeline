from do_import import *
from import_utils import *

from model_import_actions import model_import_actions

def test(engine,job_dir = None):
        
    NUM_TEST_ROWS = int(os.environ.get("NUM_TEST_ROWS", 200))
#    if job_dir is not None:
#        setup_loggers(job_dir)
        
    def output(msg):
        if job_dir is not None:
            log_output(msg)
        else:
            print(msg)
            
    output(f"testing {os.environ.get('SCHEMA_NAME')}, assembly {os.environ.get('SET_VAR_ASSEMBLY')}")
    output(f"num test rows: {NUM_TEST_ROWS}")
    output(f"tsv folder is {os.environ.get('PIPELINE_OUTPUT_PATH', '/fixtures')}")
    
    
    if isinstance(schema, str) and len(schema) > 0:
        metadata.reflect(bind=engine, schema=schema)
    else:
        metadata.reflect(bind=engine)
        
    did_pass = True
    num_pass = 0
    def fail(msg, msg2 = None):
        nonlocal did_pass
        did_pass = False
        if (msg2 is not None):
            output(F"❌ {msg}\n    {msg2}")
        else:
            output(F"❌ {msg}")
#        cleanup(None, None)
#        exit()
        
    def get_random_tsv_file(model_folder):
        path = os.path.join(rootDir, model_folder)
        files = [os.path.join(path,f) for f in os.listdir(path) if not f.startswith(".")]
        file = np.random.choice(files, 1)[0]
        return inspectTSV(file), file
        
    def get_random_tsv_rows(model_folder, n):
        random_file_info, random_file = get_random_tsv_file(model_folder)
        
        types = model_import_actions[model_folder].get("tsv_types") or {}
        
        df = readTSV(random_file, random_file_info, dtype=types)
        df_rowcount = len(df)
        if df_rowcount < n:
            n = df_rowcount
        return df.sample(n)
    
    def testmodel(model, select_tables, join_fn, where_fn, data_cols, checks=[], skip=False, null_fk_case=None):
        output(f"testing {model}...")
        
        if skip:
            return
        nonlocal num_pass
        with engine.connect() as connection:
            tsv_rows = get_random_tsv_rows(model, NUM_TEST_ROWS)
            for _, row in tsv_rows.iterrows():
                tsv_row = row.to_dict()
                
                tsv_row_filters = model_import_actions[model].get("filters") or {}
                for col, filter in tsv_row_filters.items():
                    tsv_row[col] = filter(tsv_row[col])
                    
                statement = select(*select_tables)
                statement = join_fn(statement)
                statement = where_fn(statement, tsv_row)
                
                db_rows = connection.execute(statement).fetchall()
                
                if len(db_rows) == 0 and null_fk_case is not None:
                    db_rows = connection.execute(null_fk_case(tsv_row)).fetchall()
                
                if len(db_rows) == 0:
                    fail(f"{model} row not found", tsv_row)
                elif len(db_rows) > 1:
                    fail(f"{model} multiple rows found", tsv_row)
                else:
                    row_dict = db_rows[0]._mapping
                    for col in data_cols:
                        db_val = row_dict[col]
                        tsv_val = tsv_row[col]
                        
                        if isinstance(db_val, Decimal) or isinstance(db_val, float) or isinstance(db_val, int):
                            tsv_float = float(tsv_val)
                            if math.isclose(db_val, tsv_float, rel_tol=1e-7, abs_tol=1e-7):
                                continue
                            else:
                                fail(f"{model} column {col} numerical mismatch: db's {db_val} != tsv's {tsv_val}", tsv_row)
                            
                        elif isinstance(db_val, str):
                            if db_val == tsv_val:
                                continue
                            else:
                                fail(f"{model} column {col} string mismatch: db's {db_val} != tsv's {tsv_val}", tsv_row)
                        elif db_val is None and tsv_val is None:
                            continue
                        else:
                            print(f"unhandled type {type(db_val)}")
                            quit()
                    for check in checks:
                        checkFailMsg = check(row_dict, tsv_row)
                        if checkFailMsg is not None:
                            fail(checkFailMsg, tsv_row)
        num_pass += 1
        
    genes = get_table("genes")
    
    variants = get_table("variants")      
    variants_transcripts = get_table("variants_transcripts")
    transcripts = get_table("transcripts")
    genomic_ibvl_frequencies = get_table("genomic_ibvl_frequencies")
    mt_ibvl_frequencies = get_table("mt_ibvl_frequencies")
    genomic_gnomad_frequencies = get_table("genomic_gnomad_frequencies")
    mt_gnomad_frequencies = get_table("mt_gnomad_frequencies")
    mts = get_table("mts")
    snvs = get_table("snvs")
    variants_annotations = get_table("variants_annotations")
    variants_consequences = get_table("variants_consequences")
    
    testmodel("variants", 
            [variants], 
            join_fn=lambda stmt: stmt, 
            where_fn=lambda stmt, tsv_r: stmt.where(
                variants.c.variant_id == tsv_r["variant_id"], 
                variants.c.assembly == set_var_assembly), 
            data_cols=['var_type']
            )
    
    def checkTranscriptGene(db_row, tsv_row):
        tsv_gene = tsv_row["gene"]
        db_gene = db_row.get("short_name")
        if db_gene == tsv_gene:
            return None
        else:
            return f"Gene_id mismatch: db {db_gene} != tsv {tsv_gene}"

    testmodel("transcripts", 
            [transcripts, genes], 
            join_fn=lambda stmt: stmt.join(genes, transcripts.c.gene == genes.c.id), 
            where_fn=lambda stmt, tsv: stmt.where(transcripts.c.transcript_id == tsv["transcript_id"]),
            data_cols = ['transcript_type', 'tsl'],
            checks = [
                checkTranscriptGene
            ],
            null_fk_case = lambda tsv: select(transcripts).where(transcripts.c.transcript_id == tsv["transcript_id"])
            )
    testmodel("variants_transcripts", 
            [variants_transcripts, variants, transcripts], 
            join_fn=lambda stmt: stmt.join(variants, variants_transcripts.c.variant == variants.c.id).join(transcripts, variants_transcripts.c.transcript == transcripts.c.id), 
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"], 
                                                    variants.c.assembly == set_var_assembly, 
                                                    transcripts.c.transcript_id == tsv_r["transcript"]), 
            data_cols=['hgvsc']
            )
    
    testmodel("snvs", 
            [snvs, variants], 
            join_fn=lambda stmt: stmt.join(variants, snvs.c.variant == variants.c.id), 
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"], 
                                                    variants.c.assembly == set_var_assembly), 
            data_cols=['type', 'length', 'chr', 'pos', 'ref', 'alt', 'cadd_score', 'cadd_intr', 'dbsnp_id', 'dbsnp_url', 'ucsc_url', 'ensembl_url', 'clinvar_url', 'gnomad_url', 'clinvar_vcv', 'splice_ai'],
            )
    
    testmodel("mts",
            [mts, variants],
            join_fn=lambda stmt: stmt.join(variants, mts.c.variant == variants.c.id),
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"],
                                                    variants.c.assembly == set_var_assembly),
            data_cols=['pos', 'ref', 'alt', 'ucsc_url', 'mitomap_url', 'gnomad_url', 'dbsnp_id', 'dbsnp_url', 'clinvar_url', 'clinvar_vcv'],
            )
    
    testmodel("genomic_ibvl_frequencies", 
            [genomic_ibvl_frequencies, variants], 
            join_fn=lambda stmt: stmt.join(variants, genomic_ibvl_frequencies.c.variant == variants.c.id), 
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"], 
                                                    variants.c.assembly == set_var_assembly), 
            data_cols=['af_tot', 'af_xx', 'af_xy', 'ac_tot', 'ac_xx', 'ac_xy', 'an_tot', 'an_xx', 'an_xy', 'hom_tot', 'hom_xx', 'hom_xy', 'quality']
            )
    
    
    testmodel("genomic_gnomad_frequencies",
            [genomic_gnomad_frequencies, variants],
            join_fn=lambda stmt: stmt.join(variants, genomic_gnomad_frequencies.c.variant == variants.c.id),
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"],
                                                    variants.c.assembly == set_var_assembly),
            data_cols=['af_tot', 'ac_tot', 'an_tot', 'hom_tot']
            )
    
    testmodel("mt_ibvl_frequencies",
            [mt_ibvl_frequencies, variants],
            join_fn=lambda stmt: stmt.join(variants, mt_ibvl_frequencies.c.variant == variants.c.id),
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"],
                                                    variants.c.assembly == set_var_assembly),
            #variant	an	ac_hom	ac_het	af_hom	af_het	hl_hist	max_hl
            data_cols = ['an', 'ac_hom', 'ac_het', 'af_hom', 'af_het', 'hl_hist', 'max_hl']
            )

    testmodel("mt_gnomad_frequencies",
            [mt_gnomad_frequencies, variants],
            join_fn=lambda stmt: stmt.join(variants, mt_gnomad_frequencies.c.variant == variants.c.id),
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"],
                                                    variants.c.assembly == set_var_assembly),
            # an	ac_hom	ac_het	af_hom	af_het	max_hl
            data_cols = ['an', 'ac_hom', 'ac_het', 'af_hom', 'af_het', 'max_hl']
            )

    with engine.connect() as connection:
        def get_variant_transcript(tsv_row):
            statement = select(variants_transcripts, variants, transcripts)
            statement = statement.join(variants, variants_transcripts.c.variant == variants.c.id).join(transcripts, variants_transcripts.c.transcript == transcripts.c.id)
            statement = statement.where(variants.c.variant_id == tsv_row["variant"], variants.c.assembly == set_var_assembly, transcripts.c.transcript_id == tsv_row["transcript"])
            db_rows = connection.execute(statement).fetchall()
            
            if len(db_rows) == 0:
                fail(f"(variant_transcript) transcript not found", tsv_row)
                return None
            elif len(db_rows) > 1:
                fail(f"(variant_transcript) multiple transcripts found", tsv_row)
                return None
            else:
                return db_rows[0]._mapping
        
        # testing variants_consequences
        tsv_rows = get_random_tsv_rows("variants_consequences", NUM_TEST_ROWS)
        for _, row in tsv_rows.iterrows():
            tsv_row = row.to_dict()
            
            transcript = get_variant_transcript(tsv_row)
            
            if transcript is None:
                continue
            
            tsv_row_filters = model_import_actions["variants_consequences"].get("filters") or {}
            for col, filter in tsv_row_filters.items():
                tsv_row[col] = filter(tsv_row[col])
                
            statement = select(variants_consequences)
            statement = statement.where(variants_consequences.c.variant_transcript == transcript["id"])
            
            db_rows = connection.execute(statement).fetchall()
            if len(db_rows) == 0:
                fail(f"no variant consequences in db matching", tsv_row)
            
            found = False
            for db_row in db_rows:
                if db_row._mapping["severity"] == tsv_row["severity"]:
                    found = True
                    break
            if not found:
                fail(f"all matching variant consequences have wrong severity", tsv_row)
        
        # testing variants_annotations
        tsv_rows = get_random_tsv_rows("variants_annotations", NUM_TEST_ROWS)
        for _, row in tsv_rows.iterrows():
            tsv_row = row.to_dict()
            
            transcript = get_variant_transcript(tsv_row)
            
            if transcript is None:
                continue
            
            tsv_row_filters = model_import_actions["variants_annotations"].get("filters") or {}
            for col, filter in tsv_row_filters.items():
                tsv_row[col] = filter(tsv_row[col])
                
            statement = select(variants_annotations)
            statement = statement.where(variants_annotations.c.variant_transcript == transcript["id"])
            
            db_rows = connection.execute(statement).fetchall()
            if len(db_rows) == 0:
                fail(f"No variant annotations in db matching:", tsv_row)
            elif len(db_rows) > 1:
                fail(f"Multiple variant annotations found in db matching:", tsv_row)
            else:
                db_row_dict = db_rows[0]._mapping
                ##hgvsp	sift	polyphen
                for col in ['hgvsp', 'sift', 'polyphen']:
                    if db_row_dict[col] != tsv_row[col]:
                        fail(f"variant annotation mismatch: db's {db_row_dict[col]} != tsv's {tsv_row[col]}", tsv_row)

    if (did_pass):
        output("✅✅✅ \n all tests passed\n✅✅✅\n")
        
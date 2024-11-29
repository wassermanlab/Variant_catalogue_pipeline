import json
import math
from decimal import Decimal
from natsort import natsorted
from sqlalchemy import (
    MetaData,
    Table,
    Integer,
    String,
    select,
    func,
    Float,
)
from sqlalchemy.exc import DataError, IntegrityError, ProgrammingError
import pandas as pd
import numpy as np
import signal
import sys
import os
from datetime import datetime
from sqlalchemy.orm import sessionmaker
from dotenv import load_dotenv

from import_utils import *
from model_import_actions import model_import_actions


load_dotenv()

# get command line arguments
rootDir = os.environ.get("PIPELINE_OUTPUT_PATH")
dir_containing_jobs = "jobs"
current_dir = os.path.dirname(os.path.realpath(__file__))
dir_containing_jobs = os.path.join(current_dir, dir_containing_jobs)
print(f"dir containing jobs: {dir_containing_jobs}. absolute: {os.path.abspath(dir_containing_jobs)}")
chunk_size = int(os.environ.get("CHUNK_SIZE"))
# verbose = os.environ.get("VERBOSE") == "true"
dbConnectionString = os.environ.get("DB")
isDevelopment = os.environ.get("ENVIRONMENT") != "production"
schema = os.environ.get("SCHEMA_NAME")
dry_run = os.environ.get("DRY_RUN") == "true"
update = os.environ.get("UPDATE") == "true"
set_var_assembly = os.environ.get("SET_VAR_ASSEMBLY", None)

start_at_model = (
    os.environ.get("START_AT_MODEL") if os.environ.get("START_AT_MODEL") != "" else None
)
start_at_file = (
    os.environ.get("START_AT_FILE") if os.environ.get("START_AT_FILE") != "" else None
)

jobs_dir = os.path.abspath(dir_containing_jobs)

db_row_counts = {"before": {}, "after": {}}

if rootDir == None:
    rootDir = os.path.join(current_dir, "fixtures")
    print("No root directory specified. default ", rootDir)
#    exit()

print("importing", rootDir)

engine = None

metadata = MetaData()

job_dir = ""
# map of import functions

depends_on_maps = {}
current_model_existing_map = {}

last_chromosome = None
current_chromosome = None

def separate_cache_by_chromosome(action):
    return action["name"] in ["variants","variants_transcripts", "variants_annotations", "variants_consequences"]


def populate_maps(action, chromosome=None):
    global current_model_existing_map, depends_on_maps
    model = action["name"]
    
    def make_existing_map(modelName, chromosome=None):
        table = get_table(modelName)
        model_action = model_import_actions[modelName]

        with engine.connect() as connection:
            if (modelName == "variants_transcripts"):
                variants = get_table("variants")
                transcripts = get_table("transcripts")
                statement = select(
                    table.c["id"],
                    variants.c["variant_id", "assembly"], 
                    transcripts.c["transcript_id"]
                ).join(variants).join(transcripts).where(variants.c.assembly == set_var_assembly)
                if chromosome is not None:
                    statement = statement.where(variants.c.variant_id.startswith(f"{chromosome}_"))
                
            else:
                cols = [table.c[col] for col in ["id"] + model_action["pk_lookup_col"] ]
                
                if "variant" in model_action["fk_map"]:
                    variants = get_table("variants")
                    cols.append(variants.c["variant_id", "assembly"])
                    statement = select(*cols).join(variants).where(variants.c.assembly == set_var_assembly)
                    if chromosome is not None:
                        statement = statement.where(variants.c.variant_id.startswith(f"{chromosome}_"))
                elif modelName == "variants":
                    cols.append(table.c["assembly"])
                    statement = select(*cols).where(table.c.assembly == set_var_assembly)
                    if chromosome is not None:
                        statement = statement.where(table.c.variant_id.startswith(f"{chromosome}_"))
                else:
                    statement = select(*cols)
                    
            result = connection.execute(statement)
            return {
                model_action["map_key_expression"](row): row.id for row in result
            }
        
    if model in ["variants_annotations", "variants_consequences"]:
        depends_on_maps["variants_transcripts"] = make_existing_map("variants_transcripts", chromosome)
        reversed = {v: k for k, v in depends_on_maps["variants_transcripts"].items()}
        tenative_existing_map = make_existing_map(model, chromosome)
        current_model_existing_map = {reversed.get(k): v for k, v in tenative_existing_map.items() if reversed.get(k) is not None }
        
    else:
        referenced_models = action.get("fk_map").values()
        current_model_existing_map = make_existing_map(action["name"], chromosome)
        for model in referenced_models:

            depends_on_maps[model] = make_existing_map(model,chromosome)

            log_output(
                "loaded map for "
                + model
                + ". number of records: "
                + str(len(depends_on_maps[model]))
            )
    log_output(f"done populating maps for {action['name']} chromosome {chromosome}")
#    existing_json = json.dumps(current_model_existing_map)
#    existing_dependson_json = json.dumps(depends_on_maps)
#    log_output(f"existing map for {action['name']} chromosome {chromosome}: {existing_json}")
#    log_output(f"depends on maps for {action['name']} chromosome {chromosome}: {existing_dependson_json}")

#    except FileNotFoundError:
#        pk_maps[modelName] = {}
#        pass


def append_to_map(modelName, key, value):
    global depends_on_maps
    if modelName not in depends_on_maps:
        log_error(f"trying to append to map but was not populated ({modelName} {key} {value})")
    if key not in depends_on_maps[modelName]:
        depends_on_maps[modelName][key] = value


def persist_and_unload_maps():
    depends_on_maps.clear()


def get_table(model):
    
    if model not in model_import_actions:
        table_name = model
    else:
        table_name = model_import_actions[model]["table"]
        
    if isinstance(schema, str) and len(schema) > 0:
        return Table(table_name, metadata, schema=schema)
    else:
        return Table(table_name, metadata, autoload_with=engine)

def import_file(file, file_info, action):
    global current_chromosome,last_chromosome
    model = action.get("name")
    fk_map = action.get("fk_map")
    pk_lookup_col = action.get("pk_lookup_col")
    filters = action.get("filters") or {}
    
#    cache_group = file_info('cache_group')
#    print(f"cache group is {file_info['cache_group']}")

    successCount = 0
    failCount = 0
    duplicateCount = 0
    missingRefCount = 0
    updatedCount = 0
    successful_chunks = 0
    fail_chunks = 0
    
    table = get_table(model)
    types_dict = {}
    for column in table.columns:
        # convert sql types to pandas types
        if isinstance(column.type, Integer):
            types_dict[column.name] = "Int64"
        elif isinstance(column.type, String):
            types_dict[column.name] = "str"
        elif isinstance(column.type, Float):
            types_dict[column.name] = "float64"
        else:
            pass

    df = readTSV(file, file_info, dtype=types_dict)
    

    data_insert_list = []
    data_update_list = []
    for _, row in df.iterrows():
        data = row.to_dict()
        
        if (model == "variants" and set_var_assembly is not None):
            data["assembly"] = set_var_assembly
        
        if separate_cache_by_chromosome(action):
        
            current_chromosome = data.get("variant",data.get("variant_id")).split("_")[0]
            if current_chromosome != last_chromosome:
                last_chromosome = current_chromosome
                persist_and_unload_maps()
                populate_maps(action, current_chromosome)
        
        

        skip = False
        record_map_key = action.get("tsv_map_key_expression")(data)
        for col, filter in filters.items():
            data[col] = filter(data[col])
        for depended_model_col, depended_model in fk_map.items():
            depended_map_key = None
            resolved_pk = None
            
            if model in ["variants_annotations", "variants_consequences"]:
                depended_map_key = (data["variant"], data["transcript"])
            else:
                
                if isinstance(data[depended_model_col], str):
                    depended_map_key = data[depended_model_col]
                if depended_model == "genes":
                    depended_map_key = depended_map_key.upper()
            if depended_map_key == "NA":
                    data[depended_model_col] = None
            else:
                resolved_pk = depends_on_maps.get(depended_model).get( depended_map_key)
            if resolved_pk is not None:
                data[depended_model_col] = resolved_pk
            else:
                log_data_issue(
                    f"Missing {depended_model_col} {depended_map_key} in {depended_model}"
                    if depended_model_col is not None
                    else (
                        f"None {depended_map_key}"
                        if depended_map_key is not None
                        else (
                            f"None referenced from {model}"
                            if model is not None
                            else "None"
                        )
                    ),
                    model,
                )
                log_data_issue(data, model)
                missingRefCount += 1
                skip = True
        if skip:
            continue
#        # this block replaces None with empty string for string columns (good for django)
#        for col in data:
#            if (
#                col in table.columns
#                and isinstance(table.columns[col].type, String)
#                and data[col] is None
#            ):
#                data[col] = ""
#        # this block fills missing columns with None or empty string, depending on the column type
#        for table_col in table.columns:
#            if table_col.name not in data:
#                if isinstance(table_col.type, String):
#                    data[table_col.name] = ""
#                    log_data_issue("filled missing col " + table_col.name + " with empty string", name)
#                else:
#                    data[table_col.name] = None
#                    log_data_issue("filled missing col " + table_col.name + " with None", name)

        if record_map_key is None:
            print(f"record_map_key is None for {model} {data}")
            quit()
        existing_id = current_model_existing_map.get(record_map_key)
        if existing_id is not None:
            # record is already in the DB
            if update:
                data["id"] = existing_id
                data_update_list.append(data)
            else:
                duplicateCount += 1
                successCount += 1
                log_data_issue("Duplicate " + str(record_map_key), model)
                continue
        else:
            # the record is NOT in the db, so add it
            data_insert_list.append(data)
#        data_list.append(data)

    # dispose of df to save ram
    del df
    with engine.connect() as connection:

        def rowOperation(row, updating=False):
            nonlocal successCount, failCount, duplicateCount, updatedCount, successful_chunks, fail_chunks
            
            did_succeed = False
            try:
                if updating:
                    connection.execute(
                        table.update().where(table.c.id == row["id"]), row
                    )
                    connection.commit()
                    successCount += 1
                    updatedCount += 1
                    did_succeed = True
                else:
                    connection.execute(table.insert(), row)
                    connection.commit()
                    successCount += 1
                    did_succeed = True

            except DataError as e:
                log_data_issue(e, model)
                failCount += 1
            except IntegrityError as e:
                msg = str(e)
                if "Duplicate" in msg or "ORA-00001" in msg:
                    duplicateCount += 1
                    successCount += 1
#                    log_data_issue(e, model)
                    if updating:
                        updatedCount += 1
                else:
                    failCount += 1
                    log_data_issue(e, model)
            except Exception as e:

                log_data_issue(e, model)
                failCount += 1
            
            if not did_succeed:
                connection.rollback()
                
        def chunkOperation(chunk, updating=False):
            nonlocal successCount, failCount, duplicateCount, updatedCount, successful_chunks, fail_chunks
            if dry_run:
                successCount += len(chunk)
                return
            try:
                if updating:
                    for row in chunk:
                        rowOperation(row, updating)
                else: 
                    connection.execute(table.insert(), chunk)
                    # commit
                    connection.commit()
                    # chunk worked
                    successful_chunks += 1
                    successCount += len(chunk)
            except Exception as e:
                #                print(e)
                connection.rollback()
                fail_chunks += 1
                for row in chunk:
                    rowOperation(row, updating)
                        
        for chunk in chunks(data_insert_list, chunk_size):
            chunkOperation(chunk)
        for chunk in chunks(data_update_list, chunk_size):
            chunkOperation(chunk, updating=True)
        

    return {
        "success": successCount,
        "fail": failCount,
        "missingRef": missingRefCount,
        "duplicate": duplicateCount,
        "updated": updatedCount,
        "successful_chunks": successful_chunks,
        "fail_chunks": fail_chunks,
        "rowcount": file_info["total_rows"],
    }


def cleanup(sig, frame):
    global engine, depends_on_maps, metadata
    print("terminating, cleaning up ...")
    persist_and_unload_maps()
    engine.dispose()
    # garbage collect
    del depends_on_maps
    del metadata
    print("done")
    sys.exit(0)


signal.signal(signal.SIGINT, cleanup)


def start(db_engine):

    arrived_at_start_model = False
    arrived_at_start_file = False
    global job_dir, engine, schema, current_chromosome, last_chromosome
    engine = db_engine

    if isinstance(schema, str) and len(schema) > 0:
        metadata.reflect(bind=engine, schema=schema)
    else:
        metadata.reflect(bind=engine)
    Session = sessionmaker(bind=engine)

    os.makedirs(jobs_dir, exist_ok=True)
    os.makedirs(os.path.join(jobs_dir, "1"), exist_ok=True)

    without_hidden = [f for f in os.listdir(jobs_dir) if not f.startswith(".")]
    last_job = int(natsorted(without_hidden)[-1])
    if os.listdir(os.path.join(jobs_dir, str(last_job))) == []:
        job_dir = os.path.join(jobs_dir, str(last_job))
    else:
        job_dir = os.path.join(jobs_dir, str(last_job + 1))
    os.makedirs(job_dir, exist_ok=True)
    os.chmod(job_dir, 0o777)  # Set read and write permissions for the directory
    setup_loggers(job_dir)

    now = datetime.now()
    counts = {}
    counts["success"] = 0
    counts["fail"] = 0
    counts["missingRef"] = 0
    counts["duplicate"] = 0
    counts["updated"] = 0
    counts["successful_chunks"] = 0
    counts["fail_chunks"] = 0
    counts["rowcount"] = 0
    
    
    severitiesFile = os.path.join(current_dir, "severities.tsv")
    
    file_info = inspectTSV(severitiesFile)
    log_output(
        "\nimporting severities "
        + " ("
        + severitiesFile.split("/")[-1]
        + "). Expecting "
        + str(file_info["total_rows"])
        + " rows..."
    )
    # log_output(targetFile)
    if file_info["total_rows"] == 0:
        log_output("Skipping empty file")
    import_file(
        severitiesFile,
        file_info,
        {"name":"severities", "fk_map":{}, "pk_lookup_col":None, "tsv_map_key_expression": lambda row: row["severity_number"], "filters":{}},
    )
    log_output("done importing severities")


    for modelName, action_info in model_import_actions.items():
        model_counts = {}
        model_counts["success"] = 0
        model_counts["fail"] = 0
        model_counts["missingRef"] = 0
        model_counts["duplicate"] = 0
        model_counts["updated"] = 0
        model_counts["successful_chunks"] = 0
        model_counts["fail_chunks"] = 0
        model_counts["rowcount"] = 0
        model_directory = os.path.join(rootDir, modelName)
        
        with engine.connect() as connection:
            table = get_table(modelName)
            num_rows = connection.execute(select(func.count()).select_from(table)).scalar()
            db_row_counts["before"][modelName] = num_rows

        if (
            isinstance(start_at_model, str)
            and modelName != start_at_model
            and not arrived_at_start_model
        ):
            log_output("Skipping " + modelName + ", until " + start_at_model)
            continue

        if isinstance(start_at_model, str) and modelName == start_at_model:
            arrived_at_start_model = True

        ######### added in v2. handles case when the pipeline output directory
        # is not a directory of directories of tsv files (ie, per chromosome), but a single directory of tsv files,
        # with one tsv file per model
        large_model_file_tsv = os.path.join(rootDir, modelName + ".tsv")
        large_model_file_tsv_exists = os.path.isfile(large_model_file_tsv)

        if action_info.get("skip") or not os.path.isdir(model_directory):
            if large_model_file_tsv_exists:
                log_output("using large model tsv file " + large_model_file_tsv)
            else:
                log_output(
                    "Skipping " + modelName + " (expected dir: " + model_directory + ")"
                )
                continue
        
        if separate_cache_by_chromosome(action_info):
            pass #populate maps instead happens per 1 chromosome
        else:
            populate_maps(action_info)
        modelNow = datetime.now()

        if large_model_file_tsv_exists:
            sorted_files = [modelName + ".tsv"]
        else:
            sorted_files = natsorted(
                [f for f in os.listdir(model_directory) if not f.startswith(".")],
            )

        for file in sorted_files:
            if file.endswith(".tsv"):

                if (
                    isinstance(start_at_file, str)
                    and file != start_at_file
                    and not arrived_at_start_file
                ):
                    log_output("Skipping " + file + ", until " + start_at_file)
                    continue
                if isinstance(start_at_file, str) and file == start_at_file:
                    arrived_at_start_file = True

                ######## added in v2. large tsv file handling as explained above.
                if large_model_file_tsv_exists:
                    targetFile = large_model_file_tsv
                else:
                    targetFile = model_directory + "/" + file
                
                file_info = inspectTSV(targetFile)
                log_output(
                    "\nimporting "
                    + modelName
                    + " ("
                    + targetFile.split("/")[-1]
                    + "). Expecting "
                    + str(file_info["total_rows"])
                    + " rows..."
                )
                # log_output(targetFile)
                if file_info["total_rows"] == 0:
                    log_output("Skipping empty file")
                    continue
                results = import_file(
                    targetFile,
                    file_info,
                    action_info,
                )
                if results["success"] == 0:
                    log_output("No rows were imported.")

                for key in [
                    "success",
                    "fail",
                    "missingRef",
                    "duplicate",
                    "updated",
                    "successful_chunks",
                    "fail_chunks",
                    "rowcount"
                ]:
                    model_counts[key] += results[key]
                    counts[key] += results[key]

                report_counts(results)

        log_output(
            "\nFinished importing "
            + modelName
            + ". Took this much time: "
            + str(datetime.now() - modelNow)
        )
        current_chromosome = None
        last_chromosome = None
        report_counts(model_counts)
        this_model_index = list(model_import_actions.keys()).index(modelName)
        if this_model_index + 1 < len(model_import_actions.keys()):
            leftover_models = list(model_import_actions.keys())[this_model_index + 1 :]
            log_output("\nmodels left still: " + str(leftover_models) + "\n")

        persist_and_unload_maps()
        
        
        with engine.connect() as connection:
            table = get_table(modelName)
            num_rows = connection.execute(select(func.count()).select_from(table)).scalar()
            db_row_counts["after"][modelName] = num_rows
            
    log_output(f"finished importing IBVL. Time Taken: {str(datetime.now() - now)}. was job {job_dir}")
    report_counts(counts)
    
    delta = {}
    for beforeafter, counts in db_row_counts.items():
        for modelName, count in counts.items():
            if beforeafter == "before":
                delta[modelName] = count
            else:
                delta[modelName] = count - delta[modelName]
                log_output(f"DB row count {modelName} {count} ( grew by {delta[modelName]})")
    
    log_output("testing...")
    
    did_pass = True
    def fail(msg):
        did_pass = False
        log_error(F"❌❌❌\n{msg} \n❌❌❌\n\n")
#        cleanup(None, None)
#        exit()
        
    def get_random_tsv_file(model_folder):
        path = os.path.join(rootDir, model_folder)
        files = natsorted(
            [os.path.join(path,f) for f in os.listdir(path) if not f.startswith(".")],
        )
        file = np.random.choice(files, 1)[0]
        return inspectTSV(file), file
        
    def get_random_tsv_rows(model_folder, n):
        random_file_info, random_file = get_random_tsv_file(model_folder)
        df = readTSV(random_file, random_file_info)
        df_rowcount = len(df)
        if df_rowcount < n:
            n = df_rowcount
        return df.sample(n)
    
    
    

    
    def testmodel(model, select_tables, join_fn, where_fn, data_cols, checks=[]):
        with engine.connect() as connection:
            table = get_table(model)
            tsv_rows = get_random_tsv_rows(model, 10)
            for _, row in tsv_rows.iterrows():
                tsv_row = row.to_dict()
                
                tsv_row_filters = model_import_actions[model].get("filters") or {}
                for col, filter in tsv_row_filters.items():
                    tsv_row[col] = filter(tsv_row[col])
                    
                statement = select(*select_tables)
                statement = join_fn(statement)
                statement = where_fn(statement, tsv_row)
                
                db_rows = connection.execute(statement).fetchall()
                if len(db_rows) == 0:
                    fail(f"Row not found in db: {tsv_row}")
                if len(db_rows) > 1:
                    fail(f"Multiple rows found in db: {tsv_row}")
                row_dict = db_rows[0]._mapping
                for col in data_cols:
                    if isinstance(row_dict[col], Decimal) or isinstance(row_dict[col], float):
                        if math.isclose(row_dict[col], tsv_row[col], rel_tol=1e-7, abs_tol=1e-7):
                            continue
                        else:
                            fail(f"{model} column {col} numerical mismatch: db's {row_dict[col]} != tsv's {tsv_row[col]}")
                    elif isinstance(row_dict[col], str) or isinstance(row_dict[col], int):
                        if row_dict[col] == tsv_row[col]:
                            continue
                        else:
                            fail(f"{model} column {col} string mismatch: db's {row_dict[col]} != tsv's {tsv_row[col]}")
                    elif row_dict[col] is None and tsv_row[col] is None:
                        continue
                    else:
                        print(f"unhandled type {type(row_dict[col])}")
                        quit()
                for check in checks:
                    checkFailMsg = check(row_dict, tsv_row)
                    if checkFailMsg is not None:
                        fail(checkFailMsg)
        
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

    testmodel("transcripts", 
            [transcripts, genes], 
            join_fn=lambda stmt: stmt.join(genes), 
            where_fn=lambda stmt, tsv: stmt.where(transcripts.c.transcript_id == tsv["transcript_id"]),
            data_cols = ['transcript_type', 'tsl'],
            checks = [
                lambda db_row, tsv_row: None if db_row["short_name"] == tsv_row["gene"].upper() else f"Gene_id mismatch: {db_row['short_name']} != {tsv_row['gene'].upper()}"
            ]
            )
    testmodel("variants_transcripts", 
            [variants_transcripts, variants, transcripts], 
            join_fn=lambda stmt: stmt.join(variants).join(transcripts), 
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"], 
                                                    variants.c.assembly == set_var_assembly, 
                                                    transcripts.c.transcript_id == tsv_r["transcript"]), 
            data_cols=['hgvsc']
            )
    
    testmodel("snvs", 
            [snvs, variants], 
            join_fn=lambda stmt: stmt.join(variants), 
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"], 
                                                    variants.c.assembly == set_var_assembly), 
            data_cols=['type', 'length', 'chr', 'pos', 'ref', 'alt', 'cadd_score', 'cadd_intr', 'dbsnp_id', 'dbsnp_url', 'ucsc_url', 'ensembl_url', 'clinvar_url', 'gnomad_url', 'clinvar_vcv', 'splice_ai'],
            )
    
    testmodel("mts",
            [mts, variants],
            join_fn=lambda stmt: stmt.join(variants),
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"],
                                                    variants.c.assembly == set_var_assembly),
            data_cols=['pos', 'ref', 'alt', 'ucsc_url', 'mitomap_url', 'gnomad_url', 'dbsnp_id', 'dbsnp_url', 'clinvar_url', 'clinvar_vcv'],
            )
    
    testmodel("genomic_ibvl_frequencies", 
            [genomic_ibvl_frequencies, variants], 
            join_fn=lambda stmt: stmt.join(variants), 
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"], 
                                                    variants.c.assembly == set_var_assembly), 
            data_cols=['af_tot', 'af_xx', 'af_xy', 'ac_tot', 'ac_xx', 'ac_xy', 'an_tot', 'an_xx', 'an_xy', 'hom_tot', 'hom_xx', 'hom_xy', 'quality']
            )
    
    
    testmodel("genomic_gnomad_frequencies",
            [genomic_gnomad_frequencies, variants],
            join_fn=lambda stmt: stmt.join(variants),
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"],
                                                    variants.c.assembly == set_var_assembly),
            data_cols=['af_tot', 'ac_tot', 'an_tot', 'hom_tot']
            )
    
    testmodel("mt_gnomad_frequencies",
            [mt_gnomad_frequencies, variants],
            join_fn=lambda stmt: stmt.join(variants),
            where_fn=lambda stmt, tsv_r: stmt.where(variants.c.variant_id == tsv_r["variant"],
                                                    variants.c.assembly == set_var_assembly),
            # an	ac_hom	ac_het	af_hom	af_het	max_hl
            data_cols = ['an', 'ac_hom', 'ac_het', 'af_hom', 'af_het', 'max_hl']
            )

    if (did_pass):
        log_output("✅✅✅ \n all tests passed\n✅✅✅\n")
        
    cleanup(None, None)



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
import sys
import os
from datetime import datetime

from dotenv import load_dotenv

from import_utils import *
from model_import_actions import model_import_actions


load_dotenv()

# get command line arguments
rootDir = os.environ.get("PIPELINE_OUTPUT_PATH")
current_dir = os.path.dirname(os.path.realpath(__file__))

chunk_size = int(os.environ.get("CHUNK_SIZE"))
# verbose = os.environ.get("VERBOSE") == "true"

schema = os.environ.get("SCHEMA_NAME")
dry_run = os.environ.get("DRY_RUN") == "true"
update = os.environ.get("UPDATE") == "true"
set_var_assembly = os.environ.get("SET_VAR_ASSEMBLY", None)
set_var_assembly = int(set_var_assembly) if set_var_assembly is not None else None

if (set_var_assembly is None):
    log_error("SET_VAR_ASSEMBLY is not set.")
    quit()

start_at_model = (
    os.environ.get("START_AT_MODEL") if os.environ.get("START_AT_MODEL") != "" else None
)
start_at_file = (
    os.environ.get("START_AT_FILE") if os.environ.get("START_AT_FILE") != "" else None
)


db_row_counts = {"before": {}, "after": {}}

if rootDir == None:
    rootDir = os.path.join(current_dir, "fixtures")

engine = None

metadata = MetaData()

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
                ).join(variants, table.c.variant == variants.c.id).join(transcripts, table.c.transcript == transcripts.c.id).where(variants.c.assembly == set_var_assembly)
                if chromosome is not None:
                    statement = statement.where(variants.c.variant_id.startswith(f"{chromosome}_"))
                
            else:
                cols = [table.c[col] for col in ["id"] + model_action["pk_lookup_col"] ]
                
                if "variant" in model_action["fk_map"]:
                    variants = get_table("variants")
                    cols.append(variants.c["variant_id", "assembly"])
                    statement = select(*cols).join(variants, table.c.variant == variants.c.id).where(variants.c.assembly == set_var_assembly)
                    if chromosome is not None:
                        statement = statement.where(variants.c.variant_id.startswith(f"{chromosome}_"))
                elif modelName == "transcripts":
                    cols.append(table.c["assembly"])
                    statement = select(*cols).where(table.c.assembly == set_var_assembly)
                elif modelName == "variants":
                    cols.append(table.c["assembly"])
                    statement = select(*cols).where(table.c.assembly == set_var_assembly)
                    if chromosome is not None:
                        statement = statement.where(table.c.variant_id.startswith(f"{chromosome}_"))
                else:
                    statement = select(*cols)
                    
            result = connection.execute(statement)
            existing = {
                model_action["map_key_expression"](row): row.id for row in result
            }
            num_in_existing_map = len(existing.keys())
            log_output(f"    caching {str(num_in_existing_map)} {modelName} (chr:{chromosome})")
            return existing
        
    if model == "variants_annotations": # unique variant_transcript per annotation
        depends_on_maps["variants_transcripts"] = make_existing_map("variants_transcripts", chromosome)
        reversed = {v: k for k, v in depends_on_maps["variants_transcripts"].items()}
        tenative_existing_map = make_existing_map(model, chromosome)
        current_model_existing_map = {reversed.get(k): v for k, v in tenative_existing_map.items() if reversed.get(k) is not None }
        
    elif model == "variants_consequences": # non unique variant_transcripts per consequence
        depends_on_maps["variants_transcripts"] = make_existing_map("variants_transcripts", chromosome)
        reversed = {v: k for k, v in depends_on_maps["variants_transcripts"].items()}
        tentative_existing_map = make_existing_map(model, chromosome)
        actual_existing_map = {}
        for vt_id_tuple, id in tentative_existing_map.items():
            vt_id, descriminator = vt_id_tuple
            variant_id, transcript_id = reversed.get(vt_id, (None,None))
            if variant_id is not None and transcript_id is not None and descriminator is not None:
                actual_existing_map[(variant_id, transcript_id, descriminator)] = id
        current_model_existing_map = actual_existing_map #{reversed.get(k): v for k, v in tentative_existing_map.items() if reversed.get(k) is not None }
        
    else:
        referenced_models = action.get("fk_map").values()
        current_model_existing_map = make_existing_map(action["name"], chromosome)
        for m in referenced_models:

            depends_on_maps[m] = make_existing_map(m,chromosome)
    for key in depends_on_maps.keys():
        log_output(f"    depends-on map {key} has: {len(depends_on_maps[key])} items")
    log_output(f"    existing map {model} has: {len(current_model_existing_map)} items")


#    log_output(f"done populating maps for {action['name']} chromosome {chromosome}")
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
    filters = action.get("filters") or {}
    
    table = get_table(model)
    
#    cache_group = file_info('cache_group')
#    print(f"cache group is {file_info['cache_group']}")

    successCount = 0
    failCount = 0
    duplicateCount = 0
    missingRefCount = 0
    insertCount = 0
    updatedCount = 0
    successful_chunks = 0
    fail_chunks = 0
        
    action_types = action.get("tsv_types" ,{})

    df = readTSV(file, file_info, dtype=action_types)
    

    data_insert_list = []
    data_update_list = []
    for _, row in df.iterrows():
        data = row.to_dict()
        
        if (model in ["variants","transcripts"]):
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
            fk = None
            
            if model in ["variants_annotations", "variants_consequences"]:
                depended_map_key = (data["variant"], data["transcript"])
            else:
                
                if isinstance(data[depended_model_col], str):
                    depended_map_key = data[depended_model_col]
            if depended_map_key == "NA":
                
                    data[depended_model_col] = None
            elif depended_map_key is None:
                fk = None
            else:
                fk = depends_on_maps.get(depended_model).get( depended_map_key)
                
            if fk is not None:
                data[depended_model_col] = fk
            else:
                if depended_model_col in action.get("null_ok", []) and data[depended_model_col] is None:
                    pass
                else:
                    log_data_issue(
                        f"Missing {depended_model_col} {depended_map_key} in {depended_model}",
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

    # dispose of df to save ram
    del df
    with engine.connect() as connection:

        def rowOperation(row, updating=False):
            nonlocal successCount, failCount, duplicateCount, insertCount, updatedCount, successful_chunks, fail_chunks
            
            did_succeed = False
            try:
                if updating:
                    
                    id = row["id"]
                    del row["id"]
                    connection.execute(
                        table.update().where(table.c.id == id), row
                    )
                    connection.commit()
                    successCount += 1
                    updatedCount += 1
                    did_succeed = True
                else:
                    connection.execute(table.insert(), row)
                    connection.commit()
                    successCount += 1
                    insertCount += 1
                    did_succeed = True

            except DataError as e:
                log_data_issue(e, model)
                failCount += 1
            except IntegrityError as e:
                msg = str(e)
                if "Duplicate" in msg or "ORA-00001" in msg:
                    duplicateCount += 1
                    successCount += 1
                    log_data_issue(e, model)
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
            nonlocal successCount, failCount, duplicateCount, insertCount, updatedCount, successful_chunks, fail_chunks
            if dry_run:
                successCount += len(chunk)
                return
            if updating:
                
                ids = [row["id"] for row in chunk]
                existing_full_records = connection.execute(select(table).where(table.c.id.in_(ids))).fetchall()
                existing_map = {row.id: row._mapping for row in existing_full_records}
                
                filtered_update_list = []
                
                for row in chunk:
                    existing_row = existing_map.get(row["id"])
                    if existing_row is None:
                        log_data_issue(f"updating, but row with id {row['id']} was removed??", model)
                        failCount += 1
                        continue
                    else:
                        def is_equal(a, b):
                            if a == b:
                                return True
                            if isinstance(a, (float, Decimal)) and isinstance(b, (float, Decimal)):
                                return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)
                            return False
                        
                        if any(not is_equal(existing_row.get(col), row[col]) for col in existing_row.keys()):
                            filtered_update_list.append(row)
                        else:
                            successCount += 1
                            
                for row in filtered_update_list:
                    rowOperation(row, updating)
            else:
                try: 
                    connection.execute(table.insert(), chunk)
                    # commit
                    connection.commit()
                    # chunk worked
                    successful_chunks += 1
                    successCount += len(chunk)
                    insertCount += len(chunk)
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
        "inserted": insertCount,
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




def start(db_engine):

    print("importing", rootDir)
    arrived_at_start_model = False
    arrived_at_start_file = False
    global engine, schema, current_chromosome, last_chromosome
    engine = db_engine

    if isinstance(schema, str) and len(schema) > 0:
        metadata.reflect(bind=engine, schema=schema)
    else:
        metadata.reflect(bind=engine)

    job_dir = get_job_dir()
    os.makedirs(job_dir, exist_ok=True)
    os.chmod(job_dir, 0o777)  # Set read and write permissions for the directory
    setup_loggers(job_dir)

    now = datetime.now()
    counts = {}
    counts["success"] = 0
    counts["fail"] = 0
    counts["missingRef"] = 0
    counts["duplicate"] = 0
    counts["inserted"] = 0
    counts["updated"] = 0
    counts["successful_chunks"] = 0
    counts["fail_chunks"] = 0
    counts["rowcount"] = 0
    
    with engine.connect() as connection:
        table = get_table("severities")
        num_rows = connection.execute(select(func.count()).select_from(table)).scalar()
        if num_rows == 0:
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
        model_counts["inserted"] = 0
        model_counts["updated"] = 0
        model_counts["successful_chunks"] = 0
        model_counts["fail_chunks"] = 0
        model_counts["rowcount"] = 0
        model_directory = os.path.join(rootDir, modelName)

        if (
            isinstance(start_at_model, str)
            and modelName != start_at_model
            and not arrived_at_start_model
        ):
            log_output("Skipping " + modelName + ", until " + start_at_model)
            continue

        if isinstance(start_at_model, str) and modelName == start_at_model:
            arrived_at_start_model = True
        
        log_output(f"*** import {modelName} ***")
        with engine.connect() as connection:
            table = get_table(modelName)
            num_rows = connection.execute(select(func.count()).select_from(table)).scalar()
            db_row_counts["before"][modelName] = num_rows

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
                    "inserted",
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
            
    log_output(f"\n\nfinished importing IBVL. Time Taken: {str(datetime.now() - now)}. was job {job_dir}")
    report_counts(counts)
    log_output("\n\n")
    
    delta = {}
    for beforeafter, counts in db_row_counts.items():
        for modelName, count in counts.items():
            if beforeafter == "before":
                delta[modelName] = count
            else:
                delta[modelName] = count - delta[modelName]
                log_output(f"DB row count {modelName} {count} ( grew by {delta[modelName]})")
    return job_dir
#    cleanup(None, None)

import json
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
copy_maps_from_job = os.environ.get("COPY_MAPS_FROM_JOB")
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

if rootDir == None:
    rootDir = os.path.join(current_dir, "fixtures")
    print("No root directory specified. default ", rootDir)
#    exit()

print("importing", rootDir)

engine = None

metadata = MetaData()

job_dir = ""
maps_load_dir = ""
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
                statement = select(table.c["id"],variants.c["variant_id"], transcripts.c["transcript_id"]).join(variants).join(transcripts)
                if chromosome is not None:
                    statement = statement.where(variants.c.variant_id.startswith(f"{chromosome}_"))
                
            else:
                cols = [table.c[col] for col in ["id"] + model_action["pk_lookup_col"] ]
                
                if "variant" in model_action["fk_map"]:
                    variants = get_table("variants")
                    cols.append(variants.c["variant_id"])
                    statement = select(*cols).join(variants)
                    if chromosome is not None:
                        statement = statement.where(variants.c.variant_id.startswith(f"{chromosome}_"))
                elif modelName == "variants" and chromosome is not None:
                    statement = select(*cols).where(table.c.variant_id.startswith(f"{chromosome}_"))
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
    log_output("cleared the pk maps")


def get_table(model):
    
    if model not in model_import_actions:
        table_name = model
    else:
        table_name = model_import_actions[model]["table"]
        
    if isinstance(schema, str) and len(schema) > 0:
        return Table(table_name, metadata, schema=schema)
    else:
        return Table(table_name, metadata, autoload_with=engine)


#todo add cache group to this and to append_to_map
def inject(model, data, map_key):
    # need to dynamically inject the single obj that was missing from original data
    id = None
    table = get_table(model)
    with engine.connect() as connection:
        try:
            result = connection.execute(table.insert(), data)
            connection.commit()
            id = result.inserted_primary_key[0]
            
            append_to_map(model, map_key, id)
            log_data_issue(f"dynamically added {data}", model)
        except IntegrityError as e:
            log_data_issue("a dynamically injected obj had an integrity error.", model)
            log_data_issue(e, model)
        #                                quit() # LATER: comment this out
        except Exception as e:
            log_data_issue("a dynamically injected obj had an error.", model)
            log_data_issue(e, model)
    #                                quit() # LATER: comment this out?
    return id

def import_file(file, file_info, action):
    global current_chromosome,last_chromosome
    model = action.get("name")
    fk_map = action.get("fk_map")
    pk_lookup_col = action.get("pk_lookup_col")
    filters = action.get("filters") or {}
    
#    cache_group = file_info('cache_group')
#    print(f"cache group is {file_info['cache_group']}")

    missingRefCount = 0
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
    
    df.replace(np.nan, None, inplace=True)

    data_list = []
    for _, row in df.iterrows():
        data = row.to_dict()
        
        if separate_cache_by_chromosome(action):
        
            current_chromosome = data.get("variant",data.get("variant_id")).split("_")[0]
            if current_chromosome != last_chromosome:
                last_chromosome = current_chromosome
                persist_and_unload_maps()
                populate_maps(action, current_chromosome)
        
        

        skip = False
        for col, filter in filters.items():
            data[col] = filter(data[col])
        for fk_col, fk_model in fk_map.items():
            map_key = None
            resolved_pk = None
            debug_row = None
            
            if model in ["variants_annotations", "variants_consequences"]:
                map_key = (data["variant"], data["transcript"])
            else:
                
                if isinstance(data[fk_col], str):
                    map_key = data[fk_col].upper()
            if map_key == "NA":
                    data[fk_col] = None
            else:
                resolved_pk = depends_on_maps.get(fk_model).get( map_key)
                ## resolved PK was not found from maps, so.. if it's a gene, we could dynamically inject
                if (
                    resolved_pk == None
                    and fk_col == "gene"
                    and model == "transcripts"
                    and map_key != None
                ):
                    resolved_pk = inject("genes", {"short_name": map_key}, map_key)
                elif (
                ## resolved PK was not found from maps, so.. if it's a variant, we could dynamically inject
                    resolved_pk == None
                    and fk_col == "variant"
                    and model in ["sv_consequences", "svs", "snvs", "mts"]
                ):

                    if model == "sv_consequences" or model == "svs":
                        var_type = "SV"
                    elif model == "snvs":
                        var_type = "SNV"
                    elif model == "mts":
                        var_type = "MT"
                    resolved_pk = inject(
                        "variants",
                        {"variant_id": map_key, "var_type": var_type},
                        map_key,
                    )
            if resolved_pk is not None:
                data[fk_col] = resolved_pk
            else:
                log_data_issue(
                    "Missing " + fk_col
                    if fk_col is not None
                    else (
                        "None" + " " + map_key
                        if map_key is not None
                        else (
                            "None" + " referenced from " + model
                            if model is not None
                            else "None"
                        )
                    ),
                    model,
                )
                if debug_row is not None:
                    log_data_issue(debug_row, model)
                else:
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

        data_list.append(data)

    # dispose of df to save ram
    del df
    with engine.connect() as connection:
        successCount = 0
        failCount = 0
        duplicateCount = 0
        successful_chunks = 0
        fail_chunks = 0

        for chunk in chunks(data_list, chunk_size):
            if dry_run:
                successCount += len(chunk)
                break
            try:
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
                    did_succeed = False
                    try:
                        connection.execute(table.insert(), row)
                        connection.commit()
                        successCount += 1
                        did_succeed = True

                    except DataError as e:
                        log_data_issue(e, model)
                        failCount += 1
                    #                        quit()
                    except IntegrityError as e:
                        msg = str(e)
                        if "Duplicate" in msg or "ORA-00001" in msg:
                            duplicateCount += 1
                            successCount += 1
                        else:
                            failCount += 1
                            log_data_issue(e, model)
                    #                            quit()
                    except Exception as e:

                        log_data_issue(e, model)
                        failCount += 1
                    ####### added in v2
                    if not did_succeed:
                        connection.rollback()

    return {
        "success": successCount,
        "fail": failCount,
        "missingRef": missingRefCount,
        "duplicate": duplicateCount,
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
    global job_dir, maps_load_dir, engine, schema, current_chromosome, last_chromosome
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

    if copy_maps_from_job is not None and copy_maps_from_job != "":
        maps_load_dir = os.path.join(jobs_dir, copy_maps_from_job)
    else:
        maps_load_dir = job_dir
    print("using job dir " + maps_load_dir)

    now = datetime.now()
    counts = {}
    counts["success"] = 0
    counts["fail"] = 0
    counts["missingRef"] = 0
    counts["duplicate"] = 0
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
        {"name":"severities", "fk_map":{}, "pk_lookup_col":None},
    )
    log_output("done importing severities")


    for modelName, action_info in model_import_actions.items():
        model_counts = {}
        model_counts["success"] = 0
        model_counts["fail"] = 0
        model_counts["missingRef"] = 0
        model_counts["duplicate"] = 0
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
    print("finished importing IBVL. Time Taken: " + str(datetime.now() - now))
    report_counts(counts)
    cleanup(None, None)

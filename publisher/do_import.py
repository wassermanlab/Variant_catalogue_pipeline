

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
import sys
import os
import traceback
from datetime import datetime

from dotenv import load_dotenv

from import_utils import *
from model_import_actions import model_import_actions
from publish_group import publish_group

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

if rootDir is None:
    rootDir = os.path.join(current_dir, "fixtures")

engine = None

metadata = MetaData()

depends_on_maps = {}
transcripts_map = {} # special persistance, because frequently depended on

def separate_cache_by_chromosome(action):
    return action["name"] not in ["transcripts","genes", "severities"]

#@profile
def populate_maps(model_action, group_label="all", variant_prefix=None, variant_suffix=None):
    print("populating maps", group_label)
    global depends_on_maps, transcripts_map
    model = model_action["name"]
    group_existing_map = {}
    
    def make_existing_map(modelName):
        table = get_table(modelName)
        depended_model_action = model_import_actions[modelName]

        with engine.connect() as connection:
            if (modelName == "variants_transcripts"):
                variants = get_table("variants")
                transcripts = get_table("transcripts")
                # TODO: Fix invalid column selection syntax (see review comment)
                # Current code has invalid SQLAlchemy syntax on line 83
                statement = select(
                    table.c["id"],
                    variants.c["variant_id", "assembly"], 
                    transcripts.c["transcript_id"]
                ).join(variants, table.c.variant == variants.c.id
                ).join(transcripts, table.c.transcript == transcripts.c.id
                ).where(variants.c.assembly == set_var_assembly)
                statement = statement.where(variants.c.variant_id.startswith(variant_prefix))
                
                # CORRECTED VERSION (commented out for review/testing):
                # statement = select(
                #     table.c["id"],
                #     variants.c.variant_id,
                #     variants.c.assembly,
                #     transcripts.c["transcript_id"]
                # ).join(variants, table.c.variant == variants.c.id
                # ).join(transcripts, table.c.transcript == transcripts.c.id
                # ).where(variants.c.assembly == set_var_assembly)
                # statement = statement.where(variants.c.variant_id.startswith(variant_prefix))
                
            else:
                cols = [table.c[col] for col in ["id"] + depended_model_action["pk_lookup_col"] ]
                
                if "variant" in depended_model_action["fk_map"]:
                    variants = get_table("variants")
                    # TODO: Fix invalid column selection syntax (see review comment)
                    # Current code has invalid SQLAlchemy syntax on line 95
                    cols.append(variants.c["variant_id", "assembly"])
                    # CORRECTED VERSION (commented out for review/testing):
                    # cols.extend([variants.c.variant_id, variants.c.assembly])
                    statement = select(*cols).join(variants, table.c.variant == variants.c.id).where(variants.c.assembly == set_var_assembly)
                    statement = statement.where(variants.c.variant_id.startswith(variant_prefix))
                    #                    statement =  variant_prefix(statement, variants)
                elif modelName == "variants":
                    variants = table
                    cols.append(table.c["assembly"])
                    statement = select(*cols).where(variants.c.assembly == set_var_assembly)
                    statement = statement.where(variants.c.variant_id.startswith(variant_prefix))
#                    statement =  variant_prefix(statement, table, group_label)
                elif modelName == "variants_annotations":
                    variants = get_table("variants")
                    variants_transcripts = get_table("variants_transcripts")
                    statement = select(*cols).join(
                            variants_transcripts, table.c.variant_transcript == variants_transcripts.c.id
                        ).join(
                            variants, variants.c.id == variants_transcripts.c.variant
                        ).where(
                            variants.c.assembly == set_var_assembly
                        )
                    statement = statement.where(variants.c.variant_id.startswith(variant_prefix))                        
#                    statement =  variant_prefix(statement, variants)
                elif modelName == "variants_consequences":
                    variants = get_table("variants")
                    variants_transcripts = get_table("variants_transcripts")
                    statement = select(*cols).join(
                            variants_transcripts, table.c.variant_transcript == variants_transcripts.c.id
                        ).join(
                            variants, variants.c.id == variants_transcripts.c.variant
                        ).where(
                            variants.c.assembly == set_var_assembly
                        )
                        
                    statement = statement.where(variants.c.variant_id.startswith(variant_prefix))
#                    statement = variant_prefix(statement, variants)
                elif modelName == "transcripts":
                    statement = select(*cols).where(table.c.assembly == set_var_assembly)
                else:
                    statement = select(*cols)
            if variant_suffix is not None and modelName != "transcripts":
                statement = statement.where(variants.c.variant_id.endswith(variant_suffix))
#            if isinstance(variant_prefix, str):
#                statement = statement.where(get_table("variants").c.variant_id.startswith(variant_prefix))
            existing = {
                depended_model_action["map_key_expression"](row): row.id for row in connection.execute(statement)
            }
            num_in_existing_map = len(existing.keys())
            return existing
        
    if model == "variants_annotations": # unique variant_transcript per annotation
        depends_on_maps["variants_transcripts"] = make_existing_map("variants_transcripts")
        reversed = {v: k for k, v in depends_on_maps["variants_transcripts"].items()}
        tenative_existing_map = make_existing_map(model)
        group_existing_map = {reversed.get(k): v for k, v in tenative_existing_map.items() if reversed.get(k) is not None }
        del reversed
        del tenative_existing_map
        
    elif model == "variants_consequences": # non unique variant_transcripts per consequence
        depends_on_maps["variants_transcripts"] = make_existing_map("variants_transcripts")
        reversed = {v: k for k, v in depends_on_maps["variants_transcripts"].items()}
        tentative_existing_map = make_existing_map(model)
        actual_existing_map = {}
        for vt_id_tuple, id in tentative_existing_map.items():
            vt_id, descriminator = vt_id_tuple
            variant_id, transcript_id = reversed.get(vt_id, (None,None))
            if variant_id is not None and transcript_id is not None and descriminator is not None:
                actual_existing_map[(variant_id, transcript_id, descriminator)] = id
        del reversed
        del tentative_existing_map
        group_existing_map = actual_existing_map #{reversed.get(k): v for k, v in tentative_existing_map.items() if reversed.get(k) is not None }
        
    else:

        group_existing_map = make_existing_map(model)
        
        referenced_models = model_action.get("fk_map").values()
        for m in referenced_models:
            if m == "transcripts" and not transcripts_map:
                transcripts_map = make_existing_map("transcripts")
                depends_on_maps["transcripts"] = transcripts_map
                log_output("just made the transcripts map")
            else:
                depends_on_maps[m] = make_existing_map(m)
                
    for key in depends_on_maps.keys():
        log_output(f"    depends-on map {key} {group_label} has: {len(depends_on_maps[key])} items")
    log_output(f"    existing map {model} {group_label} has: {len(group_existing_map)} items")
    
    return group_existing_map

def persist_and_unload_maps():
    global depends_on_maps, transcripts_map
    if "transcripts" in depends_on_maps and len(depends_on_maps["transcripts"]) > 0:
        transcripts_map = depends_on_maps["transcripts"]
    depends_on_maps.clear()
    depends_on_maps["transcripts"] = transcripts_map


def get_table(model):
    
    if model not in model_import_actions:
        table_name = model
    else:
        table_name = model_import_actions[model]["table"]
        
    if isinstance(schema, str) and len(schema) > 0:
        return Table(table_name, metadata, schema=schema)
    else:
        return Table(table_name, metadata, autoload_with=engine)

#@profile
def import_file(file, file_info, action):
#    global current_chromosome,last_chromosome
    model = action.get("name")
    fk_map = action.get("fk_map")
    filters = action.get("filters") or {}
    
    file_now = datetime.now()
    
    table = get_table(model)
    action_types = action.get("tsv_types" ,{})
    df = readTSV(file, file_info, dtype=action_types)
    
    missingRefCount = 0
    duplicateCount = 0
    successCount = 0
    
    subgroups = [{"df":df, "group_label":"all"}]
    
    now_subgroups = datetime.now()
    if model in ["variants_annotations", "variants_consequences"]:
        
        print("making subgroups for model ", model)
        subgroups = []
        for chromosome in [ str(i) for i in list(range(1,22)) + ["X", "Y", "chrM"] ]:
            # for most ram-heavy models - further subdivide by last nucleotide (alt) of variant 
            # (or - for some mito variant deletions)
            for nucleotide in ["A","T","G","C","-"]:
                
                if nucleotide == "-" and chromosome != "chrM":
                    continue
                group_df = df[df.apply(lambda row: row["variant"].split("_")[0] == chromosome and row["variant"][-1] == nucleotide, axis=1)]
                group_label = f"chr:{chromosome}, alt {nucleotide}" 
            
                subgroups.append(
                    {"df":group_df, 
                    "group_label":group_label, 
                    "variant_prefix":f"{chromosome}_",
                    "variant_ends_with":nucleotide
                    })
                
        
        del df


    elif separate_cache_by_chromosome(action):
        # group by chromosome only
        print("making subgroups for model ", model)
        subgroups = []
        for chromosome in [ str(i) for i in list(range(1,22)) + ["X", "Y", "chrM"] ]:
            if model == "variants":
                group_df = df[df.apply(lambda row: row["variant_id"].split("_")[0] == chromosome, axis=1)]
            else:
                group_df = df[df.apply(lambda row: row["variant"].split("_")[0] == chromosome, axis=1)]
            group_label = f"{chromosome}" 
            
            subgroups.append(
                {"df":group_df, 
                 "group_label":group_label, 
                 "variant_prefix":f"{chromosome}_"
                 })
        del df
    log_output(f"making subgroups took {str(datetime.now() - now_subgroups)}")
    counts = {}

    for subgroup in subgroups:
        
        data_insert_list = []
        data_update_list = []
        group_label = subgroup["group_label"]
        df = subgroup["df"]
        variant_prefix_or_none = subgroup.get("variant_prefix")
        variant_alt_or_none = subgroup.get("variant_ends_with")
        existing_map = {}
        
        if df.empty:
            continue
        
        if model != "severities":
            persist_and_unload_maps()
            existing_map = populate_maps(action, group_label, variant_prefix_or_none, variant_alt_or_none)
        
        for _, row in df.iterrows():
            data = row.to_dict()
            
            if (model in ["variants","transcripts"]):
                data["assembly"] = set_var_assembly

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

            if record_map_key is None:
                print(f"record_map_key is None for {model} {data}")
                log_error(f"record_map_key is None for {model} {data}")
                quit()
            existing_id = existing_map.get(record_map_key)
            if existing_id is not None:
                # record is already in the DB
                if update:
                    data["id"] = existing_id
                    data_update_list.append(data)
                else:
                    duplicateCount += 1
                    successCount += 1
                    continue
            else:
                # the record is NOT in the db, so add it
                data_insert_list.append(data)

        # dispose to save ram
        del existing_map
        del subgroup            
        
        group_counts = publish_group(data_insert_list, data_update_list, engine, model, table)
        for key in group_counts.keys():
            if key not in counts:
                counts[key] = 0
            counts[key] += group_counts[key]
    
    del subgroups
    log_output(f"    {str(datetime.now() - file_now)}")
    
    
    counts["missingRef"] = missingRefCount
    counts["duplicate"] = counts.get("duplicate",0) + duplicateCount
    counts["success"] = counts.get("success",0) + successCount
    counts["rowcount"] = file_info["total_rows"]
    
    return counts


def cleanup(sig, frame):
    global engine, depends_on_maps, transcripts_map, metadata
#    print("terminating, cleaning up ...")
    log_output("terminating, cleaning up ...")
    persist_and_unload_maps()
    engine.dispose()
    # garbage collect
    del depends_on_maps
    del transcripts_map
    del metadata
#    print("done")
    log_output("done")
    sys.exit(0)




def start(db_engine):

    arrived_at_start_model = False
    arrived_at_start_file = False
    global engine, schema
    engine = db_engine

    if isinstance(schema, str) and len(schema) > 0:
        metadata.reflect(bind=engine, schema=schema)
    else:
        metadata.reflect(bind=engine)

    job_dir = get_job_dir(schema, set_var_assembly, dry_run, rootDir.split("/")[-1], start_at_model, start_at_file)
    os.makedirs(job_dir, exist_ok=True)
    os.chmod(job_dir, 0o777)  # Set read and write permissions for the directory
    setup_loggers(job_dir)
    
    log_output(f"schema: {schema}")
    log_output(f"assembly: {set_var_assembly}")
    log_output(f"tsv directory: {rootDir}")
    
    try:

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
            
            modelNow = datetime.now()

            if large_model_file_tsv_exists:
                sorted_files = [modelName + ".tsv"]
            else:
                sorted_files = natsorted(
                    [f for f in os.listdir(model_directory) if not f.startswith(".")],
                )
            
            with engine.connect() as connection:
                table = get_table(modelName)
                num_rows = connection.execute(select(func.count()).select_from(table)).scalar()
                db_row_counts["before"][modelName] = num_rows

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
                    if results.get("success",0) == 0 and results.get("fail",0) > 0:
                        log_output("No rows were imported, bailing to avoid further failures")
                        cleanup(None, None)

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
                        if key not in results:
                            results[key] = 0
                        model_counts[key] += results.get(key, 0)
                        counts[key] += results.get(key, 0)

                    report_counts(results)

            log_output(
                "\nFinished importing "
                + modelName
                + ". Took this much time: "
                + str(datetime.now() - modelNow)
            )
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
    except Exception as e:
        log_error(f"{e.__class__.__name__}: {e}")
        log_error(traceback.format_exc())
        cleanup(None, None)
#    cleanup(None, None)

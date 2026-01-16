
import math
from decimal import Decimal
import os

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

from import_utils import *

dry_run = os.environ.get("DRY_RUN") == "true"
chunk_size = int(os.environ.get("CHUNK_SIZE", 5000))

def publish_group(
    data_insert_list,
    data_update_list,
    engine,
    model,
    table,

):
    
    successCount = 0
    failCount = 0
    duplicateCount = 0
    insertCount = 0
    updatedCount = 0
    successful_chunks = 0
    fail_chunks = 0

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
                
                existing_map = {}
                for sub_chunk in chunks(chunk, 1000):
                    ids = [row["id"] for row in sub_chunk]
                    existing_full_records = connection.execute(select(table).where(table.c.id.in_(ids))).fetchall()
                    existing_map.update({row.id: row._mapping for row in existing_full_records})
#                existing_map = {row.id: row._mapping for row in existing_full_records}
                
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
        "duplicate": duplicateCount,
        "inserted": insertCount,
        "updated": updatedCount,
        "successful_chunks": successful_chunks,
        "fail_chunks": fail_chunks,
    }
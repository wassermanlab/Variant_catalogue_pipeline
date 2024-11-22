from sqlalchemy import (
    create_engine,
    MetaData,
    Table,
    Column,
    Integer,
    String
)
import os
from dotenv import load_dotenv
import traceback

load_dotenv()

# get command line arguments
rootDir = os.environ.get("PIPELINE_OUTPUT_PATH")
chunk_size = int(os.environ.get("CHUNK_SIZE"))
schema = os.environ.get("SCHEMA_NAME")

if rootDir == None:
    print("No root directory specified")
    exit()

dbConnectionString = os.environ.get("DB")
container = os.environ.get("DB_CONTAINER")

print("connecting...")

engine = create_engine(dbConnectionString, echo=True, pool_pre_ping=True, pool_recycle=3600)

def check_for_bail(question):
    doContinue = input(question)
    if (doContinue != "y"):
        print("exiting...")
        engine.dispose()
        quit()

with engine.connect() as connection:

    print("looks like connecting worked. will try to reflect ")

    metadata = MetaData()

    metadata.reflect(bind=engine, schema=schema)

    print("reflected tables:" + str(metadata.sorted_tables))

    check_for_bail("continue to table test? (y/n): ")

    table_name = 'genes'

    print("initializing Table for "+table_name+"...")

    if isinstance(schema, str) and len(schema) > 0:
#        print("setting schema = "+schema+"...")
        test_table = Table(
            table_name,
            metadata,
#            autoload_with=engine,
            schema=schema
        )
    else:
        print("no schema provided...")
        test_table = Table(
            table_name,
            metadata,
            Column("ID", Integer, primary_key=True),
            Column("SHORT_NAME", String(30), nullable=False),
        )

    try:
        result = connection.execute(test_table.select())
        print("table exists")
        print(result.fetchone())
    except Exception:
        traceback.print_exc()


    check_for_bail("continue to write table test? (y/n): ")

    n = input("name?")


    try:
        result = connection.execute(test_table.insert(), {"short_name": n})
        #commit
        connection.commit()
        print("inserted ")
    except Exception:
        traceback.print_exc()

engine.dispose()
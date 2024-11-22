from do_import import start
from dotenv import load_dotenv
import os
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import (
    create_engine,
    event
)

load_dotenv()

dbConnectionString = os.environ.get("DB")
verbose = os.environ.get("VERBOSE") == "true" or os.environ.get("VERBOSE") == "True"
isDevelopment = os.environ.get("ENVIRONMENT") != "production"

if isDevelopment:
    # create the database if it doesn't exist
    engine = create_engine(dbConnectionString, echo=False)
    if database_exists(engine.url):
        #assume already has structure
        pass
    else:
        create_database(dbConnectionString)
        import tables
    engine.dispose()


print("connecting...")
engine = create_engine(dbConnectionString, echo=False, pool_pre_ping=True, pool_recycle=3600)

start(engine)   

exit()


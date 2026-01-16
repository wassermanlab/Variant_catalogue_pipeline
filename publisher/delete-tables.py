from do_import import start
from dotenv import load_dotenv
import os
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import (
    create_engine,
    MetaData,
    Table,
    text,
    event
)
from sqlalchemy.orm import sessionmaker

load_dotenv()

dbConnectionString = os.environ.get("DB")
verbose = os.environ.get("VERBOSE") == "true" or os.environ.get("VERBOSE") == "True"

print("connecting...")
engine = create_engine(dbConnectionString, echo=False, pool_pre_ping=True, pool_recycle=3600)
print("connected")

metadata = MetaData()

# Reflect the database
metadata.reflect(bind=engine)

# Create a session
Session = sessionmaker(bind=engine)
session = Session()

# Disable foreign key checks
session.execute(text("SET FOREIGN_KEY_CHECKS = 0;"))
session.commit()

# Truncate each table
for table in metadata.sorted_tables:
    session.execute(text(f"TRUNCATE TABLE {table.name};"))
    #reset autoincrement
    session.execute(text(f"ALTER TABLE {table.name} AUTO_INCREMENT = 1;"))
    session.commit()

# Re-enable foreign key checks
session.execute(text("SET FOREIGN_KEY_CHECKS = 1;"))
session.commit()

# Close the session
session.close()

print("all tables deleted (truncated)")
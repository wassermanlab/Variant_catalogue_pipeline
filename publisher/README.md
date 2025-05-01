
## Publisher: upload the finished pipeline data into the portal

How to run:

  0) (re)create the database (eg, for a mySQL db: `mysql -u root -e "DROP DATABASE IF EXISTS ibvltest; CREATE DATABASE ibvltest;"`
  1) copy the `.env-sample` file to `.env` and set values appropriately
  2) (optional - for development purposes) run `python tables.py` to create the tables (database should be empty before this)
  3) `python publish.py` will kick off the migration

## Notes:

The script creates a directory called "jobs", and additional subdirectories every time it is run. Each of these job folders has log, error output and an error / warning list for each model.

You can run delete-tables.py in between runs while developing to flush out all the data.

model_import_actions.py defines the list of tables and lambda functions related to routine operations done during the import. Most of the custom model functionality is here, but are still some "magic" operations happening outside these (for example: the If statements in https://github.com/wassermanlab/Variant_catalogue_pipeline/blob/publisher-dev/publisher/do_import.py#L67

test-db.py can be used to verify connectivity with an Oracle DB

## Import environment vars
  - `PIPELINE_OUTPUT_PATH` - the full path to the directory containing pipeline output files ( optional - defaults to test/fixures )
  - `SCHEMA_NAME` - for an Oracle destination db, the schema name (database name) goes here. ( for non-Oracle, probably just use database name )
  - `START_AT_MODEL` - to pick up after a previous migration run left off, you can enter the model name here, and the script will skip to that model (it runs in the order of keys as defined in the `model_import_actions` map)
  - (`START_AT_FILE`) - for convenience, you can also skip to a particular file in the first model dir imported, using natural sorting. Be very careful if using this in production as it will lead to false duplicates unless the primary key for new row insertions is corrected.

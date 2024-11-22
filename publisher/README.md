
## Publisher: upload the finished pipeline data into the portal

How to run an import:
  1) copy the `import/.env-sample` file to `import/.env` and set values appropriately
  2) (optional) if you need to, run `python tables.py` to create the tables (database should be empty before this)
  3) `python orchestrate.py` will kick off the migration

The script creates a directory called "jobs", and a directory inside that called "1" the first time, "2" the second time, eg. 

Each of these job folders has working data for the migration and two output logs (one for errors, one for progress). The working data is just (for each model) a file with the latest primary key, and a reverse lookup map for entity id (eg gene or variant or transcript id) to primary key.

### Import environment vars
  - `PIPELINE_OUTPUT_PATH` - the full path to the directory containing pipeline output files
  - `JOBS_PATH` - relative path (from execution path) to folder to contain jobs.
  - `COPY_MAPS_FROM_JOB` - The script maintains maps in order to resolve primary keys, they are persisted as json to the job directory, named using an incrementing number. If a job fails and you want to use the maps from a previous run, enter the run's job folder number as the value of this environment variable
  - `SCHEMA_NAME` - for an Oracle destination db, the schema name goes here.
  - `START_AT_MODEL` - to pick up after a previous migration run left off, you can enter the model name here, and the script will skip to that model (it runs in the order of keys as defined in the `model_import_actions` map)
  - (`START_AT_FILE`) - for convenience, you can also skip to a particular file in the first model dir imported, using natural sorting. Be very careful if using this in production as it will lead to false duplicates unless the primary key for new row insertions is corrected.
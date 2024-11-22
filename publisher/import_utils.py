
import logging
from datetime import datetime
import pandas as pd

from dotenv import load_dotenv
import os

load_dotenv()

chunk_size = int(os.environ.get("CHUNK_SIZE"))
verbose = os.environ.get("VERBOSE") == "true" or os.environ.get("VERBOSE") == "True"
print('verbose is', verbose)

def inspectTSV(file):
    total_rows = 0
    separator = "\t"
    small_read = pd.read_csv(file, sep="\t", nrows=4, header=None)

    num_columns = small_read.shape[1]
    if num_columns == 1 and "gene" not in file:
        separator = " "
        small_read = pd.read_csv(file, sep=" ", nrows=4, header=None)
        if small_read.shape[1] == 1:
            print("Could not determine separator for " + file)
            quit()

    columns = [col.lower() for col in small_read.values.tolist()[0]]

    for chunk in pd.read_csv(file, sep="\t", chunksize=chunk_size):
        total_rows += len(chunk)

    return {
        "total_rows": total_rows,
        "num_columns": num_columns,
        "columns": columns,
        "types": {},
        "separator": separator,
    }

def readTSV(file, info, dtype={}):
#    df = pd.read_csv(file, sep=info["separator"], dtype=dtype, na_values=["NA"], keep_default_na=False)
    df = pd.read_csv(file, sep=info["separator"])
    df.rename(columns={"All_info$variant": "variant"}, inplace=True)
    df.columns = [col.lower() for col in df.columns]
    return df

def setup_loggers(job_dir):
    global data_issue_logger, output_logger  # Add global keyword
    data_issue_logger = logging.getLogger("data_issues")
    data_issue_logger.setLevel(logging.WARNING)
    data_issue_logger.propagate = False
    output_logger = logging.getLogger("output")
    output_logger.setLevel(logging.INFO)
    output_logger.propagate = False

    data_issue_handler = logging.FileHandler(os.path.join(job_dir,"data_issues.log"))
    data_issue_handler.setLevel(logging.WARNING)
    data_issue_logger.addHandler(data_issue_handler)

    output_logger_handler = logging.FileHandler(os.path.join("./",job_dir,"output.log"))
    output_logger_handler.setLevel(logging.INFO)
    output_logger.addHandler(output_logger_handler)

    output_logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    print("logging to", os.path.join(job_dir,"data_issues.log"), os.path.join(job_dir,"output.log"))

def log_data_issue(s):
    data_issue_logger.warning(s)
    if (verbose):
        print(s)
def log_output(s):
    output_logger.info(s)
    if (verbose):
        print(s)


def report_counts(counts):
    percent_success = "N/A"
    if counts["success"] + counts["fail"] > 0:
        percent_success = (
            100
            * counts["success"]
            / (counts["success"] + counts["fail"])
        )
    log_output(
        str(percent_success)
        + "% success. Count: "
        + str(counts["success"])
        + ", Total fail: "
        + str(counts["fail"])
        + ". Total missing refs: "
        + str(counts["missingRef"])
        + ". Total duplicates: "
        + str(counts["duplicate"])
        + ". Total successful chunks: "
        + str(counts["successful_chunks"])
        + ". Total fail chunks: "
        + str(counts["fail_chunks"])
    )


def chunks(l, n):
    """Yield successive n-sized chunks from list l."""
    for i in range(0, len(l), n):
        yield l[i : i + n]
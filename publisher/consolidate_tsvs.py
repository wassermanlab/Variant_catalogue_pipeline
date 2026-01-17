# consolidates the tsvs in a pipeline output folder into one single tsv per model (rather than per chromosome)
# will truncate part of the filename after the last _ to create the consolidated tsv name
# note: columns might have a different ordering, and rounding will occur

from do_import import *
from import_utils import *
from pandas import DataFrame

from model_import_actions import model_import_actions

from dotenv import load_dotenv
import os
load_dotenv()

#get folder path from arg -p
path = os.environ.get('PIPELINE_OUTPUT_PATH')#"~/sg_test_data/6chrs/transcripts/"

output_path = path+"_consolidated"

#iterate through folders in path
import os
for folder in os.listdir(path):
    
    folder_path = os.path.join(path, folder)
    if os.path.isdir(folder_path):
        print("Processing folder:", folder, folder_path)
        tsv_frame = None
        tsv_info = None
        for file in os.listdir(folder_path):
            if file.endswith(".tsv"):
                file_path = os.path.join(folder_path, file)
                print("  Found TSV file:", file_path)
                if tsv_frame is None:
                    tsv_info = inspectTSV(file_path)
                    tsv_frame = pd.DataFrame(columns=tsv_info['columns'])
                    
                new_tsv_frame = pd.read_csv(file_path, sep=tsv_info["separator"])
            tsv_frame = pd.concat([tsv_frame, new_tsv_frame], ignore_index=True)
        os.makedirs(output_path, exist_ok=True)
        os.makedirs(os.path.join(output_path, folder), exist_ok=True)
        #copy to output_path with folder name prefixed
        file_parts = file.split("_")
        if len(file_parts) < 2:
            new_name = file
        else:
            new_name = "_".join(file_parts[:-1]) + ".tsv"
        if folder == "variants":
            print("variants folder detected, adding filter column")
            tsv_frame["filter"] = ""
            
            
        
        dest_file = os.path.join(output_path, folder,new_name)
        print("    Copying to:", dest_file)
        #write tsv_frame to dest_file as TSV
        tsv_frame.to_csv(dest_file, sep=tsv_info["separator"], index=False)
        tsv_frame = None
        tsv_info = None
exit()

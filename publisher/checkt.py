from do_import import *
from import_utils import *

from model_import_actions import model_import_actions


transcripts = get_table("transcripts")
t1_file = os.environ.get('T1')#"~/sg_test_data/6chrs/transcripts/transcripts_10.tsv"
t2_file = os.environ.get('T2')#"~/sg_test_data/2024hg38_test_v7/transcripts/transcripts_chr10.tsv"

t1_info = inspectTSV(t1_file)
t1 = readTSV(t1_file, t1_info)

t2_info = inspectTSV(t2_file)
t2 = readTSV(t2_file, t2_info)


print("t1 info", t1_info)
print("t2 info", t2_info)

transcript_map = {}
transcript_map_none = {}

for i, row in t1.iterrows():
    transcript_map[row.transcript_id] = {'t1_gene': row.gene}

for i, row in t2.iterrows():
    if row.transcript_id in transcript_map:
        transcript_map[row.transcript_id]['t2_gene'] = row.gene
        if transcript_map[row.transcript_id]['t1_gene'] == transcript_map[row.transcript_id]['t2_gene']:
            del transcript_map[row.transcript_id]
        
for k,v in transcript_map.items():
    print(k, v)
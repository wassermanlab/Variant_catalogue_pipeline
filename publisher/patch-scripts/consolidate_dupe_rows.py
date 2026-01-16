import pandas as pd
import argparse
import os
from io import StringIO

pd.set_option('future.no_silent_downcasting', True)

freq_types =         {
            'variant': str,
            'af_tot': float, 
            'af_xx': float, 
            'af_xy': float,
            'ac_tot': int,
            'ac_xx': int,
            'ac_xy': int,
            'an_tot': int,
            'an_xx': int,
            'an_xy': int,
            'hom_tot': int,
            'hom_xx': int,
            'hom_xy': int,
            'quality': float
        }
num_processed = 0
resolved_with_quality = 0
same_ac_and_q = 0

def process_tsv(input_file,output_file):
    global resolved_with_quality, num_processed
    print("processing file", input_file)
    df = pd.read_csv(input_file, sep='\t', dtype=freq_types, na_values=['.'])
#    df = pd.read_csv(input_file, sep='\t')
#    df.replace('.', 0.00000, inplace=True)
    df = df.astype(freq_types)
    
    duplicates = df[df.duplicated(subset=["variant"], keep=False)]
    df_duplicates = df[df["variant"].isin(duplicates["variant"])]
    
    columns_to_sum = [col for col in df_duplicates.columns if col not in ["quality", "an_tot","an_xx","an_xy"] and pd.api.types.is_numeric_dtype(df_duplicates[col])]
    grouped = df_duplicates.groupby("variant")
    
    new_rows = []
    
    for group_name, group in grouped:
#        print("group_name", group_name)
#        print("group quality", group[["quality"]])
#        print("group numeric cols", group[columns_to_sum])
#        print("idx with highest ac_tot", group["ac_tot"].idxmax())
#        print("group quality with highest ac_tot", group.loc[group["ac_tot"].idxmax()]["quality"])
#        print("summed", group[columns_to_sum].sum())

        def resolve_group(g):
            global resolved_with_quality
            
            n = len(g)
            
            if g["ac_tot"].nunique() == n:   
                # clear case: all unique, pick max ac.
                new_row = g.loc[g["ac_tot"].idxmax()].copy()
            elif g["ac_tot"].nunique() == 1:
                resolved_with_quality += 1
                print(group_name)
                new_row = g.loc[g["quality"].idxmax()].copy()
            elif g["ac_tot"].nunique() < n:
                # remove min ac_tot
                group_without_min = g[g["ac_tot"] != g["ac_tot"].min()]
                new_row = resolve_group(group_without_min)
            else:
                print("shouldn't be here")
                raise ValueError("shouldn't be here")
            return new_row
            

            
        new_rows.append(resolve_group(group))
        num_processed += 1
#        try: -- uncomment if you want to see what rows are failing to enter
#            test_df = pd.DataFrame([new_row]).astype(freq_types)
#        except:
#            print("new row fails to enter", new_row)
#            raise
 
    if (len(new_rows) > 0):
        
        try:
            new_df = pd.DataFrame(new_rows).astype(freq_types) 
        except:
            print("unable to make dataframe from", new_rows)
            raise
#        print("new new", new_df)
    
        new_df.to_csv(output_file, sep='\t', index=False)
        print("done:", output_file)
    else:
        print("no duplicates found in", input_file)
def main():
    parser = argparse.ArgumentParser(description='Produce consolidated duplicated rows from TSV file.')
    parser.add_argument('--input_dir', type=str, help='Path to the dir containing input TSV files', default='.')
    parser.add_argument('--output_dir', type=str, help='Path to the output dir', default='./out')

    args = parser.parse_args()
    files_in_dir = os.listdir(args.input_dir)
    for file in files_in_dir:
        if file.endswith(".tsv") and not file.startswith("0_indelpatch_"):
            input_file = os.path.join(args.input_dir, file)
            output_file = os.path.join(args.output_dir, f"0_indelpatch_{file}")
            process_tsv(input_file, output_file)
    print("done creating patch TSVs. Number of rows:", num_processed)
    print("this is how many dupe indels had same allele count:", resolved_with_quality)


if __name__ == '__main__':
    main()

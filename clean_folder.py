from pathlib import Path
import shutil
import argparse

def parse1():
    parser = argparse.ArgumentParser(description="clean_folder")
   
    parser.add_argument('--output', type=str, help="Name of output file")
    
    args = parser.parse_args()
    return args

parser = parse1()
output =parser.output
#print(output)



src_dir = Path(".")
dst_dir_extra = Path("DAFT_extras")
dst_dir_results = Path("DAFT_results")

dst_dir_extra.mkdir(exist_ok=True)
dst_dir_results.mkdir(exist_ok=True)

prefixes = ("introgression", "rev", "result" , "djiNNI", "list_")
#prefixes = ("introgression", "rev", "result" , "djiNNI")


for csv in src_dir.glob("*.csv"):
    if csv.name.startswith(prefixes):
        shutil.move(csv, dst_dir_extra / csv.name)


shutil.move(
    f"DAFT_Test_{output}.txt",
    Path(dst_dir_results, f"DAFT_Test_{output}.txt")
)



files = [
    ('DAFT_Direction_'+output+'.txt', 'DAFT_Direction_'+output+'.txt'),
    ('DAFT_Test_'+output+'.xlsx', 'DAFT_Test_'+output+'.xlsx'),
    ('important_'+output+'.csv', 'important_'+output+'.csv'),
    ('branch_map.csv', 'branch_map.csv'),
]

for src, dst in files:
    try:
        shutil.move(src, dst_dir_results / dst)
    except FileNotFoundError:
        pass
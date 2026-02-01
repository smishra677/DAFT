from pathlib import Path
import shutil

src_dir = Path(".")
dst_dir_extra = Path("DAFT_extras")
dst_dir_results = Path("DAFT_results")

dst_dir_extra.mkdir(exist_ok=True)
dst_dir_results.mkdir(exist_ok=True)

prefixes = ("introgression", "rev", "result")


for csv in src_dir.glob("*.csv"):
    if csv.name.startswith(prefixes):
        shutil.move(csv, dst_dir_extra / csv.name)


shutil.move('DAFT_Significance.txt', dst_dir_results / 'DAFT_Significance.txt')
shutil.move('DAFT_Direction.txt', dst_dir_results / 'DAFT_Direction.txt')

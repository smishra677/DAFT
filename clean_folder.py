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

files = [
    ('DAFT_Direction.txt', 'DAFT_Direction.txt'),
    ('DAFT_Significance.xlsx', 'DAFT_Significance.xlsx'),
]

for src, dst in files:
    try:
        shutil.move(src, dst_dir_results / dst)
    except FileNotFoundError:
        pass
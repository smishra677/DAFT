import pandas as pd
import subprocess
import os




CHUNK_SIZE = 2000
OUTPUT_DIR = "./"
sim_dict = {
    "P": "DCMD_1",
}
DCMD_3 = [
    "python", "DAFT_Significance.py",
    "--sp","(((A,B),C),((F,G),H))",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1"
]

DCMD_2 = [
    "python", "DAFT_Significance.py",
    "--sp","(((((A,B),C),F),G),H);",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1"
]

DCMD_1 = [
    "python", "DAFT_Significance.py",
    "--sp", "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1"
]
for key, val in sim_dict.items():
    if val == "DCMD_3":
        DAFT_CMD = DCMD_3
    elif val == "DCMD_2":
        DAFT_CMD = DCMD_2
    else:
        DAFT_CMD = DCMD_1

    INPUT_CSV = "./list_" + key + ".csv"
    print(INPUT_CSV)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)

    num_chunks = (len(df) + CHUNK_SIZE - 1) // CHUNK_SIZE

    for i in range(6,7):
        start = i * CHUNK_SIZE
        end = start + CHUNK_SIZE
        chunk = df.iloc[start:end]

        chunk_csv = f"list_gt_1{i+1}.csv"
        chunk.to_csv(chunk_csv, index=False)

        out_prefix = os.path.join(OUTPUT_DIR, f"{i+1}")

        cmd = DAFT_CMD + [
            "--gt", chunk_csv,
            "--output", str(i+1)
        ]

        print(f"Running batch {i+1}/{num_chunks}...")
        subprocess.run(cmd, check=True)

    #subprocess.run(["mv", "./overall_result.csv", "./DAFT_results/overall_result.csv"], check=True)
    #subprocess.run(["mv", "./Summary.csv", "./DAFT_results/Summary.csv"], check=True)
    
    #subprocess.run(["python", "./excel_direction.py"], check=True)
    
    #base_dir = "/N/u/samishr/Quartz/Desktop/DAFT_03_23/DAFT/RESULTS/RESULTS/MARCH_25"
    #os.makedirs(base_dir, exist_ok=True)

    #dest_dir = f"{base_dir}/DAFT_results_{key}"

    #subprocess.run(["mv","./DAFT_results",dest_dir], check=True)
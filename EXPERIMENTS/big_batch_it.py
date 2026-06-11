import pandas as pd
import subprocess
import os
import shutil
import glob
import csv

def merge_summary():
    input_files = [f"./DAFT_extras_{i}/Summary_{i}.csv" for i in range(1, 11)]
    output_file = "Summary.csv"

    standard_file = input_files[0]

    if not os.path.exists(standard_file):
        raise FileNotFoundError(f"Standard header file not found: {standard_file}")

    
    
    with open(standard_file, "r", newline="") as fin:
        reader = csv.reader(fin)
        header = next(reader)

    with open(output_file, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)

        for file in input_files:
            if not os.path.exists(file):
                print(f"Warning: {file} does not exist. Skipping.")
                continue

            with open(file, "r", newline="") as fin:
                reader = csv.DictReader(fin)
                file_header = reader.fieldnames

                if file_header != header:
                    print(f"Warning: Header mismatch in {file}. Reordering by column name.")

                    missing_cols = [col for col in header if col not in file_header]
                    extra_cols = [col for col in file_header if col not in header]

                    if missing_cols:
                        print(f"  Missing columns in {file}: {missing_cols}")
                    if extra_cols:
                        print(f"  Extra columns in {file}: {extra_cols}")

                for row in reader:
                    writer.writerow([row.get(col, "") for col in header])

    print(f"Created {output_file} using header from {standard_file}")



CHUNK_SIZE = 2000
OUTPUT_DIR = "./"

sim_dict = {
    "A": "DCMD_3",
    "B": "DCMD_2",
    "C": "DCMD_1",
    "D": "DCMD_3",
    "E": "DCMD_2",
    "F": "DCMD_1",
    "G": "DCMD_3",
    "H": "DCMD_2",
    "I": "DCMD_1",
    "J": "DCMD_3",
    "K": "DCMD_2",
    "L": "DCMD_1",
    "M": "DCMD_3",
    "N": "DCMD_2",
    "O": "DCMD_2",
    "P": "DCMD_1",
    "Q": "DCMD_3",
    "R": "DCMD_2",
    "S": "DCMD_2",
    "T": "DCMD_1",
    "U": "DCMD_3",
    "V": "DCMD_2",
    "W": "DCMD_2",
    "X": "DCMD_1",
    "Y": "DCMD_1",
    "Z": "DCMD_1",
    "AA": "DCMD_2",
    "BB": "DCMD_2",
    "CC": "DCMD_3",
    "DD": "DCMD_2",
    "EE": "DCMD_2",
    "FF": "DCMD_1",
    "GG": "DCMD_3",
    "HH": "DCMD_2",
    "II": "DCMD_3",
    "JJ": "DCMD_2",
    "KK": "DCMD_1",
    "LL": "DCMD_1",
    "MM": "DCMD_1",
    "NN": "DCMD_1",
    "OO": "DCMD_1",
    "PP": "DCMD_1",
    "QQ": "DCMD_1",
    "RR": "DCMD_3",
    "SS": "DCMD_2",
    "TT": "DCMD_1",
    "Q_1": "DCMD_3",
    "R_1": "DCMD_2",
    "S_1": "DCMD_2",
    "T_1": "DCMD_1",
    "T_2": "DCMD_1",
}

DCMD_3 = [
    "daft-test",
    "--sp","(((A,B),C),((F,G),H))",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1",
    "--allow_inconsistent_rooting", "1",
                
]

DCMD_2 = [
    "daft-test",
    "--sp","(((((A,B),C),F),G),H);",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1",
    "--allow_inconsistent_rooting", "1",
]

DCMD_1 = [
    "daft-test",
    "--sp", "(A,(((W,V),(B,U)),((((K,L),(M,N)),(J,((E,(D,H)),((C,I),(G,F))))),(T,(S,(O,((Q,R),P)))))));",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1",
    "--allow_inconsistent_rooting", "1",
    
]

# Gante et al.
DCMD_4 = [
    "daft-test",
    "--sp", "(On,(Mz,(Ma,(Gr,(Br,(Pu,Ol))))));",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1",
    "--correct", "1",
    "--allow_inconsistent_rooting", "1",
    
]

# Suvorov clade 9

'''
  species_to_letters  = {
    "D_albomicans": "A",
    "D_nasuta": "B",
    "D_kepulauana": "C",
    "D_neonasuta": "D",
    "D_sulfurigaster_albostrigata": "E",
    "D_pulaua": "F",
    "D_sulfurigaster_bilimbata": "G",
    "D_sulfurigaster_sulfurigaster": "H",
    "D_neohypocausta": "I",
    "D_immigrans_kari17": "K",
    "D_immigrans": "J",
    "D_pruinosa": "L",
    "D_arawakana": "M",
    "D_dunni": "N",
    "D_cardini": "O",
    "D_ornatifrons": "P",
    "D_subbadia": "Q",
    "D_pallidipennis": "R",
    "D_funebris": "S",
    "D_guttifera": "T",
    "D_innubila": "U",
    "D_mush_saotome": "V",
    "D_quadrilineata": "W"  }
'''
DCMD_5 = [
    "daft-test",
    "--sp", "((((((((A,B),C),((D,E),(F,(G,H)))),I),(J,K)),L),(((((M,N),O),(P,Q)),R),(S,(T,(U,V))))),W);",
    "--sibling", "1",
    "--excel", "1",
    "--direction", "1",
    "--correct", "1",
    
]


for key, val in sim_dict.items():
    if val == "DCMD_3":
        DAFT_CMD = DCMD_3
    elif val == "DCMD_2":
        DAFT_CMD = DCMD_2
    
    elif val == "DCMD_4":
        DAFT_CMD = DCMD_4
    elif val == "DCMD_5":
        DAFT_CMD = DCMD_5
    else:
        DAFT_CMD = DCMD_1

    INPUT_CSV = "./DATA_SET/list_" + key + ".csv"
    print(INPUT_CSV)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)

    num_chunks = (len(df) + CHUNK_SIZE - 1) // CHUNK_SIZE

    for i in range(num_chunks):
        start = i * CHUNK_SIZE
        end = start + CHUNK_SIZE
        chunk = df.iloc[start:end]

        chunk_csv = f"list_gt_1{i+1}.csv"
        chunk.to_csv(chunk_csv, index=False)

        out_prefix = os.path.join(OUTPUT_DIR, f"{i+1}")
        
        cmd = DAFT_CMD + [
                "--gt", chunk_csv,
                "--demography", './DATA_SET/demography_' + key,
                "--output", str(i+1)
            ]
        if key in ['G','H','I','Q','R','S','T','T_6','T_7']:

            cmd = DAFT_CMD + [
                "--gt", chunk_csv,
                "--demography", './DATA_SET/demography_' + key,
                "--output", str(i+1),
                "--correct", "1",
                
            ]
    

        print(f"Running batch {i+1}/{num_chunks}...")
        subprocess.run(cmd, check=True)
        folder_name = key

        old_cache_dir = "./djiNNI_cache"

        extra_dir = "./DAFT_extras"
        os.makedirs(extra_dir, exist_ok=True)

        new_cache_dir = os.path.join(
            extra_dir,
            f"djiNNI_cache{i}_{folder_name}"
        )

        if os.path.exists(old_cache_dir):
            if os.path.exists(new_cache_dir):
                shutil.rmtree(new_cache_dir)

            shutil.move(old_cache_dir, new_cache_dir)
            print(f"Moved {old_cache_dir} to {new_cache_dir}")
        else:
            print(f"No cache folder found: {old_cache_dir}")
        
    #exit()

    merge_summary()
    #subprocess.run(["mv", "./overall_result.csv", "./DAFT_results/overall_result.csv"], check=True)
    os.makedirs("./DAFT_results", exist_ok=True)
    base_dir = "./RESULTS/RESULTS/APRIL_19"
    os.makedirs(base_dir, exist_ok=True)
    
    dest_dir = f"{base_dir}/DAFT_results_{key}"

    prefixes = ["DAFT_results", "DAFT_extras", "DAFT_log"]

    for prefix in prefixes:
        for i in range(1, 11):
            src = f"./{prefix}_{i}"
            dst = f"./{dest_dir}/{prefix}_{i}"

            if os.path.exists(src):
                shutil.move(src, dst)
            else:
                print(f"Warning: {src} does not exist")

    subprocess.run(
        ["mv", "./Summary.csv", f"./{dest_dir}/Summary.csv"],
        check=True
    )

  




    '''
    subprocess.run([
        "mv",
        "./DAFT_results",
        dest_dir
    ], check=True)
    
    '''
    
    
    
    

'''
sim_dict = {
    "J": "DCMD_3",
    "K": "DCMD_2",
    "L": "DCMD_1",
    "CC": "DCMD_3",
    "DD": "DCMD_2",
    "EE": "DCMD_2",
    "FF": "DCMD_1",
    
    
sim_dict = {
    "suvorov9": "DCMD_5"
}
}

'''
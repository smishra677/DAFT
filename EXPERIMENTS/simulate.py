import msprime
import sys
import os
sys.path.append("./DAFT_utils")
sys.path.append("./DAFT_utils/reconcILS")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd


num_gene_trees = 20000
sequence_length = 1
data_path = "./DATA_SET"


tree_small_balanced = "(((A:300.0,B:300.0)D:1200,C:1500)E:1200,((F:300.0,G:300.0)I:1100,H:1400)J:1300)K"

tree_small_unbalanced = "(((((A:300.0,B:300.0)D:1200,C:1500)E:1200,F:2700)I:1200,G:3900.0)J:1200,H:5100)K"

tree_test = "(L:4500,(((A:300.0,B:300.0)D:300,C:600)E:300,((F:300.0,G:300.0)I:300,H:600)J:300)K:3600)M"

tree_big = "(A:10282.960000000001,(((W:2999.210000000001,V:2999.210000000001)N1:2999.19,(B:2999.210000000001,U:2999.210000000001)N2:2999.19)N3:2999.19,((((K:2142.29,L:2142.29)N4:2142.2799999999997,(M:2142.29,N:2142.29)N5:2142.2799999999997)N6:2142.2799999999997,(J:5141.48,((E:2570.74,(D:1285.37,H:1285.37)N7:1285.37)N8:1285.37,((C:1285.37,I:1285.37)N9:1285.37,(G:1285.37,F:1285.37)N10:1285.37)N11:1285.37)N12:1285.37)N13:1285.37)N14:1285.37,(T:6426.85,(S:5141.48,(O:3856.11,((Q:1285.37,R:1285.37)N15:1285.37,P:2570.74)N16:1285.37)N17:1285.37)N18:1285.37)N19:1285.37)N20:1285.37)N21:1285.37)N22;"


tree_dict = {
    "small_balanced": tree_small_balanced,
    "small_unbalanced": tree_small_unbalanced,
    "test": tree_test,
    "big": tree_big,
}

sample_dict = {
    "small_balanced": {"A": 1, "B": 1, "C": 1, "F": 1, "G": 1, "H": 1},
    "small_unbalanced": {"A": 1, "B": 1, "C": 1, "F": 1, "G": 1, "H": 1},
    "test": {"A": 1, "B": 1, "C": 1, "F": 1, "G": 1, "H": 1, "L": 1},
    "big": {
        "A": 1, "B": 1, "C": 1, "D": 1, "E": 1, "F": 1, "G": 1, "H": 1, "I": 1,
        "J": 1, "K": 1, "L": 1, "M": 1, "N": 1, "O": 1, "P": 1, "Q": 1, "R": 1,
        "S": 1, "T": 1, "U": 1, "V": 1, "W": 1,
    },
}


simulations = [
    {
        "out": "A",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
        ],
    },
    {
        "out": "AA",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
            [30, 5, 0, 0.02],
        ],
    },
    {
        "out": "B",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
        ],
    },
    {
        "out": "BB",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
            [270, 5, 0, 0.02],
        ],
    },
    {
        "out": "C",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
        ],
    },
    {
        "out": "D",
        "tree_key": "small_balanced",
        "population_size": 100,
        "mass_migrations": [
        ],
    },
    {
        "out": "E",
        "tree_key": "small_unbalanced",
        "population_size": 100,
        "mass_migrations": [
        ],
    },
    {
        "out": "F",
        "tree_key": "big",
        "population_size": 100,
        "mass_migrations": [
        ],
    },
    {
        "out": "G",
        "tree_key": "small_balanced",
        "population_size": 10000,
        "mass_migrations": [
        ],
    },
    {
        "out": "GG",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
            [750, 8, 2, 0.02],
        ],
    },
    {
        "out": "H",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
        ],
    },
    {
        "out": "HH",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
            [750, 6, 5, 0.02],
        ],
    },
    {
        "out": "I",
        "tree_key": "big",
        "population_size": 10000,
        "mass_migrations": [
        ],
    },
    {
        "out": "II",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
            [750, 2, 8, 0.02],
        ],
    },
    {
        "out": "I_2",
        "tree_key": "big",
        "population_size": 10000000,
        "mass_migrations": [
            [10, 5, 0, 0.02],
        ],
    },
    {
        "out": "I_3",
        "tree_key": "big",
        "population_size": 1000000,
        "mass_migrations": [
        ],
    },
    {
        "out": "JJ",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
            [750, 5, 6, 0.02],
        ],
    },
    {
        "out": "KK",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [2357, 27, 31, 0.02],
        ],
    },
    {
        "out": "LL",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [4499, 39, 34, 0.02],
        ],
    },
    {
        "out": "M",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
            [150, 3, 0, 0.02],
        ],
    },
    {
        "out": "MM",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [1285, 17, 22, 0.02],
            [1928, 31, 10, 0.02],
        ],
    },
    {
        "out": "MM_1",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
            [145, 3, 5, 0.02],
            [150, 2, 1, 0.02],
        ],
    },
    {
        "out": "MM_4",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
            [150, 2, 1, 0.02],
            [155, 3, 5, 0.02],
        ],
    },
    {
        "out": "M_4",
        "tree_key": "small_balanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 3, 0, 0.05],
        ],
    },
    {
        "out": "N",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
            [150, 5, 0, 0.02],
        ],
    },
    {
        "out": "NN",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [643, 12, 13, 0.02],
            [643, 18, 20, 0.02],
        ],
    },
    {
        "out": "O",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
            [150, 0, 5, 0.02],
        ],
    },
    {
        "out": "OO",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [642.7, 15, 9, 0.02],
            [1071, 9, 5, 0.02],
        ],
    },
    {
        "out": "OO_1",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
            [150, 3, 5, 0.2],
            [1450, 6, 9, 0.2],
        ],
    },
    {
        "out": "P",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [643, 12, 13, 0.02],
        ],
    },
    {
        "out": "PP",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [630, 17, 20, 0.02],
            [642.7, 19, 20, 0.02],
        ],
    },
    {
        "out": "Q",
        "tree_key": "small_balanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 3, 0, 0.02],
        ],
    },
    {
        "out": "QQ",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [643, 12, 13, 0.02],
            [4499, 26, 34, 0.02],
        ],
    },
    {
        "out": "QQ_1",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [642, 20, 1, 0.02],
            [3642, 23, 27, 0.02],
        ],
    },
    {
        "out": "Q_4",
        "tree_key": "small_balanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 3, 0, 0.077],
        ],
    },
    {
        "out": "Q_6",
        "tree_key": "small_balanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 3, 0, 0.1],
        ],
    },
    {
        "out": "R",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 5, 0, 0.02],
        ],
    },
    {
        "out": "RR",
        "tree_key": "small_balanced",
        "population_size": 1000,
        "mass_migrations": [
            [150, 0, 3, 0.02],
            [151, 3, 0, 0.02],
        ],
    },
    {
        "out": "R_1",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [10, 5, 0, 0.02],
        ],
    },
    {
        "out": "R_2",
        "tree_key": "small_unbalanced",
        "population_size": 10000000,
        "mass_migrations": [
            [10, 5, 0, 0.02],
        ],
    },
    {
        "out": "R_4",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 5, 0, 0.077],
        ],
    },
    {
        "out": "R_5",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 5, 0, 0.077],
        ],
    },
    {
        "out": "R_6",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 5, 0, 0.1],
        ],
    },
    {
        "out": "S",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 0, 5, 0.02],
        ],
    },
    {
        "out": "SS",
        "tree_key": "small_unbalanced",
        "population_size": 1000,
        "mass_migrations": [
            [150, 0, 5, 0.02],
            [151, 5, 0, 0.02],
        ],
    },
    {
        "out": "S_4",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 0, 5, 0.077],
        ],
    },
    {
        "out": "S_6",
        "tree_key": "small_unbalanced",
        "population_size": 10000,
        "mass_migrations": [
            [150, 0, 5, 0.1],
        ],
    },
    {
        "out": "T",
        "tree_key": "big",
        "population_size": 10000,
        "mass_migrations": [
            [643, 12, 13, 0.02],
        ],
    },
    {
        "out": "TT",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [2356, 27, 31, 0.02],
            [2357, 31, 27, 0.02],
        ],
    },
    {
        "out": "T_1",
        "tree_key": "big",
        "population_size": 10000,
        "mass_migrations": [
            [643, 12, 13, 0.5],
        ],
    },
    {
        "out": "T_3",
        "tree_key": "big",
        "population_size": 1000000,
        "mass_migrations": [
            [643, 12, 13, 0.5],
        ],
    },
    {
        "out": "T_5",
        "tree_key": "big",
        "population_size": 10000,
        "mass_migrations": [
            [643, 12, 13, 0.082],
        ],
    },
    {
        "out": "T_6",
        "tree_key": "big",
        "population_size": 10000,
        "mass_migrations": [
            [643, 12, 13, 0.1],
        ],
    },
    {
        "out": "T_7",
        "tree_key": "big",
        "population_size": 10000,
        "mass_migrations": [
            [643, 21, 13, 0.02],
        ],
    },
    {
        "out": "U",
        "tree_key": "small_balanced",
        "population_size": 100,
        "mass_migrations": [
            [150, 3, 0, 0.02],
        ],
    },
    {
        "out": "V",
        "tree_key": "small_unbalanced",
        "population_size": 100,
        "mass_migrations": [
            [150, 5, 0, 0.02],
        ],
    },
    {
        "out": "W",
        "tree_key": "small_unbalanced",
        "population_size": 100,
        "mass_migrations": [
            [150, 0, 5, 0.02],
        ],
    },
    {
        "out": "X",
        "tree_key": "big",
        "population_size": 100,
        "mass_migrations": [
            [643, 12, 13, 0.02],
        ],
    },
    {
        "out": "Y",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [643, 16, 13, 0.02],
        ],
    },
    {
        "out": "Z",
        "tree_key": "big",
        "population_size": 1000,
        "mass_migrations": [
            [643, 7, 13, 0.02],
        ],
    },
    {
        "out": "test",
        "tree_key": "test",
        "population_size": 100,
        "mass_migrations": [
            [200, 3, 0, 0.1],
        ],
    },
]


def write_sp_file(tree,out):
    red = readWrite.readWrite()
    gene_tree_string = tree.replace("e-", "0")
    sp = red.parse_bio(gene_tree_string)

    oee_sp = open(data_path + "/sp_" + out, "w+")
    oee_sp.write(str(sp.to_newick()))
    oee_sp.close()


def write_demography_file(demography,out):
    demography.sort_events()
    oee = open(data_path + "/demography_" + out, "w+")
    oee.write(str(demography))
    oee.close()


def simulate_one(item):
    out = item["out"]
    tree_key = item["tree_key"]
    tree = tree_dict[tree_key]
    samples1 = sample_dict[tree_key]

    initial_size = {}
    for key,_ in samples1.items():
        initial_size[key] = item["population_size"]

    if tree_key == "small_balanced" or tree_key == "small_unbalanced":
        for key in ["D", "E", "I", "J", "K"]:
            initial_size[key] = item["population_size"]

    if tree_key == "test":
        for key in ["D", "E", "I", "J", "K", "M"]:
            initial_size[key] = item["population_size"]

    if tree_key == "big":
        for key in [
            "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10", "N11", "N12", "N13", "N14",
            "N15", "N16", "N17", "N18", "N19", "N20", "N21", "N22",
        ]:
            initial_size[key] = item["population_size"]

    demography = msprime.Demography.from_species_tree(tree, initial_size)

    for migration in item["mass_migrations"]:
        demography.add_mass_migration(
            time=migration[0],
            source=migration[1],
            dest=migration[2],
            proportion=migration[3],
        )

    write_sp_file(tree,out)
    write_demography_file(demography,out)

    print("Running " + out)
    print(demography)
    print(samples1)

    red = readWrite.readWrite()
    list_ = {"gt": []}

    for k in range(num_gene_trees):
        tree_sequence = msprime.sim_ancestry(
            ploidy=1,
            demography=demography,
            samples=samples1,
            sequence_length=sequence_length,
            model=msprime.StandardCoalescent(),
        )

        ts = tree_sequence
        id_to_name = {i: pop.name for i, pop in enumerate(demography.populations)}

        labels = {}
        for u in ts.samples():
            pop_id = ts.node(u).population
            labels[u] = id_to_name.get(pop_id, "pop" + str(pop_id))

        t = ts.first()
        newick = t.as_newick(node_labels=labels)
        print(newick)

        gene_tree_string = newick.replace("e-", "0")
        tr = red.parse_bio(gene_tree_string)
        print(tr.to_newick())

        list_["gt"] += [tr.to_newick()]

    pd.DataFrame(list_).to_csv(data_path + "/list_" + out + ".csv", index=False)
    #pd.DataFrame(list_).to_csv

def run_all():
    if not os.path.exists(data_path):
        os.makedirs(data_path)

    for item in simulations:
        simulate_one(item)


run_all()
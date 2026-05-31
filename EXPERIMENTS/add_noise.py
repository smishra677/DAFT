
from ete3 import PhyloTree
from collections import Counter
import sys
sys.path.append("../DAFT_utils/reconcILS/")
sys.path.append("../DAFT_utils/")
from reconcILS import *
from utils_reconcILS import *
import pandas as pd 
import numpy as np 
import argparse


reco =reconcils()
red= readWrite.readWrite()


def enumerate_branch(tr):
    returni ={}
    stack=[tr]
    curr_index=0
    while stack:
        
        current_node=stack.pop()
        if current_node:
            if current_node.isRoot ==None and current_node.isLeaf ==None:
                if current_node.parent.leftChild ==current_node:

                    returni[curr_index]=[current_node,current_node.parent,'Left']
                else:
                    returni[curr_index]=[current_node,current_node.parent,'Right']
                curr_index=curr_index+1
            
            stack.append(current_node.leftChild)
            stack.append(current_node.rightChild)
    return returni

def parse1():
    parser = argparse.ArgumentParser(description="Add noise ")
    parser.add_argument('--file', type=str, help="Path to gene tree folder")
    parser.add_argument('--og', type=str, help="og data")
    parser.add_argument('--error', type=str, help="Error")
    parser.add_argument('--out', type=str, help="out")
    args= parser.parse_args()
    return(args)

parser = parse1()
file___ =parser.file
error = parser.error

dupcoal_location=file___
outfile =parser.out
og_tree =parser.og

poi = np.random.poisson(lam=float(error), size=20000)
#op_tree = open(dupcoal_location+'trees.tre').read().split('\n')
op_tree = pd.read_csv(dupcoal_location+'/'+og_tree)['gt'].to_list()

#print(dupcoal_location+"error_trees_"+str(error)+".tre")

errored_gene_trees= {'gt':[],'NNI':[]}
with open(dupcoal_location+'/'+outfile+".txt", "w") as file:
    for idx_,ij in enumerate(poi):
        #rando_branch = np.random.randint(ij,size=20)
        tree= str(op_tree[idx_])
        tree  = tree.replace('e-', '0')
        tr = red.parse(tree)
        tr.isRoot=True
        tr.label_internal()
        enum= enumerate_branch(tr)
        #print(len(enum),enum)
        #exit()
        if ij==0 or len(enum)==0 or ij>len(enum):
            errored_gene_trees['gt']+=[tr.to_newick()]
            errored_gene_trees['NNI']+=[0]
            file.write(f"{tr.to_newick()}\n")
            continue
        else:
            nni=0
            tosses = np.random.binomial(n=1, p=ij/len(enum), size=len(enum))
            for idx,toss in enumerate(tosses):
                if toss==1:
                    get_branch =enum[idx]
                    val =get_branch[1].NNI(tr,flag=get_branch[2])
                    toss_2 = np.random.binomial(n=1, p=1/2, size=1)
                    if toss_2==1:
                        tr = red.parse(val[1][1].to_newick())
                    else:
                        tr = red.parse(val[0][1].to_newick())
                    tr.isRoot=True
                    tr.label_internal()
                    enum= enumerate_branch(tr)
                    nni=nni+1
            

            errored_gene_trees['gt']+=[tr.to_newick()]
            errored_gene_trees['NNI']+=[nni]
            file.write(f"{tr.to_newick()}\n")


pd.DataFrame(errored_gene_trees).to_csv(dupcoal_location+'/'+outfile+".csv", index=False)


log_filename = outfile+"_log.txt"

with open(dupcoal_location+'/'+log_filename, "w") as log:
    log.write("\nParsed Arguments:\n")
    log.write(f"file: {file___}\n")
    log.write(f"og: {og_tree}\n")
    log.write(f"error: {error}\n")
    log.write(f"out: {outfile}\n")
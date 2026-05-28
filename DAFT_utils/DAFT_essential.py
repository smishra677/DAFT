import sys
import math
import re
sys.path.append("./DAFT_utils/reconcILS")
sys.path.append("./DAFT_utils/")
import ete3
import pandas as pd
import numpy as np
from utils_reconcILS import *
from reconcILS import *
import numpy as np
import warnings
from pathlib import Path



warnings.filterwarnings(
    "ignore",
    message="Downcasting object dtype arrays on .*fillna.*",
    category=FutureWarning,
)

warnings.filterwarnings(
    "ignore",
    message="DataFrameGroupBy.apply operated on the grouping columns.*",
    category=DeprecationWarning,
)
        
reco =reconcils()
red= readWrite.readWrite()
Il=ILS.ILS()

class daft_essential:
    def __init__(self):
        #self.reco =reconcils()
        #self.red= readWrite.readWrite()
        #self.Il=ILS.ILS()
        np.random.seed(42)


    def convert_gt_to_csv(self,gt_file):
        gt_file = Path(gt_file)

        if gt_file.suffix.lower() == ".csv":
            return str(gt_file)

        data = {"gt": []}

        with open(gt_file, "r", encoding="utf-8") as f:
            for line in f:
                tree = line.strip()

                if tree:
                    data["gt"].append(tree)

        if not data["gt"]:
            raise ValueError(f"Gene tree file is empty: {gt_file}")

        df = pd.DataFrame(data)

        # Save as CSV
        csv_file = gt_file.with_suffix(".csv")
        df.to_csv(csv_file, index=False)

        return str(csv_file)



    def find_sibling(self,gt,returni):
        stack=[gt]
        while stack:
            current_address=stack.pop()
            if current_address:
                if not current_address.isLeaf:
                    returni+=[[i.to_newick() for i in current_address.children]]
                
                stack.append(current_address.leftChild)
                stack.append(current_address.rightChild)
        
        return returni


    def find_uncle(self,gt,returni):
        stack=[gt]
        while stack:
            current_address=stack.pop()
            if current_address:
                if current_address.parent:
                    if current_address.parent.parent:
                        returni+=[[[current_address.to_newick() , i.to_newick()] for i in current_address.parent.parent.children if i!=current_address.parent ]]
                    
                stack.append(current_address.leftChild)
                stack.append(current_address.rightChild)
        

        
        return [i[0]  for i in returni]



    def project_z_score_z(self, dummy, flag, df):
        
        mask = dummy['comparison_' + flag].notna() & (dummy[flag + '_count'] == 0)

        dummy[flag + '_count'] = dummy[flag + '_count'].astype(float)
        dummy.loc[mask, flag + '_count'] += 5e-5

        dummy['total_count'] = dummy['total_count'].astype(float)
        dummy.loc[mask, 'total_count'] += 5e-5
        
        
        x1 = dummy[flag+'_count'].astype(float)
        x2 = dummy['total_count'].astype(float)
        
        n1 = dummy[flag+'_count_population'].astype(float)
        n2 = dummy['test_count_population'].astype(float)
        
        
        dummy['Z-value-'+flag] = np.nan  
        
        denom = np.sqrt(x1 + x2)

        mask = dummy['comparison_'+flag].notna() & (denom != 0)
        dummy.loc[mask, 'Z-value-'+flag] = (x1 - x2) / denom
        
        dummy['Z-value-'+flag+'_corrected_scaled_down'] = np.nan  
  
        
        x1_copy=x1
        x2_copy=x2
        
        scale_x1 = n2 < n1
        scale_x2 = ~scale_x1

        #x1.loc[scale_x1] = (x1.loc[scale_x1] / n1.loc[scale_x1]) * n2.loc[scale_x1]
        #x2.loc[scale_x2] = (x2.loc[scale_x2] / n2.loc[scale_x2]) * n1.loc[scale_x2]
        x1 = x1.where(~scale_x1, x1 / n1 * n2)
        x2 = x2.where(~scale_x2, x2 / n2 * n1)   
        
        
        dummy[flag+'_count_scaled_down'] = x1 
        dummy[flag+'_total_count_scaled_down'] = x2
        
        #print(x2)
        
        #dummy[flag+'_count'] = x1
        #dummy['total_count'] = x2
        
        #print(dummy['total_count'])
        denom = np.sqrt(x1 + x2)

        mask = dummy['comparison_'+flag].notna() & (denom != 0)
        dummy.loc[mask, 'Z-value-'+flag+'_corrected_scaled_down'] = (x1 - x2) / denom
        
        

    def project_z_score(self, dummy, flag, df):
        
        mask = dummy['comparison_' + flag].notna() & (dummy[flag + '_count'] == 0)

        dummy[flag + '_count'] = dummy[flag + '_count'].astype(float)
        dummy.loc[mask, flag + '_count'] += 5e-5

        dummy['total_count'] = dummy['total_count'].astype(float)
        dummy.loc[mask, 'total_count'] += 5e-5
        
        
        x1 = dummy[flag+'_count'].astype(float)
        x2 = dummy['total_count'].astype(float)
        
        n1 = dummy[flag+'_count_population'].astype(float)
        n2 = dummy['test_count_population'].astype(float)
        
        
        dummy['Z-value-'+flag] = np.nan  
        
        denom = np.sqrt(x1 + x2)

        mask = dummy['comparison_'+flag].notna() & (denom != 0)
        dummy.loc[mask, 'Z-value-'+flag] = (x1 - x2) / denom
        
        dummy['Z-value-'+flag+'_corrected_scaled_down'] = np.nan  
  
        
        x1_copy=x1
        x2_copy=x2
        
        scale_x1 = n2 < n1
        scale_x2 = ~scale_x1

        #x1.loc[scale_x1] = (x1.loc[scale_x1] / n1.loc[scale_x1]) * n2.loc[scale_x1]
        #x2.loc[scale_x2] = (x2.loc[scale_x2] / n2.loc[scale_x2]) * n1.loc[scale_x2]
        x1 = x1.where(~scale_x1, x1 / n1 * n2)
        x2 = x2.where(~scale_x2, x2 / n2 * n1)   
        
        
        dummy[flag+'_count_scaled_down'] = x1 
        dummy[flag+'_total_count_scaled_down'] = x2
        
        #print(x2)
        
        #dummy[flag+'_count'] = x1
        #dummy['total_count'] = x2
        
        #print(dummy['total_count'])
        denom = np.sqrt(x1 + x2)

        mask = dummy['comparison_'+flag].notna() & (denom != 0)
        dummy.loc[mask, 'Z-value-'+flag+'_corrected_scaled_down'] = (x1 - x2) / denom
        
        
        
        dummy['Z-value-'+flag+'_corrected_scaled_up'] = np.nan

        #x1.loc[scale_x1] = (x1.loc[scale_x1] / n1.loc[scale_x1]) * n2.loc[scale_x1]
        #x2.loc[scale_x2] = (x2.loc[scale_x2] / n2.loc[scale_x2]) * n1.loc[scale_x2]
        x1_copy = x1_copy.where(scale_x1, x1_copy / n1 * n2)
        x2_copy = x2_copy.where(scale_x2, x2_copy / n2 * n1)   
        
        
        dummy[flag+'_count_scaled_up'] = x1_copy 
        dummy[flag+'_total_count_scaled_up'] = x2_copy
        #print(x2)
        
        #dummy[flag+'_count'] = x1
        #dummy['total_count'] = x2
        
        #print(dummy['total_count'])
        denom = np.sqrt(x1_copy + x2_copy)

        mask = dummy['comparison_'+flag].notna() & (denom != 0)
        dummy.loc[mask, 'Z-value-'+flag+'_corrected_scaled_up'] = (x1_copy - x2_copy) / denom
        
        

    def calculate_pooled_z(self,dummy,flag,df):
        
        x1 = dummy[flag+'_count'].astype(float)
        x2 = dummy['total_count'].astype(float)
        
        n1 = dummy[flag+'_count_population'].astype(float)
        n2 = dummy['test_count_population'].astype(float)
        #print(n1,n2)
        
        
        
        df['Z-value-'+flag] = np.nan  
       
        
        #ummy.loc[mask, 'Z-value-'+flag] = (X - Y) / denom
       
        p1 = x1 / n1
        p2 = x2 / n2
        
      
        p_pooled = (x1 + x2) / (n1 + n2)
        
       
        #if p_pooled == 0 or p_pooled == 1:
        #return 0.0
            
        
        se = np.sqrt(p_pooled * (1 - p_pooled) * (1/n1 + 1/n2))
        

        z = (p1 - p2) / se
        denom = np.sqrt(x1 + x2)
        mask = dummy['comparison_'+flag].notna() & (denom != 0)
        dummy.loc[mask, 'Z-value-'+flag] =round(z, 4)


    def current_address(self,taxa,gene_tree):
            stack = [gene_tree] 
            
            if type(taxa)!=set:
                taxa={taxa}
            while stack:
                current_node = stack.pop()
                    
                if current_node:
                    current_taxa=current_node.taxa
                    if type(current_node.taxa)!=set:
                            current_taxa={current_taxa}
                    if current_taxa==taxa:
                            return current_node
            
                    stack.append(current_node.leftChild)
                    stack.append(current_node.rightChild)



    #https://www.geeksforgeeks.org/dsa/lowest-common-ancestor-binary-tree-set-1/
    # Function to find the level of a node
    def findLevel(self,root, k, level):
        if root is None:
            return -1
        if root.taxa == k:
            return level
        
        #print(level)
        
        # Recursively call function on left child
        leftLevel = self.findLevel(root.leftChild, k, level + 1)
        
        # If node is found on left, return level
        # Else continue searching on the right child
        if leftLevel != -1:
            return leftLevel
        else:
            return self.findLevel(root.rightChild, k, level + 1)

    # Function to find the lowest common ancestor 
    # and calculate distance between two nodes
    def findLcaAndDistance(self,root, a, b, d1, d2, dist, lvl):
        if root is None:
            return None
        
        if root.taxa == a:
        
            # If first node found, store level and 
            # return the node
            d1[0] = lvl
            return root
        if root.taxa == b:
        
            # If second node found, store level and 
            # return the node
            d2[0] = lvl
            return root

        # Recursively call function on left child
        left = self.findLcaAndDistance(root.leftChild, a, b, d1, d2, dist, lvl + 1)
    
        # Recursively call function on right child
        right = self.findLcaAndDistance(root.rightChild, a, b, d1, d2, dist, lvl + 1)

        if left is not None and right is not None:
        
            # If both nodes are found in different 
            # subtrees, calculate the distance
            dist[0] = d1[0] + d2[0] - 2 * lvl

        # Return node found or None if not found
        if left is not None:
            return left
        else:
            return right

    # Function to find distance between two nodes
    def findDist(self,root, a, b):
        d1 = [-1]
        d2 = [-1]
        dist = [0]
        
        # Find lowest common ancestor and calculate distance
        lca = self.findLcaAndDistance(root, a, b, d1, d2, dist, 1)
        
        if d1[0] != -1 and d2[0] != -1:
        
            # Return the distance if both nodes are found
            return dist[0]

        if d1[0] != -1:
        
            # If only first node is found, find 
            # distance to second node
            dist[0] = self.findLevel(lca, b, 0)
            return dist[0]
        
        if d2[0] != -1:
        
            # If only second node is found, find 
            # distance to first node
            dist[0] = self.findLevel(lca, a, 0)
            return dist[0]
        
        # Return -1 if both nodes not found
        return -1

    def isequal(self,tree1, tree2):
        if tree1==None or tree2==None:
            return False
        if tree1[-1] != ";":
            tree1 = tree1 + ";"
        if tree2[-1] != ";": 
            tree2 = tree2 + ";"
        firstTree = ete3.Tree(tree1)
        secondTree = ete3.Tree(tree2)
        node_set_firstTree = {node.name for node in firstTree.traverse() if node.name}
        node_set_secondTree = {node.name for node in secondTree.traverse() if node.name}
        if len(tree1) <= 4 or len(tree2) <= 4: 
            if node_set_firstTree == node_set_secondTree: 
                return True
            else:
                return False 
        else: 
            rf, _, _, _, _, _, _ = firstTree.robinson_foulds(secondTree)
            
            if rf == 0 and node_set_firstTree==node_set_secondTree: 
                return True 
            else: 
                return False 
            


    def find_all_lineage(self,tree):
            ret=[]
            stack=[tree]

            while stack:
                    current_node= stack.pop()
                    if current_node:
                        ret.append(current_node.to_newick())
                                    
                            
                        if current_node.leftChild:
                                    stack.append(current_node.leftChild)
                        if current_node.rightChild:
                                    stack.append(current_node.rightChild)
            return ret
            
    def isequal_set(self,tree1,tree2):
        if tree1==None or tree2 ==None:
            return False
        Tree1= red.parse(tree1)
        Tree1.label_internal()
        Tree2= red.parse(tree2)
        
        Tree2.label_internal()

        if Tree1.taxa==Tree2.taxa:
            return True
        else:
            return False

    def match_lineage(self,val, lineage_list):
        #print(val,lineage_list)
        for lineage in lineage_list:
            if self.isequal_set(val, lineage):
                return lineage   
        return val 


        
    def find_pairs(self,df):
        pairs=[]
        for what_moved, L1, L2,Z_1,Z_2 in zip(df['What_moved'],df['Lineage1'],df['Lineage2'],df['Z_score'],df['Z_score_sibling']):
            if self.isequal_set(what_moved,L1):
                pairs+=[[what_moved,L2]]
            else:
                pairs+=[[what_moved,L1]]
            #print(Z_1,Z_2)
            Z_1 = 0.0 if pd.isna(Z_1) else Z_1
            Z_2 = 0.0 if pd.isna(Z_2) else Z_2
            #if Z_1>=-4:
            if Z_2<-1.96:
                    pairs1 =pairs[-1]
                    pairs+=[[pairs1[1],pairs1[0]]]


        return pairs



    def write_events_network(self,tree):
            if tree.isLeaf:
                ev=''
                if len(tree.donor)==0 and len(tree.receiver)==0:
                    return ev+f"{tree.taxa}[&label={tree.id_map_tag}]"

                for idx,event in enumerate(tree.donor):
                    if idx==0:
                        ev= '('+f"{tree.taxa}[&label={tree.id_map_tag}]"+')'+event
                    else:
                        ev='('+ev+')'+event
                for idx ,event  in  enumerate(tree.receiver):
                    if idx==0 and len(tree.donor)==0:
                        ev= '('+event+','+f"{tree.taxa}[&label={tree.id_map_tag}]"+')'
                    else:
                        #print('==========>')
                        ev='('+ev+','+event+')'
                return ev  
            else:
                return ''

    
    def check_internal(self,newick,tree):
        if tree.isLeaf==None:
            for idx,event in enumerate(tree.donor):
                if idx==0:
                    newick= '('+newick+')'+event
                else:
                    newick='('+newick+')'+event

            for idx ,event  in  enumerate(tree.receiver):
                if idx==0 and len(tree.donor)==0:
                    newick= '('+event+','+newick+')'
                else:
                    newick='('+newick+','+event+')'
        
        return newick

        #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
        # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    def traverse_network(self,newick,tree):
            if tree.leftChild and not tree.rightChild:
                newick = f"(,{self.traverse_network(newick,tree.leftChild)}){self.write_events_network(tree)}"
            elif not tree.leftChild and tree.rightChild:
                newick = f"({self.traverse_network(newick,tree.rightChild)},){self.write_events_network(tree)}"
            elif tree.leftChild and tree.rightChild:
                newick = f"({self.traverse_network(newick,tree.rightChild)},{self.traverse_network(newick,tree.leftChild)}){self.write_events_network(tree)}"
            elif not tree.leftChild and not tree.rightChild :
                newick = f"{self.write_events_network(tree)}"
            else:
                pass
            newick=self.check_internal(newick,tree)
            return newick





        # Got it From StackOverflow:
        #https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
        # https://stackoverflow.com/questions/61117131/how-to-convert-a-binary-tree-to-a-newick-tree-using-python
    def to_network(self,tree):
            newick = ""
            newick = self.traverse_network(newick,tree)
            newick = f"{newick};"
            return newick



    def put_network_in_tree(self,df,sp):

        pairs= self.find_pairs(df)


        #sp = red.parse(tree)
        sp.label_internal()
        i=1
        current_event = '#H'+str(i)
        for donor, reciver in pairs:
            l1 = red.parse(donor)
            l2= red.parse(reciver)
            l1.label_internal()
            l2.label_internal()
            taxa_1=l1.taxa
            taxa_2=l2.taxa
            l1_address= self.current_address(taxa_1,sp)
            l2_address= self.current_address(taxa_2,sp)
            
            l1_address.donor+=[current_event]
            l2_address.receiver+=[current_event]


            i=i+1
            current_event = '#H'+str(i)


        return sp

    def find_dist_string(self,sp,tree1,tree2):
        Tree1= red.parse(tree1)
        Tree1.label_internal()
        Tree2= red.parse(tree2)
        
        Tree2.label_internal()

        return self.findDist(sp,Tree1.taxa,Tree2.taxa)-2


    def _is_sibling(self,lineage1,lineage2,sp_string):
        l1=red.parse(lineage1)
        l2=red.parse(lineage2)
        l1.label_internal()
        l2.label_internal()
        sp=red.parse(sp_string)
        sp.label_internal()

        curr_1=self.current_address(l1.taxa,sp)
        if curr_1:
            if curr_1.parent:
                other_child =[chil for chil in curr_1.parent.children if chil!=curr_1][0]
                if other_child.taxa==l2.taxa:
                    return True
                else:
                    return False

    def _is_in_unique(self,lineage1,lineage2,lineagesU):
        for l1,l2 in lineagesU:
            if (self.isequal_set(l1,lineage1) and self.isequal_set(l2,lineage2)) or (self.isequal_set(l2,lineage1) and self.isequal_set(l1,lineage2)):
                return True
        return False


    def rename_leaves(self,newick, rename_dict):
        newick= red.parse(newick).to_newick()
        for old, new in rename_dict.items():
            newick = re.sub(rf'(?<=\(|,){old}(?=[:,)\s])', f'{old} {new}', newick)
                    
        return newick
                
    def rename_subtrees(self,newick, rename_dict):
        newick= red.parse(newick).to_newick()
        for old, new in rename_dict.items():
            if old in newick:
                newick = newick.replace(old, f"{old} {new}")
        return newick

    def sibling_bidirectional(self,df,lineagesU,sp_string):
        df = df.copy()

        #for c in BASE_COLS:
        df['Significant_Pairs'] = pd.NA
        df['minor_sibling_Total_gene_trees'] = pd.NA
        df['minor_sibling'] = pd.NA
        df['minor_sibling_count'] = pd.NA
        df['minor_sibling_What_moved'] = pd.NA
        df['minor_sibling_To_where'] = pd.NA
        df['Z_score'] = pd.NA
        df['Z_score_sibling'] = pd.NA
        df['inUnique'] = False
        df['CouldbeBidirectional'] = False



        idxs = list(df.index)

        for i in idxs:
            L1_i = df.at[i, 'Lineage1']
            L2_i = df.at[i, 'Lineage2']
            df.at[i,'Significant_Pairs'] = L1_i+' AND '+L2_i

            if self._is_in_unique(L1_i,L2_i,lineagesU):
                df.at[i,'inUnique'] =True


            L1_count_i = df.at[i, 'Count1']
            L2_count_i = df.at[i, 'Count2']
            minZ_i=min(L1_count_i,L2_count_i)
            maxZ_i=max(L1_count_i,L2_count_i)
            
            denomi =L1_count_i+L2_count_i
            if denomi ==0:
                Z_score=np.nan
            else:
                Z_score= (minZ_i-maxZ_i)/ math.sqrt(L1_count_i+L2_count_i)
            #if Z_score<-0.4:
            #continue
            df.at[i,'Z_score'] = Z_score
            
            
            found_j = None
            found_k = None
            for j in idxs:
                if j == i:
                    #print('yess')
                    continue

                L1_j = df.at[j, 'Lineage1']
                L2_j = df.at[j, 'Lineage2']
                L1_count_j = df.at[j, 'Count1']
                L2_count_j = df.at[j, 'Count2']

                if self.isequal_set(L1_i, L1_j) and self._is_sibling(L2_i, L2_j,sp_string) and L2_count_i<=L1_count_i:
                    found_j = j
                    found_k = 1
                    break

                if self.isequal_set(L1_i, L2_j) and self._is_sibling(L2_i, L1_j,sp_string) and L2_count_i<=L1_count_i:
                    found_j = j
                    found_k = 2
                    break

                if self.isequal_set(L2_i, L1_j) and self._is_sibling(L1_i, L2_j,sp_string) and L1_count_i<=L2_count_i:
                    found_j = j
                    found_k = 1
                    break

                if self.isequal_set(L2_i, L2_j) and self._is_sibling(L1_i, L1_j,sp_string) and L1_count_i<=L2_count_i:
                    found_j = j
                    found_k = 2
                    break

            if found_j is not None:
                df.at[i,'minor_sibling_Total_gene_trees'] = df.at[found_j, 'Total_gene_trees']
                df.at[i,'minor_sibling_What_moved'] = df.at[found_j, 'What_moved']
                df.at[i,'minor_sibling_To_where'] = df.at[found_j, 'To_where']
                

                if found_k==1:
                    df.at[i,'minor_sibling'] = df.at[found_j, 'Lineage2']
                    df.at[i,'minor_sibling_count'] = df.at[found_j, 'Count2']
                    denomi = minZ_i+L2_count_j
                    
                    if denomi==0:
                        Z_score_sib=np.nan
                    else:
                        Z_score_sib= (L2_count_j-minZ_i)/ math.sqrt(minZ_i+L2_count_j)

                else:
                    df.at[i,'minor_sibling'] = df.at[found_j, 'Lineage1']
                    df.at[i,'minor_sibling_count'] = df.at[found_j, 'Count1']
                    denomi= minZ_i+L1_count_j
                    if denomi==0:
                        Z_score_sib=np.nan
                    else:
                        Z_score_sib= (L1_count_j-minZ_i)/ math.sqrt(minZ_i+L1_count_j)

                df.at[i,'Z_score_sibling']=Z_score_sib

                if Z_score_sib<=-1.96:
                    df.at[i,'CouldbeBidirectional']=True

        return df

    # ---------- styling----------
    def render_val_2(self,col, v):
        if col in ("Z_score", "Z_score_sibling"):
            return "-" if pd.isna(v) else f"{float(v):.2f}"
        elif col in ("Total_gene_trees", "Count1", "Count2",
                    "minor_sibling_Total_gene_trees", "minor_sibling_count"):
            if pd.isna(v):
                return "-"
            try:
                fv = float(v)
                return str(int(fv)) if fv.is_integer() else str(fv)
            except Exception:
                return str(v)
        else:
            return "" if v is None or (isinstance(v, float) and math.isnan(v)) else str(v)
        
    
    def render_val_1(self,col, v):
        if col in ("Z_score", "Z_score_sibling"):
            return "-" if pd.isna(v) else f"{float(v):.2f}"
        elif col in (
            "Total_gene_trees",
            "Count1",
            "Count2",
            "minor_sibling_Total_gene_trees",
            "minor_sibling_count",
        ):
            if pd.isna(v):
                return "-"
            try:
                fv = float(v)
                return str(int(fv)) if fv.is_integer() else str(fv)
            except Exception:
                return str(v)
        else:
            return "" if v is None or (isinstance(v, float) and math.isnan(v)) else str(v)


    def compute_widths(self,frame, cols):
        widths = {}
        for c in cols:
            header_len = len(c)
            if c in frame.columns:
                data_lens = frame[c].map(lambda v: len(self.render_val(c, v)))
            else:
                data_lens = pd.Series([0])
            widths[c] = max([header_len] + data_lens.tolist())
        return widths


    ##  Turn list into string
    def as_arg_list(self,x, sep=","):
        if isinstance(x, (list, tuple, set)):
            return sep.join(map(str, x))
        return str(x)

    def map_pair_uncle(self,what_moved,Elist):
        #print(Elist)
        #print('###')
        #print(what_moved)

        for pair in Elist:

            key=pair[0]
            uncle =pair[1]
            #print(key,uncle)
        


            if self.isequal_set(what_moved,key):
                #print(what_moved)
                return uncle
            #if self.isequal_set(what_moved,uncle):
            #return key
        for pair in Elist:

            key=pair[0]
            uncle =pair[1]
            #print(key,uncle)
        


            #if self.isequal_set(what_moved,key):
                #print(what_moved)
                #return uncle
            if self.isequal_set(what_moved,uncle):
                return key
        return None  

    def map_pair_sibling(self,what_moved,Elist):
        #print(Elist)
        #print('###')
        #print(what_moved)
        for pair in Elist:

            key=pair[0]
            uncle =pair[1]
            #print(key,uncle)
        


            if self.isequal_set(what_moved,key):
                #print(what_moved)
                return uncle
            if self.isequal_set(what_moved,uncle):
                return key
        return None  


    def is_in(self,tree1,tree2):
        if tree1==None or tree2 ==None:
            return False
        Tree1= red.parse(tree1)
        Tree1.label_internal()
        Tree2= red.parse(tree2)
        Tree2.label_internal()


        if Tree1.isLeaf:
            #print(Tree1.taxa)
            Tree1.taxa={Tree1.taxa}
        if Tree2.isLeaf:
            #print(Tree2.taxa)
            Tree2.taxa={Tree2.taxa}

        if Tree1.taxa.intersection(Tree2.taxa) == Tree1.taxa:
            return True
        else:
            return False

    def compute_z(self,dummy,flag,df):
        X = dummy[flag+'_count'].astype(float)
        Y = dummy['total_count'].astype(float)
        denom = np.sqrt(X + Y)
        df['Z-value-'+flag] = np.nan  
        mask = dummy['comparison_'+flag].notna() & (denom != 0)
        dummy.loc[mask, 'Z-value-'+flag] = (X - Y) / denom






    def render_val(self,col, v):
        if col in ('Z-value-uncle', 'Z-value-sibling','Z-value-sibling_corrected_scaled_down','Z-value-uncle_corrected_scaled_down'):
            return '-' if pd.isna(v) or v=='' else f"{float(v):.2f}"
        elif col in ('total_count', 'uncle_count', 'sibling_count','Z_score_sibling','Minor_moved_count'):
            if pd.isna(v) or v =='':
                return '-'
            try:
                fv = float(v)
                return str(int(fv)) if fv.is_integer() else str(fv)
            except Exception:
                return str(v)
        elif col in ('minor_sibling_count', 'minor_sibling_Total_gene_trees', 'minor_sibling'):
            if pd.isna(v) or v =='':
                return '-'
        else:
            return '' if v is None or (isinstance(v, float) and math.isnan(v)) else str(v)

    def _fix_group_nested(self,df):
        what_list = df['What_moved'].tolist()
        for idx, cu in df['comparison_uncle'].items():
            for wm in what_list:
                if (not self.isequal(cu, wm)) and self.isequal_set(cu, wm):
                    df.at[idx, 'comparison_uncle'] = wm
                    break  
        for idx, cu in df['comparison_sibling'].items():
            for wm in what_list:
                if (not self.isequal(cu, wm)) and self.isequal_set(cu, wm):
                    df.at[idx, 'comparison_sibling'] = wm
                    break  
        return df

    def _uncle_count(self,row,grouped):
        mask1 = grouped['What_moved'].apply(lambda x: self.isequal_set(x, row['comparison_uncle'])) 
        mask2 = grouped['Where_at'].apply(lambda x: self.isequal_set(x, row['Where_at']))
        
        if mask2.any() and not mask1.any():
            return pd.NA          
        return grouped.loc[mask1&mask2, 'total_count'].sum()


    def _sibling_population_count(self,row,grouped):
        siblings = row['comparison_sibling']
        mask = grouped['Where_at'].apply(lambda x: self.isequal_set(x, siblings))
        total = grouped.loc[mask, 'total_count'].sum()
        dfer =pd.read_csv('./rev_all_corrected.csv')
        
        mask1 = dfer['Topo'].apply(lambda x: self.isequal_set(x, siblings))
        total1 = dfer.loc[mask1, 'junk_count'].sum()
        #total1=0
        return (total+total1) if (total > 0) else pd.NA
    
    
    def _uncle_population_count(self,row,grouped):
        uncle = row['comparison_uncle']
        mask = grouped['Where_at'].apply(lambda x: self.isequal_set(x, uncle))
        total = grouped.loc[mask, 'total_count'].sum()
        dfer =pd.read_csv('./rev_all_corrected.csv')
        
        mask1 = dfer['Topo'].apply(lambda x: self.isequal_set(x, uncle))
        total1 = dfer.loc[mask1, 'junk_count'].sum()
        #total1=0
        return (total+total1) if (total > 0) else pd.NA
    
    def _test_population_count(self,row,grouped):
        test = row['What_moved']
        mask = grouped['Where_at'].apply(lambda x: self.isequal_set(x, test))
        total = grouped.loc[mask, 'total_count'].sum()
        
        dfer =pd.read_csv('./rev_all_corrected.csv')
        mask1 = dfer['Topo'].apply(lambda x: self.isequal_set(x, test))
        total1 = dfer.loc[mask1, 'junk_count'].sum()
        #total1=0
        return (total+total1) if (total > 0) else pd.NA
    
    
    def _sibling_count(self,row,grouped):
        mask1 = grouped['What_moved'].apply(lambda x: self.isequal_set(x, row['comparison_sibling'])) 
        mask2 = grouped['Where_at'].apply(lambda x: self.isequal_set(x, row['Where_at']))
        
        if mask2.any() and not mask1.any():
            return pd.NA          
        return grouped.loc[mask1&mask2, 'total_count'].sum()

    def _is_uncle(self,lineage1,lineage2,sp_string):
        l1=red.parse(lineage1)
        l2=red.parse(lineage2)
        l1.label_internal()
        l2.label_internal()
        sp=red.parse(sp_string)
        sp.label_internal()

        curr_1=self.current_address(l1.taxa,sp)
        if curr_1:
            if curr_1.parent:
                if curr_1.parent.parent:
                    other_child =[chil for chil in curr_1.parent.parent.children if chil!=curr_1.parent][0]
                    if other_child.taxa==l2.taxa:
                        return True
                    else:
                        return False


    def id_it(self,sp_tree):
        sp =red.parse(sp_tree)
        id_map_node={}
        id_map_branch={'From':[],'To':[],'id':[]}
        counter_node=0
        counter_branch=0
        stack=[sp]
        while stack:
            current_node = stack.pop()
            if current_node:
                current_node.id_map_tag=counter_branch
                if current_node.parent:
                    id_map_branch['From']+=[current_node.parent.to_newick()]
                    id_map_branch['To']+=[current_node.to_newick()]
                    id_map_branch['id']+=[counter_branch]
                else:
                    id_map_branch['From']+=[current_node.to_newick()]
                    id_map_branch['To']+=[current_node.to_newick()]
                    id_map_branch['id']+=[counter_branch]

                id_map_node[current_node.to_newick()]=counter_node
                counter_node=counter_node+1
                counter_branch = counter_branch+1
                stack.append(current_node.leftChild)
                stack.append(current_node.rightChild)
        
        
        return id_map_node,id_map_branch,sp

    def map_id(self,val,id_map_node):
        for key,value in id_map_node.items():
            if self.isequal(val,key):
                return str(value)
        return None

    def idfy_it(self,dfc,id_map_node):
        df = dfc.copy()
        for idx,row in df.iterrows():
            old_what_moved = row['What_moved']
            old_where_at = row['Where_at']
            old_comp_sib = row['comparison_sibling']
            old_comp_uncle = row['comparison_uncle']
            

            new_what_moved =self.map_id(old_what_moved,id_map_node)
            new_where_at = self.map_id(old_where_at,id_map_node)
            
            new_comp_sib=old_comp_sib
            new_comp_uncle= old_comp_uncle
            if old_comp_sib:
                
                new_comp_sib =self.map_id(old_comp_sib,id_map_node)
                
            if old_comp_uncle:
                new_comp_uncle =self.map_id(old_comp_uncle,id_map_node)

            
            
            
            df.loc[idx,'What_moved']=new_what_moved
            df.loc[idx,'Where_at']=new_where_at
            df.loc[idx,'comparison_sibling']=new_comp_sib
            df.loc[idx,'comparison_uncle']=new_comp_uncle
        return df
            

    def idfy_it_direction(self,dfc,id_map_node):
        df= dfc.copy()
        for idx, row in df.iterrows():

            old_lineage1 = row['Lineage1']
            old_lineage2 = row['Lineage2']
            old_what_moved = row['What_moved']
            old_to_where = row['To_where']
            
            #print(old_lineage1,old_lineage2)

            # ---- Minor sibling fields ----
            old_minor_moved = row['Minor_Moved']
            old_minor_sibling = row['minor_sibling']
            old_minor_sib_what = row['minor_sibling_What_moved']
            old_minor_sib_to_where = row['minor_sibling_To_where']

            # ---- Significant pairs ----
            old_sig = row['Significant_Pairs']
            old_receiver, old_donor = '', ''

            if pd.notna(old_sig) and 'AND' in old_sig:
                parts = old_sig.split('AND')
                old_receiver = parts[0].strip()
                old_donor = parts[1].strip()

            # ---- Helper ----
            def convert_tree(val):
                if pd.notna(val) and val != '':
                    return self.map_id(val,id_map_node)
                return val

            # ---- Convert primary ----
            new_lineage1 = convert_tree(old_lineage1)
            new_lineage2 = convert_tree(old_lineage2)
            new_what_moved = convert_tree(old_what_moved)
            new_to_where = convert_tree(old_to_where)

            # ---- Convert significant ----
            new_receiver = convert_tree(old_receiver)
            new_donor = convert_tree(old_donor)

            # ---- Convert minor ----
            new_minor_moved = convert_tree(old_minor_moved)
            new_minor_sibling = convert_tree(old_minor_sibling)
            new_minor_sib_what = convert_tree(old_minor_sib_what)
            new_minor_sib_to_where = convert_tree(old_minor_sib_to_where)

            # ---- Store back (ALWAYS write back) ----
            df.at[idx, 'Lineage1'] = new_lineage1
            df.at[idx, 'Lineage2'] = new_lineage2
            df.at[idx, 'What_moved'] = new_what_moved
            df.at[idx, 'To_where'] = new_to_where

            df.at[idx, 'Minor_Moved'] = new_minor_moved
            df.at[idx, 'minor_sibling'] = new_minor_sibling
            df.at[idx, 'minor_sibling_What_moved'] = new_minor_sib_what
            df.at[idx, 'minor_sibling_To_where'] = new_minor_sib_to_where

            df.at[idx, 'Significant_Pairs'] = str(new_receiver) + ' AND ' + str(new_donor)

        return df


    def to_newick_with_id(self,labeled_sp):
        newick = ""
        newick = self.traverse_newick_with_id(newick,labeled_sp)
        newick = f"{newick};"
        return newick

    def traverse_newick_with_id(self,newick,tree):
        if tree.leftChild and not tree.rightChild:
            newick = f"(,{self.traverse_newick_with_id(newick,tree.leftChild)}){self.write_events_with_id(tree)}"
        elif not tree.leftChild and tree.rightChild:
            newick = f"({self.traverse_newick_with_id(newick,tree.rightChild)},){self.write_events_with_id(tree)}"
        elif tree.leftChild and tree.rightChild:
            newick = f"({self.traverse_newick_with_id(newick,tree.rightChild)},{self.traverse_newick_with_id(newick,tree.leftChild)}){self.write_events_with_id(tree)}"
        elif not tree.leftChild and not tree.rightChild :
            newick = f"{self.write_events_with_id(tree)}"
        else:
            pass
        #newick=self.check_internal(newick,tree)
        return newick
    
    def write_events_with_id(self, tree):
        return f"{tree.taxa}[&label={tree.id_map_tag}]"
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
        
reco =reconcils()
red= readWrite.readWrite()
Il=ILS.ILS()

class daft_essential:
    def __init__(self):
        #self.reco =reconcils()
        #self.red= readWrite.readWrite()
        #self.Il=ILS.ILS()
        np.random.seed(42)


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
            print(Z_1,Z_2)
            Z_1 = 0.0 if pd.isna(Z_1) else Z_1
            Z_2 = 0.0 if pd.isna(Z_2) else Z_2
            if Z_1>=-4:
                if Z_2<Z_1:
                    pairs1 =pairs[-1]
                    pairs+=[[pairs1[1],pairs1[0]]]


        return pairs



    def write_events_network(self,tree):
            if tree.isLeaf:
                ev=''
                if len(tree.donor)==0 and len(tree.receiver)==0:
                    return ev+tree.taxa
                

                for idx,event in enumerate(tree.donor):
                    if idx==0:
                        ev= '('+tree.taxa+')'+event
                    else:
                        ev='('+ev+')'+event
                for idx ,event  in  enumerate(tree.receiver):
                    if idx==0 and len(tree.donor)==0:
                        ev= '('+event+','+tree.taxa+')'
                    else:
                        print('==========>')
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



    def put_network_in_tree(self,df,tree):

        pairs= self.find_pairs(df)


        sp = red.parse(tree)
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
            Z_score= (minZ_i-maxZ_i)/ math.sqrt(L1_count_i+L2_count_i)
            #if Z_score<-0.4:
            #continue
            df.at[i,'Z_score'] = Z_score
            
            
            found_j = None
            found_k = None
            for j in idxs:
                if j == i:
                    continue

                L1_j = df.at[j, 'Lineage1']
                L2_j = df.at[j, 'Lineage2']
                L1_count_j = df.at[j, 'Count1']
                L2_count_j = df.at[j, 'Count2']

                if self.isequal_set(L1_i, L1_j) and self._is_sibling(L2_i, L2_j,sp_string) and L2_count_i<L1_count_i:
                    found_j = j
                    found_k = 1
                    break

                if self.isequal_set(L1_i, L2_j) and self._is_sibling(L2_i, L1_j,sp_string) and L2_count_i<L1_count_i:
                    found_j = j
                    found_k = 2
                    break

                if self.isequal_set(L2_i, L1_j) and self._is_sibling(L1_i, L2_j,sp_string) and L1_count_i<L2_count_i:
                    found_j = j
                    found_k = 1
                    break

                if self.isequal_set(L2_i, L2_j) and self._is_sibling(L1_i, L1_j,sp_string) and L1_count_i<L2_count_i:
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
                    Z_score_sib= (L2_count_j-minZ_i)/ math.sqrt(minZ_i+L2_count_j)

                else:
                    df.at[i,'minor_sibling'] = df.at[found_j, 'Lineage1']
                    df.at[i,'minor_sibling_count'] = df.at[found_j, 'Count1']
                    Z_score_sib= (L1_count_j-minZ_i)/ math.sqrt(minZ_i+L1_count_j)

                df.at[i,'Z_score_sibling']=Z_score_sib

                if Z_score_sib<Z_score:
                    df.at[i,'CouldbeBidirectional']=True

        return df

    # ---------- styling----------
    def render_val(self,col, v):
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
        
    
    def render_val(self,col, v):
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

    def map_pair(self,what_moved,Elist):
        for pair in Elist:

            key=pair[0]
            uncle =pair[1]
            
        


            if self.isequal_set(what_moved,key):
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
        if col in ('Z-value-uncle', 'Z-value-sibling'):
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


    def _sibling_count(self,row,grouped):
        mask1 = grouped['What_moved'].apply(lambda x: self.isequal_set(x, row['comparison_sibling'])) 
        mask2 = grouped['Where_at'].apply(lambda x: self.isequal_set(x, row['Where_at']))
        
        if mask2.any() and not mask1.any():
            return pd.NA          
        return grouped.loc[mask1&mask2, 'total_count'].sum()

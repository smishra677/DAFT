import ast

class ILS:
    def find_bipartitions(self,bipartitions, subtree):
        if not subtree:
            return []

        stack = [subtree]

        while stack:
            node = stack.pop()

            if not node.isLeaf:
                bipartitions.append(node.taxa_list)
                
                if node.rightChild:
                    stack.append(node.rightChild)
                if node.leftChild:
                    stack.append(node.leftChild)

        return bipartitions

    def find_root(self,sp):
        stack = [sp] 
        parent=None
        while stack:
            current_node = stack.pop()  
            parent=current_node
            if current_node:
                if current_node.parent:
                    stack.append(current_node.parent)
            
        
        return parent

    def find_biparition_cost(self,sp,tr):

        
        merged = sp+tr



        difference= 0

        diff_parition=[]
        for i in merged:
            if i in tr and i not in sp:
                diff_parition.append(i)
                difference= difference+1
    
        
        return difference

    def get_current_sister(self,sp,Lineage):
            stack = [sp] 
                
            while stack:
                current_node = stack.pop()
                if current_node:
                    ##print(current_node.taxa, Lineage)
                    if current_node.taxa== Lineage:
                        if current_node.parent:
                            return [chil for chil in (current_node.parent.children) if chil!=current_node][0]

                
                    stack.append(current_node.leftChild)
                    stack.append(current_node.rightChild)


    def find_lineage(self,sp,Lineage):
            stack = [sp] 
                
            while stack:
                current_node = stack.pop()
                if current_node:
                    if current_node.taxa== Lineage:
                        #if current_node.parent:
                        return [current_node.parent][0]

                
                    stack.append(current_node.leftChild)
                    stack.append(current_node.rightChild)

    def get_nni_distance(self,tr,sp,Lineage):
            sp.label_internal()
            ##print('============xxxxxxxxxxxxx',Lineage)
            sister_=self.get_current_sister(sp,Lineage)
           
            
            if sister_:
                sister_.label_internal()
                ##print(sister_.taxa,type(sister_.taxa),type(sister_.taxa)!=set)
                if type(sister_.taxa)!=set:
                    sister= {str(sister_.taxa)}
                else:
                    sister= sister_.taxa
            else:
                return None
            ##print(sister)
            address =self.find_lineage(tr,Lineage)
            hops=0
            stack = [(address,hops)] 
            visited=[]   
            while stack:
                current_node,hops = stack.pop()
                if current_node:
                    current_node.label_internal()
                    ##print('434',sister,current_node.taxa)
                    ##print(current_node.taxa,sister,current_node.taxa== set(sister),hops,sister <= current_node.taxa)
                    if sister <= set(current_node.taxa):
                        if sister==current_node.taxa:
                            ##print('hops',hops)
                            return hops-1
                            ##print('Found',hops)
                        else:
                            for chil in (current_node.children):
                                if chil.isLeaf:
                                    if sister == set(chil.taxa):
                                        ##print('hops',hops)
                                        return hops
                                        #stack.append((chil,hops+1))
                                else:
                                    if sister <= chil.taxa:

                                        if chil not in visited:
                                            stack.append((chil,hops+1))
                                            visited.append(chil)
                                
                        
                    elif current_node.parent:
                        ##print(2)
                        ##print(current_node.parent.taxa)
                        if sister in current_node.parent.taxa:
                            ##print([chil for chil in (current_node.parent.children) if chil!=current_node][0].to_newick())
                            next_chil=[chil for chil in (current_node.parent.children) if chil!=current_node][0]
                            if next_chil not in visited:
                                stack.append((next_chil,hops+1))
                                visited.append(next_chil)
                            #stack.append(([chil for chil in (current_node.parent.children) if chil!=current_node][0],hops+1))
                            #visited.append([chil for chil in (current_node.parent.children) if chil!=current_node][0])
                        else:
                            ##print(3333331)
                            stack.append((current_node.parent,hops+1))
                    elif current_node.parent==None:
                        ##print([chil for chil in (current_node.children) if chil!=current_node][0].to_newick())
                        ##print(1)
                        next_chil=[chil for chil in (current_node.children) if chil!=current_node][0]
                        if next_chil not in visited:
                            stack.append((next_chil,hops+1))
                            visited.append(next_chil)
            ##print('hops',hops)
            ##print(sister,visited[-1].taxa.intersection(sister))
            #exit()
            return hops
    

    def find_number_map(self,sp):
        stack=[sp]
        total_map=0
        while stack:
            current_node=stack.pop()

            if current_node: 
                map_=len(current_node.refTo)
                if map_>1:
                    total_map+=map_

                stack.append(current_node.leftChild)

                stack.append(current_node.rightChild)
        
        return total_map


    def find_lowest(self,gene_tree,id_list,child_dic):
        

        level_dic={}

        address_list= {child_dic[id][1]:id for id in id_list}

        ##print(id_list)
        ##print(address_list)
        for add in address_list.keys():
            stack=[(add,0)]
            while stack:
                current_node,level=stack.pop()
                
                if current_node:
                    if current_node.isLeaf:
                        
                        level_dic[address_list[add]]=level

                    stack.append((current_node.leftChild,level+1))
                    stack.append((current_node.rightChild,level+1))
                
        id_max_level= min(level_dic, key=level_dic.get)
        ##print(level_dic)
        ##print([i.to_newick() for i in address_list.keys()])

        return id_max_level


        
        


    def pick_first_edge(self,child,gene_tree,tr,visited):

        if len(child)==0:
                return 
        else:
                pool={}
                tre_pool={}
                orientation={}
                super_list={}
                LCA_dic={}
                for k in range(len(child)):
                        
                        ch1=child[k]
                        
                        gene_tree.reset()

                        ch = ch1[0].deepcopy()
                        ch.reset()


                        list_tree= ch.NNI(gene_tree,ch1[2])

                        super_list[k]=list_tree
                        for li in list_tree:
                            
                               
                                
                                tr.reset()
                                li[0].reset()
                                li[0].order_gene(tr)

                                li[0].label_internal()
                                tr.label_internal()
                                bi_sp= self.find_bipartitions([],tr)
                                bi_li_0_l= self.find_bipartitions([],li[0].leftChild)
                                bi_li_0_r= self.find_bipartitions([],li[0].rightChild)
                                bi_score_left= self.find_biparition_cost(bi_sp,bi_li_0_l)
                                bi_score_right= self.find_biparition_cost(bi_sp,bi_li_0_r)

                                tr.reset()
                                tr.cost=0
                                li[0].reset()

                                li[0].leftChild.order_gene(tr)
                                li[0].leftChild.label_internal()
                                tr.label_internal()
    
                                li[0].leftChild.map_gene(tr)
                                tr.find_loss_sp(tr)
                                loss_left = tr.cost
                                number_map_left=len(tr.refTo)
                                



                                tr.reset()
                                tr.cost=0
                                li[0].rightChild.order_gene(tr)
                                li[0].rightChild.label_internal()
                                tr.label_internal()

                                
                                li[0].rightChild.map_gene(tr)
                                tr.find_loss_sp(tr)
                                loss_right = tr.cost
                                loss_score=loss_left+loss_right

                                tr.reset()
                                li[1].order_gene(tr)

                                tr.label_internal()
                                li[1].label_internal()
                                li[1].map_gene(tr)

                                #number_map= len(tr.refTo)
                                number_map_right=len(tr.refTo)
                                tr.reset()

                                ##print('##############################')
                                ##print(li[1].to_newick())
                                ##print(tr.to_newick())
                                ##print('top_0',li[0].to_newick())
                                ##print('topo_1',li[1].to_newick())
                                ##print('sp',tr.to_newick())
                                ##print('left_loss_cost',loss_left)
                                
                                ##print('right_loss_cost',loss_right)
                                
                                ##print('left_bi_cost',bi_score_left)
                                
                                ##print('right_bi_cost',bi_score_right) 
                                #firstTree = ete3.Tree(tr.to_newick())
                                #secondTree = ete3.Tree(li[1].to_newick())
                                #rf, _, _, _, _, _, _ = firstTree.robinson_foulds(secondTree)

                                #number_map=1
                                ##print((loss_score+bi_score_left+bi_score_right)*number_map)              
                                if number_map_right==0:
                                    number_map_right=1
                                if number_map_left==0:
                                    number_map_left=1

                                if k in LCA_dic:
                                    LCA_dic[k] = min(LCA_dic[k],number_map_left+number_map_right)
                                else:
                                    LCA_dic[k] = number_map_left+number_map_right

                                combined_score =(loss_left+bi_score_left)+(loss_right+bi_score_right)
                                
                                

                                if k not in pool:
                                    pool[k] = combined_score
                                    tre_pool[k] = li[1]
                                    orientation[k] = li[2]
                                else:
                                    if pool[k] > combined_score:
                                        pool[k] = combined_score
                                        tre_pool[k] = li[1]
                                        orientation[k] = li[2]



                min_key = min(pool, key=pool.get)
                min_combined_score=min(pool.values())
                min_value_keys_count = sum(1 for value in pool.values() if value == min_combined_score)


                
                if min_value_keys_count>1:
                    keys_with_min_value = [key for key, value in pool.items() if value == min_combined_score]
                    extracted_dic = {key: LCA_dic[key] for key in keys_with_min_value}
                    min_key = min(extracted_dic, key=extracted_dic.get)
                    min_combined_score_1=min(extracted_dic.values())
                    min_value_keys_count_1 = sum(1 for value in extracted_dic.values() if value == min_combined_score_1)
                    keys_with_min_value_1 = [key for key, value in extracted_dic.items() if value == min_combined_score_1]
                    if min_value_keys_count_1>1:
                        min_key=self.find_lowest(gene_tree,keys_with_min_value_1,child)

                

                        

                ##print("==============================<><()))()()()()()",pool,LCA_dic)

                return child[min_key],tre_pool[min_key],pool[min_key],orientation[min_key],super_list[min_key]

    def find_parent_child(self,root,child):
        if len(root.refTo)>1:
                for tre in root.refTo:
                    for tree1 in root.refTo:
                        
                        if (tree1 in tre.children):
                            if tree1==tre.leftChild:
                                #####print('match_left')
                                child.append([tre,tree1,'Left'])
                            else:
                                child.append([tre,tree1,'Right'])

        return child
    

    def parent_child(self,root,child):
        if root:
            if root.isLeaf:
                return []
            else:
                child= self.find_parent_child(root,child)
        return child

    def ILS(self,introgression,local_round,tracker, gene_tree, orginal_gene_tree,tr, sp_copy, cost,best_cost, visited,test_dic):
        ##print(1,cost)
        #Initial_best_cost=best_cost
        child = self.parent_child(tr, [])
   
        same_round=False
        if len(child) == 0 or cost <= 0:
            return gene_tree, cost, -1, visited,introgression,tracker,test_dic

        geneTree = gene_tree.deepcopy()
        geneTree.reset()

        #if len(child) == 1:
        if 1==1:
            ##print(introgression)
            ##print('*******************************************')
            list_tree = child[0][0].NNI(geneTree, child[0][2])
            child[0][0].label_internal()
            if local_round!=0:
                match_found=False
                for keys in introgression:
                    if len(set(ast.literal_eval(keys)).intersection(child[0][0].taxa)) >0:
                        introgression[keys]+=[[introgression[keys][0],[child[0][0].to_newick(),child[0][1].to_newick()]]]
                        match_found=True
                
                if not match_found:
                    introgression[repr(list(sorted(child[0][0].taxa)))]=[[child[0][0].to_newick(),child[0][1].to_newick()]]
        else:
            chil, trei, cos, orientation, list_tree = self.pick_first_edge(child, gene_tree, tr, visited)
            
            ##print(chil[0].to_newick(),list_tree[1][0].to_newick(),list_tree[1][1].to_newick())
            ##print(chil[0].to_newick(),list_tree[0][0].to_newick(),list_tree[0][1].to_newick())
            
            ##print(chil[1].to_newick())
            #print('=======================><<<<>><><><><><><><><><><><><><><>')

            
            if local_round!=0:
                ##print(introgression)
                #print('########################################')
                match_found=False
                chil[0].label_internal()
                '''                for keys in introgression:
                    if len(set(ast.literal_eval(keys)).intersection(chil[0].taxa)) >0:
                        introgression[keys]+=[[introgression[keys][0],[chil[0].to_newick(),chil[1].to_newick()]]]
                        match_found=True
                
                if not match_found:
                    introgression[repr(list(sorted(chil[0].taxa)))]=[[chil[0].to_newick(),chil[1].to_newick()]]
                '''

            if cos == 0:
                trei.label_internal()
                chii = self.get_child_info(chil, orientation)
                visited.append(chii)
                return trei, cost - 1, cos, visited,introgression,tracker,test_dic
            elif cos == -1:
                return gene_tree, 0, -1, visited,introgression,tracker,test_dic
            else:
                child = [chil]

        ch = child[0][1].deepcopy()
        ch.reset()

        #best_cost = cost
        improvement = False
        ''' 
        

        new_topo = geneTree.deepcopy()
               #tracker.append([child[0][0].taxa,child[0][1].taxa,local_round])
        child[0][0].label_internal()
        gene_tree.label_internal()
        other_child= [chil2 for chil2 in child[0][0].children if chil2!=child[0][1]][0]
        ne_t= self.get_nni_distance(self.find_root(orginal_gene_tree),sp_copy,other_child.taxa)
        ##print(self.find_root(gene_tree).to_newick(),child[0][1].rightChild.taxa,child[0][1].leftChild.taxa)
        ne_right=self.get_nni_distance(self.find_root(orginal_gene_tree),sp_copy,child[0][1].rightChild.taxa)
        ne_left=self.get_nni_distance(self.find_root(orginal_gene_tree),sp_copy,child[0][1].leftChild.taxa)
        

        #print('************************************',child[0][1].rightChild.taxa,ne_right,child[0][1].leftChild.taxa,ne_left)
        if ne_left==None:
            ne_left=-1
        
        if ne_right==None:
            ne_right=-1
        
        if ne_t==None:
            ne_t=-1
        if child[0][1].leftChild.to_newick() not  in test_dic:
            test_dic[child[0][1].leftChild.to_newick()]=ne_left
        
        if other_child.to_newick() not in test_dic:
            test_dic[other_child.to_newick()]=ne_t

        if child[0][1].rightChild.to_newick() not in test_dic:    
            test_dic[child[0][1].rightChild.to_newick()]=ne_right

     
        if child[0][1].leftChild.to_newick() in test_dic:
            
            if test_dic[child[0][1].leftChild.to_newick()]<ne_left:
                test_dic[child[0][1].leftChild.to_newick()]=ne_left
        
        else:
            test_dic[child[0][1].leftChild.to_newick()]=ne_left
        
        if other_child.to_newick() in test_dic:
            
            if test_dic[other_child.to_newick()]<ne_t:
                test_dic[other_child.to_newick()]=ne_t
        
        else:
            test_dic[other_child.to_newick()]=ne_t
        

        

        if child[0][1].rightChild.to_newick() in test_dic:
            

            
            if test_dic[child[0][1].rightChild.to_newick()]<ne_right:
                test_dic[child[0][1].rightChild.to_newick()]=ne_right
        else:
            test_dic[child[0][1].rightChild.to_newick()]=ne_right
        '''
        if child[0][2]=='Left':
            child_1=child[0][0].rightChild
            child_2=child[0][1].leftChild
            child_3=child[0][1].rightChild
        else:
            child_1=child[0][0].leftChild
            child_2=child[0][1].leftChild
            child_3=child[0][1].rightChild




        key_pos=repr(list(sorted(child[0][0].taxa)))
        if key_pos in tracker:
            tracker[key_pos]+=[[child_1.to_newick(),child_2.to_newick(),child_3.to_newick(),child[0][0].to_newick(),local_round]]  
        else:
            tracker[key_pos]=[[child_1.to_newick(),child_2.to_newick(),child_3.to_newick(),child[0][0].to_newick(),local_round]]  

        for i in list_tree:

            i[1].reset()
            i[0].reset()
            cop = sp_copy.deepcopy()
            cop.reset()

            i[1].order_gene(cop)
            i[1].label_internal()
            cop.label_internal()
            i[1].map_gene(cop)
            new_cost = len(cop.refTo)

            if best_cost >= new_cost and cost >= 0:
                improvement = best_cost >= new_cost
                best_cost = new_cost
                new_topo = i[1].deepcopy()

                chii = self.get_child_info(child[0], i[2])
                
                

 
                if same_round:
                    visited=visited[:-1]
                else:
                    same_round=True
                    cost -=1
                visited.append(chii)


        ##print('-------------------------------------------cxsdsdsdsdsdsd',new_topo.to_newick())
                
        #tracker.append([new_topo.taxa,new_topo.to_newick(),local_round])

 
        if cost <= 0 or not improvement:
            return new_topo, cost, -1, visited,introgression,tracker,test_dic

        new_sp = sp_copy.deepcopy()
        new_sp.reset()
        new_topo.reset()
        new_topo.order_gene(new_sp)
        new_topo.label_internal()
        new_sp.label_internal()
        new_topo.map_gene(new_sp)
        local_round=local_round+1

        return self.ILS(introgression,local_round,tracker,new_topo,orginal_gene_tree, new_sp, sp_copy, cost,best_cost, visited,test_dic)

    def get_child_info(self, chil, orientation):
        if chil[2] == 'Left':
            if orientation == 'left':
                return [chil[0].leftChild, chil[0].leftChild.leftChild]
            else:
                return [chil[0].leftChild, chil[0].leftChild.rightChild]
        else:
            if orientation == 'left':
                return [chil[0].rightChild, chil[0].rightChild.leftChild]
            else:
                return [chil[0].rightChild, chil[0].rightChild.rightChild]

#calculate single-cut-and-join-score for sampled tree 
def calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj):
    scj_dist=0 #initialize the single-cut-or-join-distance

    #traversing all nodes except the root
    for node in tree.iter_descendants("postorder"):
        #print node.name+"_"+node.up.name #this pair represent the current considered edge of the tree
        #the adjacencies of each node are either in reconstructedAdjs (internal nodes and root) or in extantAdjacenies (leaves)
        #for each node, get the set of adjs
        set_of_node_adj=set()
        if node.name in reconstructedAdj: # if the current node is an internal node ...
            set_of_node_adj=set(reconstructedAdj[node.name]) #... get the adjacencies from reconstructedAdjs
        elif node.name in extantAdjacencies_species_adj: # if the current node isn't an internal node
            set_of_node_adj=set(extantAdjacencies_species_adj[node.name]) #... get the adjacencies from extantAdjacencies_species_adj
        else:
            print "ERROR: The species' name wasn't found"
        #for each parent node of the current node get the set of adjs
        set_of_nodes_parent_adj=set()
        if node.up.name in reconstructedAdj: # if the parent of the current node is an internal node...
            set_of_nodes_parent_adj=set(reconstructedAdj[node.up.name])	 #... get the adjacencies from reconstructedAdjs
        #elif node.up.name in extantAdjacencies_species_adj:# if the current node's parent isn't an internal node (shouldn't occur)
        #     set_of_nodes_parent_adj=set(extantAdjacencies_species_adj[node.up.name]) #... get the adjacencies from extantAdjacencies_species_adj
        else:
             print "ERROR: The species' name wasn't found"
        #print set_of_nodes_parent_adj

        #build the SCJ distance from the parent node's set to the nodes set
        scj_node=len(set_of_node_adj)
        for adj in set_of_node_adj:
            mirrored_adj=(adj[1],adj[0])
           # print str(adj)+'__'+str(mirrored_adj)
            if adj in set_of_nodes_parent_adj or mirrored_adj in set_of_nodes_parent_adj:
                scj_node=scj_node-1

        #build the SCJ distance from the node's set to the parent node's set
        scj_parent=len(set_of_nodes_parent_adj)
        for adj in set_of_nodes_parent_adj:
            mirrored_adj=(adj[1],adj[0])
            #print str(adj)+'__'+str(mirrored_adj)
            if adj in set_of_node_adj or mirrored_adj in set_of_node_adj:
                scj_parent=scj_parent-1
        #calculate the SCJ distance of the current pair and add it to the total scj-distance of the whole tree
        temp=abs(scj_node)+abs(scj_parent)
        #temp=len(set_of_node_adj.difference(set_of_nodes_parent_adj))+len(set_of_nodes_parent_adj.difference(set_of_node_adj))
        print temp
        scj_dist+= temp
    print "Single-Cut-or-Join-Distance: "+str(scj_dist)
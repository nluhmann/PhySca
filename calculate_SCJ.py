

#calculate single-cut-and-join-score for sampled tree
def calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj):
    scj_dist=0 #initialize the single-cut-or-join-distance


    for node in tree.iter_descendants("postorder"):
        #the adjacencies of each node are either in reconstructedAdjs (internal nodes and root) or in extantAdjacenies (leaves)
        set_of_node_adj=set()
        if node.name in reconstructedAdj: # if the current node is an internal node ...
            set_of_node_adj=set(reconstructedAdj[node.name])
        elif node.name in extantAdjacencies_species_adj:
            set_of_node_adj=set(extantAdjacencies_species_adj[node.name])
        else:
            print "ERROR: The species' name wasn't found"

        #for each parent node of the current node get the set of adjacencies
        set_of_nodes_parent_adj=set()
        if node.up.name in reconstructedAdj:
            set_of_nodes_parent_adj=set(reconstructedAdj[node.up.name])
        else:
             print "ERROR: The species' name wasn't found"

        #compute the SCJ distance from the parent node's set to the nodes set
        scj_node=len(set_of_node_adj) #scj_node counts how many elements are only at the node
        for adj in set_of_node_adj:
            mirrored_adj=(adj[1],adj[0])
           # print str(adj)+'__'+str(mirrored_adj)
            if adj in set_of_nodes_parent_adj or mirrored_adj in set_of_nodes_parent_adj:
                scj_node=scj_node-1 #if the current adj is also at the node's parent (even if it's mirrored), deduct 1 from scj_node

        #compute the SCJ distance from the node's set to the parent node's set
        scj_parent=len(set_of_nodes_parent_adj)
        for adj in set_of_nodes_parent_adj:
            mirrored_adj=(adj[1],adj[0])
            if adj in set_of_node_adj or mirrored_adj in set_of_node_adj:
                scj_parent=scj_parent-1
        #calculate the SCJ distance of the current pair and add it to the total scj-distance of the whole tree
        temp=abs(scj_node)+abs(scj_parent)
        scj_dist+= temp
    return scj_dist
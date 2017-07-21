__author__ = 'nluhmann'
import networkx as nx
import random


def enumJointLabelings(ccs):
    print 'Enumerate joint labelings...'

    joint = {}
    fi = {}
    for cc in ccs:
        labelings = set()
        edges = cc.edges()
        nodes = cc.nodes()
        #add first label to list of joint labels, all extremities unassigned
        first = constructLabel(nodes,[])
        first.sort(key=lambda tup: tup[0])
        fi[cc] = first
        labelings.add(tuple(first))
        counter = len(edges)
        #add edges to labeling if it is valid
        for i in xrange(0,counter):
            newlabels = set()
            for labeling in labelings:
                for edg in edges:
                    #edg = sortEdge(edg)
                    if (edg[0],"-") in labeling and (edg[1],"-") in labeling:
                        newlabel = list(labeling)
                        newlabel.remove((edg[0], "-"))
                        newlabel.remove((edg[1], "-"))
                        newlabel.append((edg[0], edg))
                        newlabel.append((edg[1], edg))
                        newlabel.sort(key=lambda tup: tup[0])
                        newlabels.add(tuple(newlabel))
            labelings = labelings.union(newlabels)
        joint[cc] = labelings

    return joint,fi

def sortEdge(edge):
    first = edge[0]
    second = edge[1]
    if first < second:
        new = (first,second)
    else:
        new = (second,first)
    return new

def constructLabel(allnodes,edges):

    #edge is a tuple of two nodes, label these nodes with the edge
    label = []
    for node in allnodes:
        if len(edges) == 2:
            if node in edges:
                label.append((node,edges))
            else:
                label.append((node,"-"))
        else:
            labeled = False
            for edge in edges:
                if node in edge:
                    edge = sortEdge(edge)
                    label.append((node,edge))
                    labeled = True
            if not labeled:
                label.append((node,"-"))
    return label

def validLabels(joint,first):
    print "Check valid labels..."
    #check for which internal nodes of the tree, a label is valid (all edges in the label are marked with the internal node)
    #remove all labels that are not valid at any internal node of the tree

    valid = {}
    validAtNode = {}
    for cc in joint:

        valid[cc] = [first[cc]]
        validAtNode[cc] = {}
        validAtNode[cc]["all"] = tuple(first[cc])
        attr = nx.get_edge_attributes(cc,"species")
        labels = joint[cc]
        for label in labels:
            species = set()
            listLabel = list(label)
            for elem in listLabel:
                if not "-" in elem:
                    edge = elem[1]
                    spec = attr[edge]
                    if not species:
                        for n in spec:
                            species.add(n)
                    else:
                        species = species.intersection(spec)
                        if not species:
                            break
            if species:
                valid[cc].append((label,species))
                for sp in species:
                    try:
                        name = sp.name
                    except AttributeError:
                        name = sp

                    if name in validAtNode[cc]:
                        validAtNode[cc][name].add(label)
                    else:
                        s = set()
                        s.add(label)
                        validAtNode[cc][name] = s

    return valid, validAtNode

def useAllLabels(joint,first,tree):
    validAtNode = {}

    for cc in joint:
        validAtNode[cc] = {}
        validAtNode[cc]["all"] = tuple(first[cc])
        labels = joint[cc]
        for node in tree.traverse():
            for label in labels:
                if node.name in validAtNode[cc]:
                        validAtNode[cc][node.name].add(label)
                else:
                    s = set()
                    s.add(label)
                    validAtNode[cc][node.name] = s
    return validAtNode

def cost(parentLabel, childLabel, edges, probs, node, alpha):

    #label is a list of tuple
    parentAdj = set()
    if len(parentLabel) == 1:
        if not "-" in parentLabel:
            parentAdj.add(parentLabel[1])
    else:
        for elem in parentLabel:
            if not "-" in elem:
                parentAdj.add(elem[1])
    childAdj = set()
    if len(childLabel) == 1:
        if not "-" in childLabel:
            childAdj.add(childLabel[1])
    else:
        for elem in childLabel:
            if not "-" in elem:
                childAdj.add(elem[1])

    cost = 0
    #TODO edge Length
    edgeLength = 1
    for edge in edges:
        #edge is present in parent but lost in child
        if edge in parentAdj:
            if not edge in childAdj:
                if node.name in probs:
                    if edge in probs[node.name]:
                        weight = float(probs[node.name][edge])
                    elif (edge[1],edge[0]) in probs[node.name]:
                        weight = float(probs[node.name][(edge[1],edge[0])])
                    else:
                        weight = 0
                else:
                    weight = 0

                cost = cost + (alpha*weight) + ((1-alpha)/edgeLength)
        #edge is gained in child
        elif edge in childAdj:
            if not edge in parentAdj:
                cost += (1 - alpha) / edgeLength
        #edge is absent both in parent in child
        else:
            if node.name in probs:
                if edge in probs[node.name]:
                    weight = float(probs[node.name][edge])
                elif (edge[1],edge[0]) in probs[node.name]:
                    weight = float(probs[node.name][(edge[1],edge[0])])
                else:
                    weight = 0
            else:
                weight = 0
            cost += alpha * weight

    return cost

# SE: here, one leaf is annotated with one label, according to presence or absence of adjacency in the extant adjacencies
def annotateleaves(leaf,cc,extant, potentialExtant):

    nodes = cc.nodes()
    edges = cc.edges()
    label = []

    for node in nodes:
        added = False
        for elem in extant:
            if node in elem:
                for spec in extant[elem]:
                    if leaf.name in spec:
                        if not elem in edges:
                            elem = (elem[1],elem[0])
                        if elem in edges:
                            label.append((node, elem))
                            added = True
        for elem in potentialExtant:
            if node in elem:
                for spec in potentialExtant[elem]:
                    if leaf.name in spec:
                        if not elem in edges:
                            elem = (elem[1],elem[0])
                        if elem in edges:
                            label.append((node, elem))
                            added = True

        if not added:
            label.append((node, "-"))

    return label

def sankoff_bottomup(t,space, edges, probs, alpha):


    for node in t.traverse(strategy="postorder"):
        #get annotation for all children of node
        if not node.is_leaf():

            annoHash = {}

            #in case the empty label is the only option at an internal node
            if not node.name in space:
                space[node.name] = [space["all"]]

            #minimum cost for this label over all children of this node
            allMin = float("inf")

            minLabel = ""


            for label in space[node.name]:
                minTotal = 0

                #for all children of the current node
                for child in node.children:
                    #minimum cost over all possible labels in this child
                    minChild = float("inf")
                    #annotation for all possible labels in this child
                    anno = child.annotation
                    #for each annotation
                    for annotation in anno:
                        #compute cost from label to annotation!
                        value = anno[annotation] + cost(label,annotation,edges,probs,node,alpha)
                        #if the cost is smaller than before, update minChild value
                        if value < minChild:
                            minChild = value
                    #add minimum cost for this child to the total cost
                    minTotal += minChild
                #save minimum cost for this label over all children
                annoHash[label] = minTotal
                #if current node is the root, we need to add the cost for this label in terms of weight at the root!
                if node.is_root():
                    minTotal += computeLabelWeight(label, edges, alpha, probs, node)
                #if minimum cost for this label is smaller than before, set this as smallest label
                if minTotal < allMin:
                    minLabel = label
                    allMin = minTotal
            node.add_feature("annotation",annoHash)
            node.add_feature("minimumLabel",minLabel)
    return t

def sankoff_topdown(t,edges, probs, alpha):
    for node in t.traverse(strategy="preorder"):
        if node.is_root():
            final = node.minimumLabel
            node.add_feature("assignment",final)
        else:
            parentAssignment = node.up.assignment
            annotation = node.annotation
            minAllCost = float("inf")
            minAllLabel = ""
            for anno in annotation:
                minCost = annotation[anno] + cost(parentAssignment,anno,edges,probs,node,alpha)
                if minCost < minAllCost:
                    minAllCost = minCost
                    minAllLabel = anno
            node.add_feature("assignment",minAllLabel)
    return t

def computeLabelings(tree, ccs, validAtNode, extant, probabilities, alpha, ancientLeaves, potentialExtant):
    print "Compute ancestral labels with SR..."

    adjacencies = {}
    for cc in ccs:

        treecopy = tree.copy()
        #1) annotate each tree node except root with cost matrix

        for node in treecopy.traverse(strategy="postorder"):
            #annotate leafs by their only possible label (maximum number of adjacencies)

            if node.is_leaf():
                #annotateleaves returns a hash: key=label, value=cost (which is zero for all leaves)
                #if the leaves is extant, the hash will contain only one entry
                #otherwise it will contain all potential labels at an ancient leaf
                if not node.name in ancientLeaves:
                    lab = annotateleaves(node,cc,extant,potentialExtant)
                    node.add_feature("annotation",{tuple(lab):0})

                    #if the leaf is extant, we can already set the minimumLabel here
                    #otherwise, we leaf this open
                    node.add_feature("minimumLabel",lab)

                    #add valid labels
                    label = tuple(lab)
                    validAtNode[cc][node.name] = set()
                    validAtNode[cc][node.name].add(label)
                else:
                    #we have an ancient Leaf here! get all labels that are valid at this node for this cc
                    valid_labels = validAtNode[cc][node.name]
                    #for each validLabel, its cost are the sum of the weights of all absent adjacencies
                    anno_hash = {}
                    adjEdges = cc.edges()
                    for lab in valid_labels:
                        if not alpha == 0:
                            cost = 0
                            for adj in adjEdges:
                                inside = False
                                for n in lab:
                                    if adj in n:
                                        inside = True
                                        continue
                                if not inside:
                                    cost = cost + probabilities[node.name][sortEdge(adj)]
                            anno_hash.update({tuple(lab):int(cost)*alpha})
                        else:
                            anno_hash.update({tuple(lab):0})
                    node.add_feature("annotation",anno_hash)


            if not node.is_root():
                parent = node.up
                #find all valid labels at node and at parent
                if node.name in validAtNode[cc]:
                    labelsNode = validAtNode[cc][node.name]
                    labelsNode.add(tuple(validAtNode[cc]["all"]))
                else:
                    labelsNode = set()
                    labelsNode.add(tuple(validAtNode[cc]["all"]))
                if parent.name in validAtNode[cc]:
                    labelsParent = validAtNode[cc][parent.name]
                    labelsParent.add(tuple(validAtNode[cc]["all"]))
                else:
                    labelsParent = set()
                    labelsParent.add(tuple(validAtNode[cc]["all"]))

            elif node.is_root():
                #find all valid labels at root
                if node.name in validAtNode[cc]:
                    labelsNode = validAtNode[cc][node.name]
                    labelsNode.add(tuple(validAtNode[cc]["all"]))
                else:
                    labelsNode = set()
                    labelsNode.add(tuple(validAtNode[cc]["all"]))



        #then compute bottom-up labeling for tree
        annoTree = sankoff_bottomup(treecopy,validAtNode[cc], cc.edges(), probabilities, alpha)

        #then compute top-down labeling for tree
        annotatedTree = sankoff_topdown(annoTree, cc.edges(), probabilities, alpha)

        for node in annotatedTree.traverse():
            if not node.is_leaf() or node.name in ancientLeaves:
                label = node.assignment
                labelAdj = set()
                if len(label) == 1:
                    if not "-" in label:
                        labelAdj.add(label[1])
                else:
                    for elem in label:
                        if not "-" in elem:
                            labelAdj.add(elem[1])
                if labelAdj:
                    if node.name in adjacencies:
                        adjacencies[node.name] += list(labelAdj)
                    else:
                        adjacencies[node.name] = list(labelAdj)

    return adjacencies

def computeLabelWeight(label,edges,alpha, probs,node):
    adj = set()
    if len(label) == 1:
        if not "-" in label:
            adj.add(label[1])
    else:
        for elem in label:
            if not "-" in elem:
                adj.add(elem[1])
    cost = 0

    for edge in edges:
        if not edge in adj:
            if node.name in probs:
                if edge in probs[node.name]:
                        weight = float(probs[node.name][edge])
                elif (edge[1],edge[0]) in probs[node.name]:
                        weight = float(probs[node.name][(edge[1],edge[0])])
                else:
                        weight = 0
            else:
                weight = 0
            cost += alpha * weight

    return cost

def reconstructedAdjacencies(resolvedCCs):
    adjacencies = {}

    for cc in resolvedCCs:
        annotatedTree = resolvedCCs[cc]
        for node in annotatedTree.traverse():
            if not node.is_leaf():
                label = node.assignment
                labelAdj = set()
                if len(label) == 1:
                    if not "-" in label:
                        labelAdj.add(label[1])
                else:
                    for elem in label:
                        if not "-" in elem:
                            labelAdj.add(elem[1])
                if labelAdj:
                    if node.name in adjacencies:
                        adjacencies[node.name] += list(labelAdj)
                    else:
                        adjacencies[node.name] = list(labelAdj)

    return adjacencies

def outputReconstructedAdjacencies(reconstrAdj,out):
    file = open(out,"w")
    for node in reconstrAdj:
        file.write(">"+node+" "+str(len(reconstrAdj[node]))+"\n")
        for adj in reconstrAdj[node]:
            file.write(node+" "+str(adj)+"\n")
    file.close()

def sankoff_bottomup_sampling(t,space, edges, probs, alpha):
    for node in t.traverse(strategy="postorder"):
        #get annotation for all children of node
        if not node.is_leaf():
            annoHash = {}
            optimalsHash = {}
            #in case the empty label is the only option at an internal node
            if not node.name in space:
                space[node.name] = [space["all"]]

            #minimum cost for this label over all children of this node
            allMin = float("inf")
            #labels having the minimum cost at this node
            minLabel = []


            #for each possible label at a node
            for label in space[node.name]:
                minTotal = 0
                optSolTotal = 1

                #for all children of the current node
                for child in node.children:
                    #sum of optimal solutions in this child
                    optSolSum = 0
                    #minimum cost over all possible labels in this child
                    minChild = float("inf")
                    #labels that give minimum costs
                    childMinLabels = []

                    #annotation for all possible labels in this child
                    anno = child.annotation
                    #for each annotation
                    for annotation in anno:
                        #compute cost from label to annotation!
                        value = anno[annotation] + cost(label,annotation,edges,probs,node,alpha)
                        #if the cost is smaller than before, update minChild value
                        if value < minChild:
                            minChild = value
                            childMinLabels = [annotation]
                        elif value == minChild:
                            childMinLabels.append(annotation)
                    #add minimum cost for this child to the total cost
                    minTotal += minChild
                    #sum over optimal solutions in all annotations giving the minimum cost
                    for a in childMinLabels:
                        optSolSum = optSolSum + child.numberOfOptimals[a]
                    optSolTotal *= optSolSum
                #save minimum cost for this label over all children
                annoHash[label] = minTotal
                #if current node is the root, we need to add the cost for this label in terms of weight at the root!
                if node.is_root():
                    minTotal += computeLabelWeight(label, edges, alpha, probs, node)
                #if minimum cost for this label is smaller than before, set this as smallest label
                if minTotal < allMin:
                    minLabel = [label]
                    allMin = minTotal
                elif minTotal == allMin:
                    minLabel.append(label)
                optimalsHash[label] = optSolTotal
            node.add_feature("annotation",annoHash)
            node.add_feature("minimumLabel",minLabel)
            node.add_feature("numberOfOptimals",optimalsHash)
    return t

def sankoff_topdown_sampling(t, edges, probs, alpha):
    for node in t.traverse(strategy="preorder"):
        if node.is_root():
            #sum NumberOfOptimals for all optimal labels at root
            minLabels = node.minimumLabel
            sumOpt = 0


            for lab in minLabels:
                sumOpt = sumOpt + node.numberOfOptimals[lab]

            draw = random.randint(1, sumOpt)
            x = draw
            final = ""
            for lab in minLabels:
                x = x - node.numberOfOptimals[lab]
                if x <= 0:
                    final = lab
                    break

            node.add_feature("assignment",final)
        else:
            parentAssignment = node.up.assignment
            annotation = node.annotation
            minAllCost = float("inf")
            minAllLabel = []
            for anno in annotation:
                minCost = annotation[anno] + cost(parentAssignment,anno,edges,probs,node,alpha)
                if minCost < minAllCost:
                    minAllCost = minCost
                    minAllLabel = [anno]
                elif minCost == minAllCost:
                    minAllLabel.append(anno)
            #draw from all optimal annotations like at the root
            sumOpt = 0
            for lab in minAllLabel:
                sumOpt = sumOpt + node.numberOfOptimals[lab]

            draw = random.randint(1, sumOpt)
            x = draw
            final = ""
            for lab in minAllLabel:
                x = x - node.numberOfOptimals[lab]
                if x <= 0:
                    final = lab
                    break


            node.add_feature("assignment",final)
    return t

def sampleLabelings(tree, ccs, validAtNode, extant, probabilities, alpha, ancientLeaves, potentialExtant):
    print "Compute ancestral labels with SR..."

    adjacencies = {}
    for cc in ccs:
        treecopy = tree.copy()
        #1) annotate each tree node except root with cost matrix
        for node in treecopy.traverse(strategy="postorder"):
            #annotate leafs by their only possible label (maximum number of adjacencies)
            if node.is_leaf():
                if not node.name in ancientLeaves:
                    lab = annotateleaves(node,cc,extant, potentialExtant)
                    node.add_feature("annotation",{tuple(lab):0})
                    node.add_feature("minimumLabel",[lab])
                    node.add_feature("numberOfOptimals",{tuple(lab):1})
                    label = tuple(lab)
                    validAtNode[cc][node.name] = set()
                    validAtNode[cc][node.name].add(label)
                else:
                    #we have an ancient Leaf here! get all labels that are valid at this node for this cc
                    valid_labels = validAtNode[cc][node.name]
                    #for each validLabel, its cost are the sum of the weights of all absent adjacencies
                    anno_hash = {}
                    numberOptimals = {}
                    adjEdges = cc.edges()
                    for lab in valid_labels:
                        if not alpha == 0:
                            cost = 0
                            for adj in adjEdges:
                                inside = False
                                for n in lab:
                                    if adj in n:
                                        inside = True
                                        continue
                                if not inside:
                                    cost = cost + probabilities[node.name][sortEdge(adj)]
                            anno_hash.update({tuple(lab):int(cost)*alpha})

                        else:
                            anno_hash.update({tuple(lab):0})
                        numberOptimals.update({tuple(lab):1})
                    node.add_feature("annotation",anno_hash)
                    node.add_feature("numberOfOptimals",numberOptimals)

            if not node.is_root():
                parent = node.up
                #find all valid labels at node and at parent
                if node.name in validAtNode[cc]:
                    labelsNode = validAtNode[cc][node.name]
                    labelsNode.add(tuple(validAtNode[cc]["all"]))
                else:
                    labelsNode = set()
                    labelsNode.add(tuple(validAtNode[cc]["all"]))
                if parent.name in validAtNode[cc]:
                    labelsParent = validAtNode[cc][parent.name]
                    labelsParent.add(tuple(validAtNode[cc]["all"]))
                else:
                    labelsParent = set()
                    labelsParent.add(tuple(validAtNode[cc]["all"]))

            elif node.is_root():
                #find all valid labels at root
                if node.name in validAtNode[cc]:
                    labelsNode = validAtNode[cc][node.name]
                    labelsNode.add(tuple(validAtNode[cc]["all"]))
                else:
                    labelsNode = set()
                    labelsNode.add(tuple(validAtNode[cc]["all"]))

        #then compute bottom-up labeling for tree
        annoTree = sankoff_bottomup_sampling(treecopy,validAtNode[cc], cc.edges(), probabilities, alpha)
        #then compute top-down labeling for tree
        annotatedTree = sankoff_topdown_sampling(annoTree,cc.edges(),probabilities,alpha)
        for node in annotatedTree.traverse():
            if not node.is_leaf() or node.name in ancientLeaves:
                label = node.assignment
                labelAdj = set()
                if len(label) == 1:
                    if not "-" in label:
                        labelAdj.add(label[1])
                else:
                    for elem in label:
                        if not "-" in elem:
                            labelAdj.add(elem[1])
                if labelAdj:
                    if node.name in adjacencies:
                        adjacencies[node.name] += list(labelAdj)
                    else:
                        adjacencies[node.name] = list(labelAdj)

    return adjacencies
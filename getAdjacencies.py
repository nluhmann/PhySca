__author__ = 'Nina'


import itertools
import subprocess
import os




# double marker id to account for orientation
def doubleMarker(marker):
    if "-" in marker:
        return int(marker[1:]) * 2, (int(marker[1:]) * 2) - 1
    else:
        return (int(marker) * 2) - 1, int(marker) * 2


# double marker and find adjacencies in chromosomes
def findAdjacencies(speciesHash):
    print "collect extant adjacencies from marker file..."
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    for species in speciesHash:
        chromosomes = speciesHash[species]
        for chrom in chromosomes:
            orderList = chromosomes[chrom]
            for i in range(0,len(orderList)-1):
                first = orderList[i]
                second = orderList[i+1]
                doubleFirst = doubleMarker(first)
                doubleSecond = doubleMarker(second)
                #take right extrem of first and left extrem of second for adj tuple
                adj = (doubleFirst[1],doubleSecond[0])
                rev = (doubleSecond[0],doubleFirst[1])
                if adj in adjacencies:
                    adjacencies[adj].append((species,chrom))
                elif rev in adjacencies:
                    adjacencies[rev].append((species,chrom))
                else:
                    adjacencies[adj] = [(species,chrom)]
    return adjacencies



#read FPSAC adjacency file
def readAdjacencyFile(file):
    print "collect extant adjacencies from provided adjacency file"
    # keys: (left marker, right marker), value: [(species,chromosome),...]
    adjacencies = {}
    f = open(file,"r")
    for line in f:
        fields = line.rstrip("\n").split("#")
        adjas = fields[0].split(":")
        left = adjas[1].split(" ")[0]
        right = adjas[1].split(" ")[1]
        specieslist = []
        for extant in fields[1:]:
            species = extant.split(":")[0].split(".")[0]
            specieslist.append((species,"mockchrom"))
        specieslist = set(specieslist)
        adjacencies[(left,right)] = specieslist
    return adjacencies


def outputAdjacencies(adjacencies):
    adjOut = "extant_adjacencies"
    file = open(adjOut,"w")
    counter = 1
    for elem in adjacencies:
        file.write(">"+str(counter)+"\t")
        for spec in adjacencies[elem]:
            file.write("#"+str(spec[0])+","+str(spec[1])+" ")
        file.write(":"+str(elem[0])+"-"+str(elem[1])+"\n")
        counter = counter + 1
    file.close()


def findTreePaths(tree):
    pairs = {}
    leaves = tree.get_leaves()
    #get all pairs of leaves
    for pair in itertools.combinations(leaves,2):
        ancestor = pair[0].get_common_ancestor(pair[1])
        node = pair[0].up
        path = []
        while not node == ancestor:
            path.append(node)
            node = node.up
        path.append(ancestor)
        node = pair[1].up
        while not node == ancestor:
            path.append(node)
            node = node.up
        pairs[pair] = path
    return pairs


def assignAncestralAdjacencies(paths, adjacencies, tree):
    #initialize dict for all internal nodes
    internal = {}
    for node in tree.traverse():
         if not node.is_leaf():
             s = set()
             internal[node] = s

    adjacenciesAncestral = {}
    for adj in adjacencies:
        for pair in itertools.combinations(adjacencies[adj],2):
            one = tree&pair[0][0]
            two = tree&pair[1][0]
            if (one,two) in paths:
                tup = (one,two)
            elif (two,one) in paths:
                tup = (two,one)
            else:
                print "Houston there is a problem here!"
            for n in paths[tup]:
                internal[n].add(adj)
                if adj in adjacenciesAncestral:
                    adjacenciesAncestral[adj].add(n.name)
                else:
                    se = set()
                    se.add(n.name)
                    adjacenciesAncestral[adj] = se
    return internal,adjacenciesAncestral


def outputAncestralAdjacencies(internal,out):
    file = open(out,"w")
    for node in internal:
        file.write(">"+node.name+" "+str(len(internal[node]))+"\n")
        adjacencies = internal[node]
        for adj in adjacencies:
            file.write(str(adj)+"\n")
    file.close()


#for each extant adjacency, use declone to compute probability that the
#adjacency is present in an adjacency forest sampled randomly from a
#Boltzmann distribution. Then assign adjacency if it is above threshold to internal node of the tree.

def deCloneProbabilities(extantAdjacencies, threshold, tree, alpha, kT, treefile):
    print "Compute probabilities with DeClone..."
    nodesPerAdjacency = {}
    adjacenciesPerNode = {}
    adjacencyProbs = {}
    for adjacency in extantAdjacencies:
        #produce list of extant adjacencies for declone
        path = os.path.dirname(os.path.realpath(__file__))
        listVar = "list_"+str(threshold)+str(alpha)+str(kT)
        file = open(path+"/"+listVar,"w")
        species = extantAdjacencies[adjacency]
        for spec in species:
            file.write(spec[0]+" "+spec[0]+"\n")
            if spec in adjacencyProbs:
                adjacencyProbs[spec][adjacency] = 1
            else:
                adjacencyProbs[spec]={adjacency:1}
        file.close()
        #use declone to compute probs

        command = './DeClone -t1 '+treefile+' -t2 '+treefile+' -a '+path+"/"+listVar+' -i -kT '+str(kT)
        output = subprocess.check_output(command, shell=True, cwd=path)
        #output is just matrix with probabilities
        #each line of the output should contain max one number greater 0, save for internal nodes
        lines = output.split("\n")
        for line in lines:
            if not line == "" and not line.startswith("\t") and not line.startswith(">"):
                node = line.split("\t")[0]
                probs = line.split("\t")[1]
                probs = probs.rstrip(" ")
                probArray = probs.split(" ")
                probArrayFl = [float(x) for x in probArray]
                probability = max(probArrayFl, key=float)
                treeNode = tree.search_nodes(name=node)
                if not treeNode:
                    print "Error, no node found"
                elif not treeNode[0].is_leaf():
                    if probability > threshold and not len(species) == 1:
                        if adjacency in nodesPerAdjacency:
                            nodesPerAdjacency[adjacency].add(treeNode[0].name)
                        else:
                            nodeset = set()
                            nodeset.add(treeNode[0].name)
                            nodesPerAdjacency[adjacency] = nodeset
                        if treeNode[0].name in adjacenciesPerNode:
                            adjacenciesPerNode[treeNode[0].name].append((adjacency,probability))
                        else:
                            adjacenciesPerNode[treeNode[0].name] = [(adjacency,probability)]
                    if treeNode[0].name in adjacencyProbs:
                        adjacencyProbs[treeNode[0].name][adjacency] = probability
                    else:
                        adjacencyProbs[treeNode[0].name] = {adjacency:probability}
    return nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode

def outputAncestralAdjacenciesWithProbabilities(adjacenciesPerNode, out):
    file = open(out,"w")
    for node in adjacenciesPerNode:
        file.write(">"+node+"\t"+str(len(adjacenciesPerNode[node]))+"\n")
        adjacencies = adjacenciesPerNode[node]
        for adj in adjacencies:
            file.write(str(adj[0])+" "+str(adj[1])+"\n")
    file.close()


#node='BD'
def computeWeightProbs(weighting, extantAdjacencies, tree):
    adjacencyProbs = {}
    #read probabilities from file, assume that all assigned adjacencies at this node have a weight in this file
    for node in tree.traverse():
        adjacencyProbs[node.name] = {}

    for adjacency in extantAdjacencies:
        for key in adjacencyProbs:
            adjacencyProbs[key][adjacency] = 0

    file = open(weighting,"r")
    for line in file:
        if line.startswith("#"):
            n = line.rstrip("\n")[1:]
        else:
            fields = line.split("\t")
            left = fields[0][1:]
            right = fields[1]
            weight = fields[3]
            adjacencyProbs[n][(left,right)] = float(weight)
    file.close()

    return adjacencyProbs




def assignGAMLweights(weighting, nodesPerAdjacency2, adjacencyProbs2, adjacenciesPerNode2, threshold):
    #read probabilities from file, assume that all assigned adjacencies at this node have a weight in this file
    file = open(weighting,"r")

    threshold = float(threshold)
    for line in file:
        if line.startswith("#"):
            n = line.rstrip("\n")[1:]
            adjacenciesPerNode2[n] = []
        else:
            fields = line.split("\t")
            left = fields[0][1:]
            right = fields[1]
            weight = float(fields[3])
            adjacency = (left,right)
            rev = (right,left)
            nodes = []
            if adjacency in nodesPerAdjacency2:
                nodes = nodesPerAdjacency2[adjacency]
            elif rev in nodesPerAdjacency2:
                nodes = nodesPerAdjacency2[rev]

            if weight > threshold:
                adjacencyProbs2[n][(left,right)] = weight
                adjacenciesPerNode2[n].append((adjacency,weight))
                if nodes:
                    if not n in nodes:
                        nodesPerAdjacency2[adjacency].add(n)
                else:
                    nodesPerAdjacency2[adjacency] = set(n)
                    #nodes = nodesPerAdjacency2[adjacency]
            else:
                if nodes:
                    if n in nodes:
                        nodesPerAdjacency2[adjacency].remove(n)



    file.close()

    return nodesPerAdjacency2, adjacencyProbs2, adjacenciesPerNode2



###maybe useful one day, not used directly
def computeAdjacencyGraph(chromos):
    #keys: nodes, value: dict (keys: neighbors, value:additional infos if needed)
    adjacencygraph = {}
    for chro in chromos:
        orderList = chromos[chro]
        for i in range(0,len(orderList)-2):
            first = orderList[i]
            second = orderList[i+1]
            doubleFirst = doubleMarker(first)
            doubleSecond = doubleMarker(second)
            #add adjacency to adjacencygraph
            assert isinstance(chro, basestring)
            adjacencygraph[doubleFirst[1]] = {doubleSecond[0]: chro}
    return adjacencygraph
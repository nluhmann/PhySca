__author__ = 'Nina'
from ete2 import Tree
import sys
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse

# steps of action
# 1) reading tree - check
# 2) read order of unique universal marker and compute adjacencies - check
# 3) assign ancestral adjacencies (different possibilities!) - Dollo check
# 4) build global adjacency graph - check
# 5) identify CCs in global adjacency graph - check
# 6) run Sankoff for each CC - check
# 7) sampling


# parser = argparse.ArgumentParser(description="PhySca")
# parser.add_argument("marker", type=str, help="marker order of extant genomes")
# parser.add_argument("tree", type=str, help="tree file in newick format")
# parser.add_argument("alpha", type=float, help="alpha parameter in objective function, [0,1]")
# group = parser.add_mutually_exclusive_group(required=True)
# group.add_argument("-d", "--dollo", action="store_true", help="Assign potential adjacencies by dollo principle")
# group.add_argument("-x", type=float, help="Assign potential adjacencies by weight threshold, [0,1]")
# parser.add_argument("-kT", type=float, help="deClone constant", default=0.1)
# parser.add_argument("-s", "--sampling", type=int, help="sample X solutions for given set of parameters")
# args = parser.parse_args()


parser = argparse.ArgumentParser(description="PhySca")
parser.add_argument("-tree", type=str, help="tree file in newick format")
parser.add_argument("-alpha", type=float, help="alpha parameter in objective function, [0,1]")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-m", "--marker", type=str, help="marker order of extant genomes")
group.add_argument("-a", "--adjacencies", type=str, help="adjacencies in extant genomes")
group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument("-d", "--dollo", action="store_true", help="Assign potential adjacencies by dollo principle")
group2.add_argument("-x", type=float, help="Assign potential adjacencies by weight threshold, [0,1]")
group2.add_argument("-w", "--weight", type=str, help="weights for adjacencies at one internal nodes, adjacencies at other nodes get weight=0.5")
parser.add_argument("-gx", type=str, help="weights for adjacencies at one internal nodes, adjacencies at other nodes get weights computed by Declone, "
                                          "parameter -x must be given!")
parser.add_argument("-kT", type=float, help="deClone constant", default=0.1)
parser.add_argument("-s", "--sampling", type=int, help="sample X solutions for given set of parameters")
args = parser.parse_args()

#some checks for argument combinations:
# if args.marker:
#     if not args.dollo or not args.x:
#         parser.error('error')




# read tree file
file = open(args.tree, "r")
newick = file.read()
file.close()
tree = Tree(newick, format=1)
print tree.get_ascii()
print " "




if args.marker:
    #compute adjacencies, compute weights with DeClone

    print "Number of undoubled marker: "
    # read undoubled marker, for each species
    species_marker_order = {}
    file = open(args.marker, "r")
    species = ""
    for line in file:
        # new species
        if line.startswith(">"):
            if not species == "":
                species_marker_order[species] = chromosomes
                print species
                print markerCounter
            species = line.split("\t")[0][1:]
            chromosomes = {}
            markerCounter = 0
            # new chromosome
        elif line.startswith("#"):
            chrom = line.rstrip("\n")[2:]
        elif not line == "\n":
            order = line.rstrip("\n")[:-2].split(" ")
            markerCounter = markerCounter + len(order)
            chromosomes[chrom] = order
    species_marker_order[species] = chromosomes
    print species
    print markerCounter
    file.close()

    #compute adjacencies
    extantAdjacencies = getAdjacencies.findAdjacencies(species_marker_order)
    getAdjacencies.outputAdjacencies(extantAdjacencies)


    if args.dollo:
        #assign by dollo principle

        #compute potential ancestral adjacencies
        print "Assign potential adjacencies by Dollo principle..."
        pairPaths = getAdjacencies.findTreePaths(tree)
        adjacenciesAtNodes, nodesPerAdjacency = getAdjacencies.assignAncestralAdjacencies(pairPaths, extantAdjacencies, tree)
        getAdjacencies.outputAncestralAdjacencies(adjacenciesAtNodes,"ancestral_assigned_adjacencies")
        nodesPerAdjacencyTwo, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,0,tree,args.alpha,args.kT, args.tree)


    elif args.x != None:
        #assign by weight threshold
        print "Assign potential adjacencies that have a weight > "+str(args.x)
        nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,args.x,tree, args.alpha, args.kT, args.tree)
        getAdjacencies.outputAncestralAdjacenciesWithProbabilities(adjacenciesPerNode,"ancestral_assigned_adjacencies_with_weight")

    else:
        parser.error('error')



elif args.adjacencies:
    #adjacencies given,
    extantAdjacencies = getAdjacencies.readAdjacencyFile(args.adjacencies)

    if args.dollo:
        # compute weights with DeClone but assign by dollo principle
        #compute potential ancestral adjacencies
        print "Assign potential adjacencies by Dollo principle..."
        pairPaths = getAdjacencies.findTreePaths(tree)
        adjacenciesAtNodes, nodesPerAdjacency = getAdjacencies.assignAncestralAdjacencies(pairPaths, extantAdjacencies, tree)
        for node in adjacenciesAtNodes:
            print "Number of assigned adjacencies at "+node.name+": "+str(len(adjacenciesAtNodes[node]))
        getAdjacencies.outputAncestralAdjacencies(adjacenciesAtNodes,"ancestral_assigned_adjacencies")
        nodesPerAdjacencyTwo, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,0,tree,args.alpha,args.kT, args.tree)

    elif args.x != None:
        #compute weights with DeClone but assign weight threshold
        print "Assign potential adjacencies that have a weight > "+str(args.x)
        nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,args.x,tree, args.alpha, args.kT, args.tree)
        getAdjacencies.outputAncestralAdjacenciesWithProbabilities(adjacenciesPerNode,"ancestral_assigned_adjacencies_with_weight")
        for node in adjacenciesPerNode:
            print "Number of assigned adjacencies at "+node+": "+str(len(adjacenciesPerNode[node]))

    elif args.weight:
        #weights for one internal node given, assign 0.5 for all other adjacencies, Dollo assignment
        print "Assign other potential adjacencies by Dollo principle..."
        pairPaths = getAdjacencies.findTreePaths(tree)
        adjacenciesAtNodes, nodesPerAdjacency = getAdjacencies.assignAncestralAdjacencies(pairPaths, extantAdjacencies, tree)
        for node in adjacenciesAtNodes:
            print "number of assigned adjacencies at "+node.name+": "+str(len(adjacenciesAtNodes[node]))
        getAdjacencies.outputAncestralAdjacencies(adjacenciesAtNodes,"ancestral_assigned_adjacencies")
        #assign given weights for one node, 0.5 for all other nodes
        adjacencyProbs = getAdjacencies.computeWeightProbs(args.weight, extantAdjacencies, tree)


    elif args.gx:
        #weights for one internal node given, DeClone weights for all other adjacencies and assignment by threshold
        nodesPerAdjacency2, adjacencyProbs2, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,args.x,tree, args.alpha, args.kT, args.tree)
        nodesPerAdjacency, adjacencyProbs = getAdjacencies.assignGAMLweights(args.weight, nodesPerAdjacency2, adjacencyProbs2)


    else:
        parser.error('error, no mode given')




#nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode







#compute CCs in global adjacency graph
ccs = globalAdjacencyGraph.createGraph(extantAdjacencies,nodesPerAdjacency)
conflicts = globalAdjacencyGraph.analyseConnectedComponents(ccs)
globalAdjacencyGraph.outputConflicts(conflicts,"conflicts")

jointLabels, first = SR.enumJointLabelings(ccs)
validLabels, validAtNode = SR.validLabels(jointLabels,first)
#validAtNode = SR.useAllLabels(jointLabels,first,tree)



topDown = SR.computeLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, args.alpha)
reconstructedAdj = SR.reconstructedAdjacencies(topDown)
SR.outputReconstructedAdjacencies(reconstructedAdj,"reconstructed_adjacencies")
for node in reconstructedAdj:
    print node
    print "Number of reconstructed adjacencies: "+str(len(reconstructedAdj[node]))



scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
undoubled = scaffolding.undoubleScaffolds(scaffolds)
scaffolding.outputUndoubledScaffolds(undoubled,"undoubled_scaffolds")
scaffolding.outputScaffolds(scaffolds,"doubled_scaffolds")
scaffolding.sanityCheckScaffolding(undoubled)

for node in undoubled:
    print node
    markerCounter = 0
    print "number of scaffolds: "+str(len(undoubled[node]))
    for scaffold in undoubled[node]:
        markerCounter = markerCounter + len(scaffold)
    print "number of reconstructed undoubled marker in scaffolds: "+str(markerCounter)

if args.sampling:
    print "SAMPLING"
    for i in range(0,args.sampling):
        validLabels, validAtNode = SR.validLabels(jointLabels,first)
        print "###################### "+str(i)+" ######################"
        topDown = SR.sampleLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, args.alpha)
        reconstructedAdj = SR.reconstructedAdjacencies(topDown)
        SR.outputReconstructedAdjacencies(reconstructedAdj,"reconstructed_adjacencies_"+str(i))
        for node in reconstructedAdj:
            print node
            print "Number of reconstructed adjacencies: "+str(len(reconstructedAdj[node]))



        scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
        undoubled = scaffolding.undoubleScaffolds(scaffolds)
        scaffolding.outputUndoubledScaffolds(undoubled,"undoubled_scaffolds_"+str(i))
        scaffolding.outputScaffolds(scaffolds,"doubled_scaffolds_"+str(i))
        scaffolding.sanityCheckScaffolding(undoubled)

        for node in undoubled:
            print node
            markerCounter = 0
            print "number of scaffolds: "+str(len(undoubled[node]))
            for scaffold in undoubled[node]:
                markerCounter = markerCounter + len(scaffold)
            print "number of reconstructed undoubled marker in scaffolds: "+str(markerCounter)

__author__ = 'nluhmann'
from ete2 import Tree
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse
import time

from calculate_SCJ import calculate_SCJ


#TODO: Redesign the structure of the main with the new prprocessing script

parser = argparse.ArgumentParser(description="PhySca")
parser.add_argument("-tree", type=str, help="tree file in newick or nhx format")
parser.add_argument("-alpha", type=float, help="alpha parameter in objective function, [0,1]")
parser.add_argument("-extant",type=str,help="file with precomputed weighted adjacencies for external nodes")
parser.add_argument("-internal",type=str,help="file with precomputed weighted adjacencies for internal nodes")
#group = parser.add_mutually_exclusive_group(required=True)
#group.add_argument("-m", "--marker", type=str, help="marker order of extant genomes")
#group.add_argument("-a", "--adjacencies", type=str, help="adjacencies in extant genomes")
#group2 = parser.add_mutually_exclusive_group(required=True)
#group2.add_argument("-d", "--dollo", action="store_true", help="Assign potential adjacencies by dollo principle")
parser.add_argument("-x", type=float, help="Assign potential adjacencies by weight threshold, [0,1]",default=0.0)
#group2.add_argument("-x", type=float, help="Assign potential adjacencies by weight threshold, [0,1]")
#group2.add_argument("-w", "--weight", type=str, help="weights for adjacencies at specific internal nodes, adjacencies at other nodes get weight=0")
#parser.add_argument("-gx", type=str, help="weights for adjacencies at specific internal nodes, adjacencies at other nodes get weights computed by Declone, "
#                                          "parameter -x must be given!")
#parser.add_argument("-kT", type=float, help="deClone constant", default=0.1)
parser.add_argument("-s", "--sampling", type=int, help="sample X solutions for given set of parameters")
args = parser.parse_args()




t0 = time.time()


# read tree file
file = open(args.tree, "r")
newick = file.read()
file.close()
tree = Tree(newick, format=1)
print tree.get_ascii()
print " "




#if args.marker:
#    #compute adjacencies, compute weights with DeClone
#
#    print "Number of undoubled marker: "
#    # read undoubled marker, for each species
#    species_marker_order = {}
#    file = open(args.marker, "r")
#    species = ""
#    for line in file:
#        # new species
#        if line.startswith(">"):
#            if not species == "":
#                species_marker_order[species] = chromosomes
#                print species
#                print markerCount
#            species = line.split("\t")[0][1:]
#            chromosomes = {}
#            markerCount = 0
#            # new chromosome
#        elif line.startswith("#"):
#            chrom = line.rstrip("\n")[2:]
#        elif not line == "\n":
#            order = line.rstrip("\n")[:-2].split(" ")
#            markerCount = markerCount + len(order)
#            chromosomes[chrom] = order
#    species_marker_order[species] = chromosomes
#    print species
#    print markerCount
#    file.close()

    #compute adjacencies
#    extantAdjacencies = getAdjacencies.findAdjacencies(species_marker_order)
#    getAdjacencies.outputAdjacencies(extantAdjacencies)


    # if args.dollo:
    #     #assign by dollo principle
    #
    #     #compute potential ancestral adjacencies
    #     print "Assign potential adjacencies by Dollo principle..."
    #     pairPaths = getAdjacencies.findTreePaths(tree)
    #     adjacenciesAtNodes, nodesPerAdjacency = getAdjacencies.assignAncestralAdjacencies(pairPaths, extantAdjacencies, tree)
    #     getAdjacencies.outputAncestralAdjacencies(adjacenciesAtNodes,"ancestral_assigned_adjacencies")
    #     nodesPerAdjacencyTwo, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,0,tree,args.alpha,args.kT, args.tree)


#    if args.x != None:
#        #assign by weight threshold
#        print "Assign potential adjacencies that have a weight > "+str(args.x)
#        nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,args.x,tree, args.alpha, args.kT, args.tree)
#        getAdjacencies.outputAncestralAdjacenciesWithProbabilities(adjacenciesPerNode,"ancestral_assigned_adjacencies_with_weight")

#    else:
#        parser.error('error')



#elif args.adjacencies:
#    #adjacencies given,
#    extantAdjacencies = getAdjacencies.readAdjacencyFile(args.adjacencies)
#    #extantAdjacencies_species_adj=getAdjacencies.readAdjacencyFile_output_spec_adj(args.adjacencies)

#    # if args.dollo:
#    #     # compute weights with DeClone but assign by dollo principle
#    #     #compute potential ancestral adjacencies
#    #     print "Assign potential adjacencies by Dollo principle..."
#    #     pairPaths = getAdjacencies.findTreePaths(tree)
#    #     adjacenciesAtNodes, nodesPerAdjacency = getAdjacencies.assignAncestralAdjacencies(pairPaths, extantAdjacencies, tree)
#    #     for node in adjacenciesAtNodes:
#    #         print "Number of assigned adjacencies at "+node.name+": "+str(len(adjacenciesAtNodes[node]))
#    #     getAdjacencies.outputAncestralAdjacencies(adjacenciesAtNodes,"ancestral_assigned_adjacencies")
#    #     nodesPerAdjacencyTwo, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,0,tree,args.alpha,args.kT, args.tree)

#    if args.x != None and args.gx == None:
#        #compute weights with DeClone but assign weight threshold
#        print "Assign potential adjacencies that have a weight > "+str(args.x)
#        nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,args.x,tree, args.alpha, args.kT, args.tree)
#        getAdjacencies.outputAncestralAdjacenciesWithProbabilities(adjacenciesPerNode,"ancestral_assigned_adjacencies_with_weight")
#        for node in adjacenciesPerNode:
#            print "Number of assigned adjacencies at "+node+": "+str(len(adjacenciesPerNode[node]))

#    elif args.weight:
#        #weights for one internal node given, assign 0 for all other adjacencies

#        print "Assign potential adjacencies that have a weight > "+str(args.x)
#        x = 0
#        nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode = getAdjacencies.deCloneProbabilities(extantAdjacencies,x,tree, args.alpha, args.kT, args.tree)
#        #getAdjacencies.outputAncestralAdjacenciesWithProbabilities(adjacenciesPerNode,"ancestral_assigned_adjacencies_with_weight")
#        for node in adjacenciesPerNode:
#            print "Number of assigned adjacencies at "+node+": "+str(len(adjacenciesPerNode[node]))
#        adjacencyProbs = getAdjacencies.computeWeightProbs(args.weight, extantAdjacencies, tree)


#    elif args.gx:
#        #weights for one internal node given, DeClone weights for all other adjacencies and assignment by threshold
#        nodesPerAdjacency2, adjacencyProbs2, adjacenciesPerNode2 = getAdjacencies.deCloneProbabilities(extantAdjacencies,args.x,tree, args.alpha, args.kT, args.tree)
#        nodesPerAdjacency, adjacencyProbs, adjacenciesPerNode = getAdjacencies.assignGAMLweights(args.gx, nodesPerAdjacency2, adjacencyProbs2, adjacenciesPerNode2, args.x)
#        getAdjacencies.outputAncestralAdjacenciesWithProbabilities(adjacenciesPerNode,"ancestral_assigned_adjacencies_with_weight")

#    else:
#        parser.error('error, no mode given')

#TODO: create hash extantAdjacenies from input file weighted_extant_adjacencies

extantAdjacencies={}
f=open(args.extant,'r')
line=f.readline()
while line:
    #>species    (AdjL,AdjR)
    spec_adj_weight=line.split('\t')
    species=spec_adj_weight[0][1:]
    adj_L_R=spec_adj_weight[1].split(',')
    adj_L=adj_L_R[0][1:].replace("'","")
    adj_R=adj_L_R[1][:-1].replace("'","")
    adj=(int(adj_L),int(adj_R))
    weight=float(spec_adj_weight[2])
    if adj in extantAdjacencies:
        new_set=extantAdjacencies[adj]
        new_set.add(species)
        extantAdjacencies[adj]=new_set
    else:
        speciesSet=set()
        speciesSet.add(species)
        extantAdjacencies.update({adj:speciesSet})
    line=f.readline()
f.close()

#TODO: create hash nodesPerAdjacencies from input file weighted_internal_adjacencies
#TODO: create adjacencyProbs from input file weighted_internal_adjacencies

threshold=float(args.x)
adjacencyProbs={}
# structure Node:{(AdjL,AdjR):weight}
nodesPerAdjacency={}
#structure  (AdjL,AdjR):set(node)
f=open(args.internal,'r')
line=f.readline()
while line:
    #>species    (AdjL,AdjR)    weight
    spec_adj_weight=line.split('\t')
    species=str(spec_adj_weight[0][1:])
    adj_L_R=spec_adj_weight[1].split(',')
    adj_L=int(adj_L_R[0][1:])
    adj_R=int(adj_L_R[1][:-1])
    adj=(adj_L,adj_R)
    weight=float(spec_adj_weight[2])

    if weight > threshold:#filtering all adjacencies with a precomputed weight smaller than the threshold
        #filling nodesPerAdjaceny
        if adj in nodesPerAdjacency:
            new_set_A=nodesPerAdjacency[adj]
            new_set_A.add(species)
            nodesPerAdjacency[adj]=new_set_A
        else:
            spec=set()
            spec.add(species)
            nodesPerAdjacency.update({adj:spec})
        #filling adjacencyProbs with internal nodes
        if species in adjacencyProbs:
            adjacencyProbs[species][adj] = weight
        else:
            adjacencyProbs[species]={adj:weight}
    else:
        print(str(adj)+': '+str(weight)+' ist geringer als '+str(args.x))
    line=f.readline()
f.close()

#f=open("./mine_adjacencyProbs_pestis",'w')
#for species in adjacencyProbs:
#    for adj in adjacencyProbs[species]:
#        weight=adjacencyProbs[species][adj]
#        f.write(str(species)+'\t'+str(adj)+'\t'+str(weight)+'\n')
#f.close()


#compute CCs in global adjacency graph
ccs = globalAdjacencyGraph.createGraph(extantAdjacencies,nodesPerAdjacency)
conflicts = globalAdjacencyGraph.analyseConnectedComponents(ccs)
globalAdjacencyGraph.outputConflicts(conflicts,"conflicts")

jointLabels, first = SR.enumJointLabelings(ccs)
validLabels, validAtNode = SR.validLabels(jointLabels,first)



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
    for scaffold in undoubled[node]:
        first = scaffold[0]
        last = scaffold[-1]
        if not first == last:
            markerCounter = markerCounter + len(scaffold)
        else:
            markerCounter = markerCounter + len(scaffold)-1
    print node+" number of reconstructed undoubled marker in scaffolds: "+str(markerCounter)
    #if args.marker:
    #    notRec = markerCount - markerCounter
    #else:
    #    notRec = 2207 - markerCounter
    notRec = 2207 - markerCounter
    print node+" number of singleton scaffolds (not reconstructed marker): "+str(notRec)
    print node+" number of scaffolds: "+str(len(undoubled[node])+notRec)



print time.time() - t0, "seconds process time"




if args.sampling:
    print "SAMPLING"

    for i in range(0,args.sampling):
        t1 = time.time()
        jointLabels, first = SR.enumJointLabelings(ccs)
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
            for scaffold in undoubled[node]:
                first = scaffold[0]
                last = scaffold[-1]
                if not first == last:
                    markerCounter = markerCounter + len(scaffold)
                else:
                    markerCounter = markerCounter + len(scaffold)-1
            print node+" number of reconstructed undoubled marker in scaffolds: "+str(markerCounter)
            #if args.marker:
            #    notRec = markerCount - markerCounter
            #else:
            #    notRec = 2207 - markerCounter
            notRec = 2207 - markerCounter
            print node+" number of singleton scaffolds (not reconstructed marker): "+str(notRec)
            print node+" number of scaffolds: "+str(len(undoubled[node])+notRec)
        print time.time() - t1, "seconds process time"
        
        #calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
        
        




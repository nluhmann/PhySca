__author__ = 'nluhmann'
from ete2 import Tree
import globalAdjacencyGraph
import SR2
import scaffolding
import argparse
import time
import multiprocessing
import sys
import os
from calculate_SCJ import calculate_SCJ
import runSample
import collections

##########################
### PARSING PARAMETERS ###
##########################

# global variables
scj_path = './SCJ_distances'
statistic_path = './statistic_allSampled_ReconstructedAdjacencies'

parser = argparse.ArgumentParser(description="PhySca")
parser.add_argument("-tree", type=str, help="tree file in newick or nhx format")
parser.add_argument("-alpha", type=float, help="alpha parameter in objective function, [0,1] (default: 0)", default=0.0)
parser.add_argument("-extant", type=str, help="file with precomputed weighted adjacencies for external nodes")
parser.add_argument("-pot_extant", type=str, help="file with potential non-conflicting adjacencies for external nodes")
parser.add_argument("-internal", type=str, help="file with precomputed weighted adjacencies for internal nodes")
parser.add_argument("-x", type=float,
                    help="Filter potential adjacencies by weight threshold, [0,1] (default: no filtering)", default=-1)
parser.add_argument("-s", "--sampling", type=int, help="sample [INT] co-optimal solutions for given set of parameters")
parser.add_argument("-pN", "--processNumber", type=int,
                    help="number of processes used for parallel sampling. Max: [number of cpu] (default: 1)", default=1)
parser.add_argument("-out", "--output", type=str, help="specify output directory, current directory as default",
                    default=".")
parser.add_argument("-skip", "--skip", type=int,
                    help="skip connected components that consist of more edges than [INT] (default: no skipping)",
                    default=-1)
parser.add_argument("-sk", "--skip_first", action='store_const', const=False,
                    help="boolean, for skipping the 0th sampling")
parser.add_argument("-sc", "--start_counting", type=int,
                    help="specifies the starting number for enumerating the samples", default=0)
args = parser.parse_args()

t0 = time.time()

# create output directory if not existing
if not os.path.isdir(args.output):
    os.makedirs(args.output)

if not args.internal or not args.extant or not args.tree:
    # parser.error('Error: wrong parameter number or usage.')
    parser.print_help()
    sys.exit(1)

############################################
### GLOBAL PARAMETERS AND DATASTRUCTURES ###
############################################

threshold = float(args.x)

adjacencyProbs = {}
# structure {Node:{(AdjL,AdjR):weight,..},..}

extantAdjacencies = collections.defaultdict(set)
# structure {(AdjL,AdjR):set([node,node,..]),..}

potentialExtantAdjacencies = collections.defaultdict(set)
# structure {(AdjL,AdjR):set([node,node,..]),..}

extantAdjacencies_species_adj = collections.defaultdict(set)
# structure {node:set([(AdjL,AdjR),..]),..}

nodesPerAdjacency = collections.defaultdict(set)
# structure  (AdjL,AdjR):set(node)

# dictionary for all scj distances
dict_SCJ = {}

# dictionary for statistics of reconstructed Adjacencies
# structure: >internal node  adjacency   number of how often this adj was reconstructed at this node among all samples
allSampleReconstructionStatistic = {}

###########################
### PARSING INPUT FILES ###
###########################

# read tree file
file = open(args.tree, "r")
newick = file.read()
file.close()
tree = Tree(newick, format=1)
print tree.get_ascii()
print " "

# fill dictionaries with input from weighted_extant_adjacencies
# SE: so far, here we only store all observed extant adjacencies
f = open(args.extant, 'r')
line = f.readline()
while line:
    spec_adj = line.split('\t')
    species = spec_adj[0][1:]
    adj_L_R = spec_adj[1].strip().split(',')
    adj_L = int(adj_L_R[0][1:].replace("'", ""))
    adj_R = int(adj_L_R[1][:-1].replace("'", "").strip())
    adj = (adj_L, adj_R)

    #if weight > threshold: #filtering all adjacencies with a weight smaller than given threshold x
    extantAdjacencies[adj].add(species)

    extantAdjacencies_species_adj[species].add(adj)

    if species in adjacencyProbs:
        adjacencyProbs[species][adj] = 1.0  # set weight of external adjacencies in adjacencyProbs to 1
    else:
        adjacencyProbs[species] = {adj: 1.0}

    line = f.readline()
f.close()

if args.pot_extant:
    # SE: read all potential adjacencies into a separate hash
    f = open(args.pot_extant, 'r')
    for line in f:
        spec_adj = line.split('\t')
        species = spec_adj[0][1:]
        adj_L_R = spec_adj[1].strip().split(',')
        adj_L = adj_L_R[0][1:].replace("'", "")
        adj_R = adj_L_R[1][:-1].replace("'", "").strip()
        adj = (int(adj_L), int(adj_R))

        potentialExtantAdjacencies[adj].add(species)

        # SE: question: what weight should these adjacencies be assigned?
        # SE: CHECKPOINT
        if species in adjacencyProbs:
            adjacencyProbs[species][adj] = 1  # set weight of external adjacencies in adjacencyProbs to 1
        else:
            adjacencyProbs[species] = {adj: 1}
    f.close()


# TODO: compare how many conflicts adding potential adjacencies adds?
# TODO: two versions: treating adjacencies as missing or treating adjacencies as potential?
# TODO: check everything for input data on multiple chromosomes: if all of them are complete,
# TODO: will there be potential adjacencies merging two chromosomes falsely?


# create hash nodesPerAdjacencies from input file weighted_internal_adjacencies
# fill adjacencyProbs with input from file weighted_internal_adjacencies
# SE: weighted_internal_adjacencies will contain all adjacencies seen at any of the leaves, fragmented or not

# filteredAdjacencies={}
ancientLeaves = set()
f = open(args.internal, 'r')
line = f.readline()
while line:
    spec_adj_weight = line.split('\t')
    species = str(spec_adj_weight[0][1:])
    # check if we have an ancient leaf
    nodes = tree.search_nodes(name=species)
    #ToDo this gives an error if we read a newick tree with ancestor names!
    if nodes[0].is_leaf():
        ancientLeaves.add(species)

    adj_L_R = spec_adj_weight[1].split(',')
    adj_L = adj_L_R[0][1:]
    adj_R = adj_L_R[1][:-1]
    adj = (int(adj_L.strip()), int(adj_R.strip()))
    weight = float(spec_adj_weight[2])

    if weight > threshold:  # filtering all adjacencies with a precomputed weight smaller than the threshold
        # filling nodesPerAdjaceny
        nodesPerAdjacency[adj].add(species)

    # else:
    #     tup = (weight,adj)
    #     if species in filteredAdjacencies:
    #         filteredAdjacencies[species].append(tup)
    #     else:
    #         filteredAdjacencies[species]=[]
    #         filteredAdjacencies[species].append(tup)
    # filling adjacencyProbs with internal nodes
    if species in adjacencyProbs:
        adjacencyProbs[species][adj] = weight
    else:
        adjacencyProbs[species] = {adj: weight}

    line = f.readline()
f.close()

###############################
### COMPUTE RECONSTRUCTIONS ###
###############################

# compute CCs in global adjacency graph
# SE: this should not be influenced by ?-adjacencies, right?
ccs, removedAdjacencies, nodesPerAdjacency = globalAdjacencyGraph.createGraph(nodesPerAdjacency, args.skip)

if not args.skip_first:
    # identify conflicts in GAG for complexity analysis
    conflicts = globalAdjacencyGraph.analyseConnectedComponents(ccs)
    globalAdjacencyGraph.outputConflicts(conflicts, args.output + "/conflicts")

    # compute joint labels and check if they are valid (adjacencies conserved) at all internal nodes
    jointLabels, first = SR2.enumJointLabelings(ccs)
    #validLabels, validAtNode = SR2.validLabels(jointLabels, first)
    validAtNode = SR2.validLabels(jointLabels, first)

    ### comparing to previous extant genomes!
    #scaffolds_old = scaffolding.scaffoldAdjacencies(SR2.collectAdjacenciesPerNode(extantAdjacencies))
    #scaffolding.outputScaffolds(scaffolds_old, args.output + "/doubled_scaffolds_extant_old")

    # run Sankoff-Rousseau on all components
    reconstructedAdj,newExtant = SR2.computeLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, args.alpha,
                                            ancientLeaves, potentialExtantAdjacencies)

    SR2.outputReconstructedAdjacencies(reconstructedAdj, args.output + "/reconstructed_adjacencies")
    for node in reconstructedAdj:
        print node
        print "Number of reconstructed adjacencies: " + str(len(reconstructedAdj[node]))

    # scaffold adjacencies
    scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
    undoubled = scaffolding.undoubleScaffolds(scaffolds)
    scaffolding.outputUndoubledScaffolds(undoubled, args.output + "/undoubled_scaffolds")
    scaffolding.outputScaffolds(scaffolds, args.output + "/doubled_scaffolds")
    scaffolding.sanityCheckScaffolding(undoubled)



    # scaffold extant genomes
    scaffolds_extant = scaffolding.scaffoldAdjacencies(newExtant)
    undoubled_extant = scaffolding.undoubleScaffolds(scaffolds_extant)
    scaffolding.outputUndoubledScaffolds(undoubled_extant, args.output + "/undoubled_scaffolds_extant")
    scaffolding.outputScaffolds(scaffolds_extant, args.output + "/doubled_scaffolds_extant")
    scaffolding.sanityCheckScaffolding(undoubled_extant)



    # compute some stats
    # reconstruct marker pairs out of extantAdjacencies
    reconstructedMarker = set()
    for adj in extantAdjacencies:
        # each adjacency equals to markerpairs
        adj_list = [adj[0], adj[1]]
        for adjpart in adj_list:
            if adjpart % 2 == 0:
                markerId = adjpart / 2
            else:
                markerId = (adjpart + 1) / 2
            reconstructedMarker.add(markerId)

    for node in undoubled:
        markerCounter = 0
        for scaffold in undoubled[node]:
            first = scaffold[0]
            last = scaffold[-1]
            if not first == last:
                markerCounter += len(scaffold)
            else:
                markerCounter = markerCounter + len(scaffold) - 1
        print node + " number of reconstructed undoubled marker in scaffolds: " + str(markerCounter)
        # number of reconstructed markerIds
        reconstructedMarkerCount = len(reconstructedMarker)
        # singleton scaffolds number / number of not reconstructed marker
        notReconstructedMarkerCount = reconstructedMarkerCount - markerCounter
        # number of all scaffolds
        allScaffoldCount = len(undoubled[node]) + notReconstructedMarkerCount
        print node + " number of singleton scaffolds (not reconstructed marker): " + str(notReconstructedMarkerCount)
        print node + " number of scaffolds: " + str(allScaffoldCount)

    # calculate SCJ-distance for unsampled solution
    scj_unsampled = calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
    dict_SCJ.update({'Unsampled': scj_unsampled})
    print "Single-Cut-or-Join-Distance: " + str(scj_unsampled) + '\n'

print time.time() - t0, "seconds process time"
t_sampling = time.time()

#########################################
### SAMPLE CO-OPTIMAL RECONSTRUCTIONS ###
#########################################

# only if a number of samples is given and if the sript is called as standalone (not imported) script
if args.sampling and __name__ == '__main__':
    print "SAMPLING"
    # reconstruct marker pairs out of extantAdjacencies
    reconstructedMarker = set()
    for adj in extantAdjacencies:
        # each adjacency equals to markerpairs
        adj_list = [adj[0], adj[1]]
        for adjpart in adj_list:
            if adjpart % 2 == 0:
                markerId = adjpart / 2
            else:
                markerId = (adjpart + 1) / 2
            reconstructedMarker.add(markerId)
    reconstructedMarkerCount = len(reconstructedMarker)
    samplesize = args.processNumber
    # limiting processNumber on number of cpus
    cpuCount = multiprocessing.cpu_count()
    if samplesize > cpuCount:
        samplesize = cpuCount
    print "Using " + str(samplesize) + " parallel processes for sampling."
    print "Making " + str(args.sampling) + " samples"
    # create a pool of workerprocesses
    pool = multiprocessing.Pool(processes=samplesize)
    # create args.sampling tasks
    tasks = ((ccs, tree, extantAdjacencies, adjacencyProbs, args.alpha, i,
              extantAdjacencies_species_adj, args.output, reconstructedMarkerCount, ancientLeaves,
              potentialExtantAdjacencies) for i in range(args.start_counting, args.sampling + args.start_counting))
    # execute the sampling tasks
    results = pool.map_async(runSample.runSample, tasks)
    # close pool so no more tasks can be handed over to the workers
    pool.close()
    # let main program wait for every workerprocess
    pool.join()
    # receive the results
    output = results.get()

    # parse output and update allSampleReconstrutionStatistic and dit_SCJ with the results
    for tuple in output:
        tempRS = tuple[0]
        tempSCJ = tuple[1]
        tempOutLog = tuple[2]
        for key in tempRS:
            if key in allSampleReconstructionStatistic:
                allSampleReconstructionStatistic[key] += tempRS[key]
            else:
                allSampleReconstructionStatistic.update({key: tempRS[key]})
        for key in tempSCJ:
            dict_SCJ[key] = tempSCJ[key]
        print tempOutLog

###################################
### OUTPUT RECONSTRUCTION STATS ###
###################################

# write all SCJ distances to output file
f = open(args.output + "/" + scj_path, 'w')
for sample in dict_SCJ:
    f.write('>' + str(sample) + '\t' + str(dict_SCJ[sample]) + '\n')
f.close()

# write statistic about reconstructed adjacencies to output file
f = open(args.output + "/" + statistic_path, 'w')
f.write('>Header:' + '\t' + str(args.sampling) + '\n')

for node_adj in sorted(allSampleReconstructionStatistic.keys()):
    node = node_adj[0]
    adj = node_adj[1]
    number = allSampleReconstructionStatistic[node_adj]
    f.write('>' + str(node) + '\t' + str(adj) + '\t' + str(number) + '\n')

f.close()
print time.time() - t_sampling, "seconds process time"
print time.time() - t0, "seconds process time complete"

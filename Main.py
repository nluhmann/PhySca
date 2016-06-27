__author__ = 'nluhmann'
from ete2 import Tree
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse
import time
import multiprocessing

from calculate_SCJ import calculate_SCJ
import runSample

#global variables
scj_path='./SCJ_distances'
statistic_path='./statistic_allSampled_ReconstructedAdjacencies'

parser = argparse.ArgumentParser(description="PhySca")
parser.add_argument("-tree", type=str, help="tree file in newick or nhx format")
parser.add_argument("-alpha", type=float, help="alpha parameter in objective function, [0,1]",default=0.0)
parser.add_argument("-extant",type=str,help="file with precomputed weighted adjacencies for external nodes")
parser.add_argument("-internal",type=str,help="file with precomputed weighted adjacencies for internal nodes")
parser.add_argument("-x", type=float, help="Assign potential adjacencies by weight threshold, [0,1]",default=0.0)
parser.add_argument("-s", "--sampling", type=int, help="sample X solutions for given set of parameters")
parser.add_argument("-pN", "--processNumber", type=int, help="number of processes used for sampling. Max: [number of cpu]",default=1)
parser.add_argument("-out", "--output", type=str, help="specify output directory, current directory as default", default=".")
parser.add_argument("-sk", "--skip_first", action='store_const', const=True, help="boolean, for skipping the 0th sampling")
parser.add_argument("-sc", "--start_counting", type=int, help="specifies the starting number for enumerating the samples",default=0)
args = parser.parse_args()

t0 = time.time()

if not args.internal or not args.extant or not args.tree :
    parser.error('Error: wrong parameter number or usage.')

# read tree file
file = open(args.tree, "r")
newick = file.read()
file.close()
tree = Tree(newick, format=1)
print tree.get_ascii()
print " "
# from input file weighted_extant_adjacencies:
# create dict extantAdjacenies, containing external nodes
# create dict adjacencyProbs from input file
# create dict  extantAdjacencies_species_adj
threshold=float(args.x)

adjacencyProbs={}
# structure {Node:{(AdjL,AdjR):weight,..},..}

extantAdjacencies={}
#structure {(AdjL,AdjR):set([node,node,..]),..}

extantAdjacencies_species_adj={}
#structure {node:set([(AdjL,AdjR),..]),..}

#fill dictionaries with input from weighted_extant_adjacencies
f=open(args.extant,'r')
line=f.readline()
while line:
    #>species    (AdjL,AdjR)
    spec_adj=line.split('\t')
    species=spec_adj[0][1:]
    adj_L_R = spec_adj[1].strip().split(',')
    adj_L=adj_L_R[0][1:].replace("'","")
    adj_R=adj_L_R[1][:-1].replace("'","").strip()
    adj=(int(adj_L),int(adj_R))
    #if weight > threshold: #filtering all adjacencies with a weight smaller than given threshold x
    if adj in extantAdjacencies:
            new_set=extantAdjacencies[adj]
            new_set.add(species)
            extantAdjacencies[adj]=new_set
    else:
            speciesSet=set()
            speciesSet.add(species)
            extantAdjacencies.update({adj:speciesSet})
    if species in extantAdjacencies_species_adj:
            new_set=extantAdjacencies_species_adj[species]
            new_set.add(adj)
            extantAdjacencies_species_adj[species]=new_set
    else:
            adj_set=set()
            adj_set.add(adj)
            extantAdjacencies_species_adj[species]=adj_set
    if species in adjacencyProbs:
            adjacencyProbs[species][adj] = 1.0#set weight of external adjacencies in adjacencyProbs to 1
    else:
            adjacencyProbs[species]={adj:1.0}
    line=f.readline()
f.close()

# create hash nodesPerAdjacencies from input file weighted_internal_adjacencies
#fill adjacencyProbs with input from file weighted_internal_adjacencies

nodesPerAdjacency={}
#structure  (AdjL,AdjR):set(node)
f=open(args.internal,'r')
line=f.readline()
while line:
    #>species    (AdjL,AdjR)    weight
    spec_adj_weight=line.split('\t')
    species=str(spec_adj_weight[0][1:])
    adj_L_R=spec_adj_weight[1].split(',')
    adj_L=adj_L_R[0][1:]
    adj_R=adj_L_R[1][:-1]
    adj=(int(adj_L.strip()),int(adj_R.strip()))
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

    line=f.readline()
f.close()


#dictionary for all scj distances
dict_SCJ={}

#compute CCs in global adjacency graph
ccs = globalAdjacencyGraph.createGraph(extantAdjacencies,nodesPerAdjacency)
if (not args.skip_first):
    conflicts = globalAdjacencyGraph.analyseConnectedComponents(ccs)
    globalAdjacencyGraph.outputConflicts(conflicts,args.output+"/conflicts")

    jointLabels, first = SR.enumJointLabelings(ccs)
    validLabels, validAtNode = SR.validLabels(jointLabels,first)

    topDown = SR.computeLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, args.alpha)

    reconstructedAdj = SR.reconstructedAdjacencies(topDown)
    SR.outputReconstructedAdjacencies(reconstructedAdj,args.output+"/reconstructed_adjacencies")
    for node in reconstructedAdj:
        print node
        print "Number of reconstructed adjacencies: "+str(len(reconstructedAdj[node]))

    scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
    undoubled = scaffolding.undoubleScaffolds(scaffolds)
    scaffolding.outputUndoubledScaffolds(undoubled,args.output+"/undoubled_scaffolds")
    scaffolding.outputScaffolds(scaffolds,args.output+"/doubled_scaffolds")
    scaffolding.sanityCheckScaffolding(undoubled)

    # reconstruct marker pairs out of extantAdjacencies
    reconstructedMarker = set()
    for adj in extantAdjacencies:
        #each adjacency equals to markerpairs
        adj_list=[ adj[0], adj[1] ]
        for adjpart in adj_list:
            if (adjpart % 2 ==0):
                markerId=adjpart/2
            else:
                markerId=(adjpart+1)/2
            reconstructedMarker.add(markerId)


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
        # number of reconstructed markerIds
        reconstructedMarkerCount = len(reconstructedMarker)
        # singleton scaffolds number / number of not reconstructed marker
        notReconstructedMarkerCount = reconstructedMarkerCount - markerCounter
        # number of all scaffolds
        allScaffoldCount = len(undoubled[node]) + notReconstructedMarkerCount
        print node + " number of singleton scaffolds (not reconstructed marker): " + str(notReconstructedMarkerCount)
        print node + " number of scaffolds: " + str(allScaffoldCount)

        #calculate SCJ-distance for unsampled solution
        scj_unsampled=calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
        dict_SCJ.update({'Unsampled':scj_unsampled})

print time.time() - t0, "seconds process time"
t_sampling=time.time()

#dictionary for statistics of reconstructed Adjacencies
#structure: >internal node  adjacency   number of how often this adj was reconstructed at this node among all samples
allSampleReconstructionStatistic={}

#Sampling
#only if a number of samples is given and if the sript is called as standalone (not imported) script
if args.sampling and  __name__ == '__main__':
    print "SAMPLING"
    # reconstruct marker pairs out of extantAdjacencies
    reconstructedMarker = set()
    for adj in extantAdjacencies:
        # each adjacency equals to markerpairs
        adj_list = [adj[0], adj[1]]
        for adjpart in adj_list:
            if (adjpart % 2 == 0):
                markerId = adjpart / 2
            else:
                markerId = (adjpart + 1) / 2
            reconstructedMarker.add(markerId)
    reconstructedMarkerCount = len(reconstructedMarker)
    samplesize=args.processNumber
    #limiting processNumber on number of cpus
    cpuCount = multiprocessing.cpu_count()
    if samplesize > cpuCount:
        samplesize = cpuCount
    print "Using " + str(samplesize) + " parallel processes for sampling."
    print "Making " +str(args.sampling) + " samples"
    #create a pool of workerprocesses
    pool = multiprocessing.Pool(processes=samplesize)
    #create args.sampling tasks
    tasks = ((ccs, tree, extantAdjacencies, adjacencyProbs, args.alpha, i,
              extantAdjacencies_species_adj, args.output,reconstructedMarkerCount) for i in range(args.start_counting, args.sampling+args.start_counting))
    #execute the sampling tasks
    results = pool.map_async(runSample.runSample, tasks)
    #close pool so no more tasks can be handed over to the workers
    pool.close()
    #let main programm wait for every workerprocess
    pool.join()
    #receive the results
    output = results.get()

    #parse output and update allSampleReconstrutionStatistic and dit_SCJ with the results
    for tuple in output:
        tempRS = tuple[0]
        tempSCJ = tuple[1]
        tempOutLog=tuple[2]
        for key in tempRS:
            if key in allSampleReconstructionStatistic:
                allSampleReconstructionStatistic[key] += tempRS[key]
            else:
                allSampleReconstructionStatistic.update({key: tempRS[key]})
        for key in tempSCJ:
            dict_SCJ[key] = tempSCJ[key]
        print tempOutLog


#write all SCJ distances to output file
f=open(args.output+"/"+scj_path,'w')
for sample in dict_SCJ:
    f.write('>'+str(sample)+'\t'+str(dict_SCJ[sample])+'\n')
f.close()

#write statistic about reconstructed adjacencies to output file

f=open(args.output+"/"+statistic_path,'w')
f.write('>Header:'+'\t'+str(args.sampling)+'\n')

for node_adj in sorted(allSampleReconstructionStatistic.keys()):
    node=node_adj[0]
    adj=node_adj[1]
    number=allSampleReconstructionStatistic[node_adj]
    f.write('>'+str(node)+'\t'+str(adj)+'\t'+str(number)+'\n')

f.close()
print time.time() - t_sampling, "seconds process time"
print time.time() - t0, "seconds process time complete"


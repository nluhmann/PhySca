__author__ = 'nluhmann'
from ete2 import Tree
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse
import time

from calculate_SCJ import calculate_SCJ

parser = argparse.ArgumentParser(description="PhySca")
parser.add_argument("-tree", type=str, help="tree file in newick or nhx format")
parser.add_argument("-alpha", type=float, help="alpha parameter in objective function, [0,1]",default=0.0)
parser.add_argument("-extant",type=str,help="file with precomputed weighted adjacencies for external nodes")
parser.add_argument("-internal",type=str,help="file with precomputed weighted adjacencies for internal nodes")
parser.add_argument("-x", type=float, help="Assign potential adjacencies by weight threshold, [0,1]",default=0.0)
parser.add_argument("-s", "--sampling", type=int, help="sample X solutions for given set of parameters")
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

#structure
f=open(args.extant,'r')
line=f.readline()
while line:
    #>species    (AdjL,AdjR)
    #spec_adj_weight=line.split('\t')
    #species=spec_adj_weight[0][1:]
    #adj_L_R=spec_adj_weight[1].split(',')
    #adj_L=adj_L_R[0][1:].replace("'","").strip()
    #adj_R=adj_L_R[1][:-1].replace("'","").strip()
    #adj=(adj_L,adj_R)
    #weight=float(spec_adj_weight[2])
    spec_adj=line.split('\t')
    species=spec_adj[0][1:]
    adj_L_R = spec_adj[1].strip().split(',')
    adj_L=adj_L_R[0][1:].replace("'","")
    adj_R=adj_L_R[1][:-1].replace("'","").strip()
    adj=(adj_L,adj_R)
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
    adj=(adj_L.strip(),adj_R.strip())
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
    else:
        print(str(adj) + ': ' + str(weight) + ' ist geringer als ' + str(args.x))
    if species in adjacencyProbs:
            adjacencyProbs[species][adj] = weight
    else:
            adjacencyProbs[species]={adj:weight}

    line=f.readline()
f.close()




#f=open('single_leaf_adjacencies','r')
#line=f.readline()
#while line:
#    adj_spec_weight=line.split('\t')
#    adj=adj_spec_weight[0].split(',')
#    adj_L=adj[0][1:]
#    adj_R=adj[1][:-1]
#    adj=(adj_L,adj_R)
#    species=adj_spec_weight[1].split(',')[0][1:].replace("'","")
#    weight=float(adj_spec_weight[2])
#    if species in adjacencyProbs:
#        adjacencyProbs[species][adj] = weight
#    else:
#        adjacencyProbs[species] = {adj: weight}
#    line=f.readline()
#f.close()

f=open("./mine_adjacencyProbs_pestis",'w')
for species in adjacencyProbs:
    for adj in adjacencyProbs[species]:
        weight=adjacencyProbs[species][adj]
        f.write(str(species)+'\t'+str(adj)+'\t'+str(weight)+'\n')
f.close()

f=open("./mine_extantAdj_pestis",'w')
for adj in extantAdjacencies:
    for spec in extantAdjacencies[adj]:
        f.write(str(adj)+'\t'+str(spec)+'\n')
f.close()

f=open("./mine_nodesPerAdjacency_pestis",'w')
for adj in nodesPerAdjacency:
    for spec in nodesPerAdjacency[adj]:
        f.write(str(adj)+'\t'+str(spec)+'\n')
f.close()

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
        
        calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
        
        




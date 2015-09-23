__author__ = 'Nina'

import networkx as nx
import itertools
import math

# extant:
# keys: (left marker, right marker), value: [(species,chromosome),...]
# ancestral:
# key: adjacency, value: list of internal nodes
# graph: #keys: nodes, value: dict (keys: neighbors, value:additional infos if needed)
def createGraph(extant,ancestral):
    # test = nx.Graph()
    # for adja in extant:
    #     test.add_edge(adja[0],adja[1])
    # #connected components as subgraphs
    # testgraphs = [c for c in sorted(nx.connected_component_subgraphs(test), key=len, reverse=True)]
    # #print nx.all_pairs_shortest_path(graphs[0])
    # print "biggest CC, number of nodes extant: "+str(nx.number_of_nodes(testgraphs[0]))
    # print "biggest CC, number of edges extant: "+str(nx.number_of_edges(testgraphs[0]))

    print "Create global adjacency graph..."
    G=nx.Graph()

    for adj in ancestral:
        ancestralSpecies = ancestral[adj]
        G.add_edge(adj[0],adj[1], species=ancestralSpecies)
        #note to Nina: if edges gets more than one annotation in species, they are put in a set, everything okay!

    # counter = 0
    # for exAdj in extant:
    #     if G.has_edge(exAdj[0],exAdj[1]):
    #         speciesSet = G[exAdj[0]][exAdj[1]]["species"]
    #     else:
    #         G.add_edge(exAdj[0],exAdj[1])
    #         counter = counter +1
    #         speciesSet = set()
    #     specs = extant[exAdj]
    #     for elem in specs:
    #         speciesSet.add(elem[0])
    #     G[exAdj[0]][exAdj[1]]["species"] = speciesSet

    print " "
    print "Global adjacency graph:"
    print "Number of edges: "+str(G.number_of_edges())
    print "Number of nodes: "+str(G.number_of_nodes())
    print "Number of connected components: "+str(nx.number_connected_components(G))
    length = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]


    #connected components as subgraphs
    graphs = [c for c in sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)]

    #compute complexity for biggest CC in global adjacency graph:
    numberOfLabels = 1
    for node in nx.nodes(graphs[0]):
        neighbors = nx.neighbors(graphs[0],node)
        numberOfLabels = numberOfLabels * len(neighbors)

    complex = 7 * math.pow(numberOfLabels,2)
    print "Number of Labels: "+str(numberOfLabels)
    print "Complexity: "+str(complex)
    if len(graphs) > 10:
        for i in range(0,10):
            #print str(i)+" CC, number of nodes: "+str(nx.number_of_nodes(graphs[i]))
            print str(nx.number_of_edges(graphs[i]))
    return graphs



def outputConnectedComponents(graphs):
    x=1


def analyseConnectedComponents(graphs):
    conflicts = {}
    for cc in graphs:
        nodes = nx.nodes(cc)
        for node in nodes:
            if nx.degree(cc,node) > 1:
                edges = nx.edges(cc,node)
                for pair in itertools.combinations(edges,2):
                    leftSet = cc[pair[0][0]][pair[0][1]]["species"]
                    left = set()
                    for elem in leftSet:
                        left.add(elem)
                    rightSet = cc[pair[1][0]][pair[1][1]]["species"]
                    right = set()
                    for elem in rightSet:
                        right.add(elem)
                    #identify conflicts
                    if not left.isdisjoint(right):
                        conflict = left.intersection(right)
                        for elem in conflict:
                            if elem in conflicts:
                                conflicts[elem].append(pair)
                            else:
                                conflicts[elem] = [pair]
    #combineConflicts(conflicts)
    print " "
    for species in conflicts:
        print "Number of conflicts in "+species+" : "+str(len(conflicts[species]))
    return conflicts


def outputConflicts(conflicts,out):
    file = open(out,"w")
    for confl in conflicts:
        file.write(">"+confl+" "+str(len(conflicts[confl]))+"\n")
        for elem in conflicts[confl]:
            file.write(str(elem)+"\n")


    file.close()

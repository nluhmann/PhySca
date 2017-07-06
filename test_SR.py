from unittest import TestCase
import networkx as nx
import SR2


__author__ = 'Nina'


class TestSR(TestCase):
    graph = nx.Graph()
    graph.add_edge("a","c",species=set("1"))
    graph.add_edge("a","b",species=set("4"))
    graph.add_edge("a","d",species=set("2"))
    graph.add_edge("b","d",species=set("3"))
    graph.add_edge("b","e",species=set("5"))
    graph.add_edge("f","e",species=set("1"))


    # def test_enumJointLabelings(self):
    #     #print nx.maximal_matching(self.graph)
    #     #print nx.max_weight_matching(self.graph)
    #     joint,first = SR2.enumJointLabelings([self.graph])
    #     TestCase.assertEqual(self,len(joint[self.graph]),15)
    #     valid, validAtNode = SR2.validLabels(joint,first)
    #     TestCase.assertEqual(self,len(valid[self.graph]),8)


    # def test_cost(self):
    #     label1 = [('a', '-'), ('b', '-'), ('c', '-'), ('d', '-'), ('e', '-'), ('f', '-')]
    #     label2 = [('a', ('a', 'c')), ('b', '-'), ('c', ('a', 'c')), ('d', '-'), ('e', ('e', 'f')), ('f', ('e', 'f'))]
    #     label3 = [('a', ('a', 'c')), ('b', '-'), ('c', ('a', 'c')), ('d', '-'), ('e', '-'), ('f', '-')]
    #     label4 = [('a', ('a', 'd')), ('b', '-'), ('c', '-'), ('d', ('a', 'd')), ('e', '-'), ('f', '-')]
    #     cost1 = SR.cost(label1,label2,self.graph.edges())
    #     #TestCase.assertEqual(self,cost1,3)
    #     cost2 = SR.cost(label2,label3,self.graph.edges())
    #     #TestCase.assertEqual(self,cost2,3)
    #     cost3 = SR.cost(label3,label4,self.graph.edges())
    #     #TestCase.assertEqual(self,cost3,3.5)
    #     labelx = [('8', ('8','9')), ('9', ('8','9'))]
    #     labely = [('8', ('8','9')), ('9', ('8','9'))]
    #     cost4 = SR.cost(labelx,labely,[('8','9')])
    #     TestCase.assertEqual(self,cost4,0)

    # def test_computeLabeling(self):
    #     chrom = {}
    #     chrom["one"] = ["3","4","5","6","7","8"]
    #     chrom["two"] = ["1","2"]
    #     species = {}
    #     species["A"] = chrom
    #     chrom = {}
    #     chrom["one"] = ["-4","5","-3","6","7","8"]
    #     chrom["two"] = ["1","2"]
    #     species["B"] = chrom
    #     chrom = {}
    #     chrom["one"] = ["3","4","5","6","7","8"]
    #     chrom["two"] = ["1","2"]
    #     species["E"] = chrom
    #     chrom = {}
    #     chrom["one"] = ["4","5","-3","6","7","8"]
    #     chrom["two"] = ["1","2"]
    #     species["D"] = chrom
    #
    #     adj = getAdjacencies.findAdjacencies(species)
    #     tree = Tree("(A:1,(B:1,(E:1,D:1)Int_1:0.5[&&NHX:ancient=1])Int_2:0.5[&&NHX:ancient=0])Root;", format=1)
    #     paths = getAdjacencies.findTreePaths(tree)
    #     internal,adjacenciesAncestral = getAdjacencies.assignAncestralAdjacencies(paths,adj,tree)
    #     graphs = globalAdjacencyGraph.createGraph(adj,adjacenciesAncestral)
    #     jointLabels, first = SR.enumJointLabelings(graphs)
    #     validLabels, validAtNode = SR.validLabels(jointLabels,first)
    #
    #     resolvedCCs = SR.computeLabelings(tree, graphs, validAtNode, adj)
    #     reconstructedAdj = SR.reconstructedAdjacencies(resolvedCCs)
    #     print reconstructedAdj


    # def test_sampleLabelings(self):
    #     tree = Tree("(A:1,(B:1,(C:1,(E:1,D:1)Int_1:0.5[&&NHX:ancient=1])Int_2:0.5[&&NHX:ancient=0])Int_3:1)Root;", format=1)
    #     chrom = {}
    #     chrom["one"] = ["3","4"]
    #     species = {}
    #     species["C"] = chrom
    #     chrom = {}
    #     chrom["one"] = ["3","4"]
    #     species["D"] = chrom
    #     chrom = {}
    #     chrom["one"] = []
    #     species["E"] = chrom
    #     chrom = {}
    #     chrom["one"] = []
    #     species["A"] = chrom
    #     chrom = {}
    #     chrom["one"] = []
    #     species["B"] = chrom
    #
    #     adj = getAdjacencies.findAdjacencies(species)
    #     paths = getAdjacencies.findTreePaths(tree)
    #     internal,adjacenciesAncestral = getAdjacencies.assignAncestralAdjacencies(paths,adj,tree)
    #     graphs = globalAdjacencyGraph.createGraph(adj,adjacenciesAncestral)
    #     jointLabels, first = SR2.enumJointLabelings(graphs)
    #     probs={"Int_1":{(6, 7):0.1},"Int_2":{(6, 7):0.1},"Int_3":{(6, 7):0.1},"Root":{(6, 7):0.1}}
    #     for i in range(0,10):
    #         validLabels, validAtNode = SR2.validLabels(jointLabels,first)
    #         resolvedCCs = SR2.sampleLabelings(tree, graphs, validAtNode, adj,probs, alpha=0)
    #         reconstructedAdj = SR2.reconstructedAdjacencies(resolvedCCs)
    #         print reconstructedAdj
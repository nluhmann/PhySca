from unittest import TestCase
import scaffolding

__author__ = 'Nina'


class TestScaffolding(TestCase):
    def test_getOtherExtremity(self):
        other = scaffolding.getOtherExtremity("3")
        TestCase.assertEqual(self,other,4)
        other = scaffolding.getOtherExtremity(4)
        TestCase.assertEqual(self,other,3)


    def test_scaffoldAdjacencies(self):
        reconstructedAdj = {'Test':[(16,17),(22,24),(2,3),(28,29),(8,9),(11,13),(18,20),(19,21),(4,5),(26,27)]}
        scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
        res = {'Test':{(25, 30): [25, 26, 27, 28, 29, 30], (12, 14): [12, 11, 13, 14], (1, 6): [1, 2, 3, 4, 5, 6], (7, 10): [7, 8, 9, 10], (15, 23): [15, 16, 17, 18, 20, 19, 21, 22, 24, 23]}}
        TestCase.assertEqual(self,scaffolds,res)

    def test_mergeScaffolds(self):
        scaffold = {(1,8):[1,2,3,4,5,6,7,8],(9,7):[9,10,8,7],(14,10):[14,13,9,10],(11,13):[11,12,14,13],(2,16):[2,1,15,16]}
        merged = scaffolding.mergeScaffolds(scaffold)
        print merged
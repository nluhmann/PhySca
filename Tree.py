__author__ = 'Nina'
from ete2 import Tree




# You can also specify the newick format. For instance, for named internal nodes we will use format 1.



# for n in t.traverse(strategy="postorder"):
#     if hasattr(n,"ancient"):
#         print n.name
#         print n.ancient




def annotateleaves(t,labels):
    """

    :rtype : tree
    """
    for leaf in t.traverse():
        if leaf.is_leaf():
           leaf.add_feature("annotation",labels[leaf.name])
    return t




def annotationTest(t,key):
    for n in t.traverse():
        if hasattr(n,key):
            print n.name
            print getattr(n,key)
            print "---------"

def costMatrixToTree(t,space):
    matrix = {}
    for elem in space:
        matrix[elem] = {}
        for another in space:
            if not elem == another:
                matrix[elem][another] = 1
    for node in t.traverse():
        if not node.is_leaf():
            node.add_feature("matrix",matrix)

    return t

def sankoff_bottomup(t,space):
    """

    :rtype : tree
    """
    for node in t.traverse(strategy="postorder"):
        #get annotation for all children of node
        if not node.is_leaf():
            annoHash = {}
            for label in space:
                minTotal = 0
                allMin = float("inf")
                minLabel = ""
                for child in node.children:
                    minChild = float("inf")
                    anno = child.annotation
                    for annotation in anno:
                        #compute cost from label to annotation!
                        if not annotation == label:
                            value = anno[annotation] + node.matrix[label][annotation]
                        else:
                            value = anno[annotation]
                        if value < minChild:
                            minChild = value
                    minTotal = minTotal + minChild
                annoHash[label] = minTotal
                if minTotal < allMin:
                    minLabel = label
            node.add_feature("annotation",annoHash)
            node.add_feature("minimumLabel",minLabel)
    return t


def sankoff_topdown(t):
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
                minCost = float("inf")
                if anno == parentAssignment:
                    minCost = annotation[parentAssignment]
                else:
                    minCost = annotation[anno] + node.up.matrix[anno][parentAssignment]
                if minCost < minAllCost:
                    minAllCost = minCost
                    minAllLabel = anno
            node.add_feature("assignment",minAllLabel)
    return t



def test():
    tree = Tree("(A:1,(B:1,(E:1,D:1)Int_1:0.5[&&NHX:ancient=1])Int_2:0.5[&&NHX:ancient=0])Root;", format=1)
    print tree.get_ascii()

    # 1) annotate each leaf with its label, leafs do not need to be annotated with all possible labels, we just assume that
#    each label that is not present in an annotation has cost infinity

    leafLabels = {
            "A":{"A":0},
            "B":{"Y":0},
            "E":{"X":0},
            "D":{"Y":0}}

    #all labels possible
    labelSpace = ["X","Y","A"]

    assert isinstance(tree, Tree)
    assert isinstance(leafLabels,dict)
    annotatedTree = annotateleaves(tree,leafLabels)
    matrixTree = costMatrixToTree(annotatedTree,labelSpace)
    bottomTree = sankoff_bottomup(matrixTree,labelSpace)

    topTree = sankoff_topdown(bottomTree)

    annotationTest(topTree, "annotation")

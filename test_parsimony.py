# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2015-08-05 12:49
# modified : 2015-08-05 12:49
"""
<++>
"""

__author__="Johann-Mattis List"
__date__="2015-08-05"


from lingpyd import *
from newick import *
from parsimony import *

patterns = [['b'],['c'],['a'],['a'],['b']]
taxa = ['A','B','C','D','E']
tree = LingPyTree('(((A,B),(C,D)),E);')
tree2 = LingPyTree('((A,E),(B,C),D);')
transitions = [
         # a b c ø
        [0, 2, 4, 1], # a
        [2, 0, 6, 0], # b
        [4, 6, 0, 10], # c
        [1, 0, 10, 0], # ø
        ]
characters = ['a','b','c',"ø"]

W = sankoff_parsimony_up(patterns, taxa, tree, transitions, characters)
S = sankoff_parsimony_down(
        W, patterns, taxa, tree, transitions, characters)
    
a,b,c = sankoff_parsimony(patterns, taxa, tree, transitions, characters,
        pprint=True)


patterns = [['ABB'],['ABB'],['ABB'],['ABB'],['BBB','BCB','BCC','BBC','CCC','CCB','CBC']]
chars = []
for p in patterns:
    chars += p
chars = sorted(set(chars))
matrix = [[0 for i in chars] for j in chars]
taxa = ['Cone','Zhonggu','Thebo','Purik','Shigatse']

def fill_matrix(chars, matrix):
    for i,a in enumerate(chars):
        for j,b in enumerate(chars):
            if i < j:
                d = edit_dist(a,b)
                matrix[i][j] = d
                matrix[j][i] = d

fill_matrix(chars,matrix)

a,b = best_tree_brute_force(
        [patterns, patterns],
        taxa,
        [matrix, matrix],
        [chars, chars]
        )

print(Tree(a[0][0]).asciiArt())
print(b)
print(len(a))


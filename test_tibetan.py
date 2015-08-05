# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2015-08-05 17:37
# modified : 2015-08-05 17:37
"""
<++>
"""

__author__="Johann-Mattis List"
__date__="2015-08-05"

from lingpy import *
import re
import itertools

from newick import *
from parsimony import *

# prepare the data
data = csv2list('data-tibetan.tsv')


allchars = [''.join(x) for x in itertools.product('ABC','ABC','ABC')]

def fill_matrix(chars, matrix):
    for i,a in enumerate(chars):
        for j,b in enumerate(chars):
            if i < j:
                d = edit_dist(a,b)
                matrix[i][j] = d
                matrix[j][i] = d

header = data[0]
data = data[1:]

taxa = header[2:]
matrices = []
concepts = []

D = {}
C = {}

# now let's make the matrix
for line in data:
    
    tmp = line[2:]
    chars = []

    C[line[0]] = line[1]

    for char in tmp:

        if char != '-':
            out = []
            bracket = False
            for c in char:
                if c == '(':
                    bracket = True
                    out += ['']
                elif c == ')':
                    bracket = False

                if not bracket and c not in ")('":
                    out += [c]
                elif bracket and c not in '()':
                    out[-1] += c
                elif c == "'":
                    out[-1] += c

            # check for general bracket 
            if '(' in char:
                new_chars = [''.join(x) for x in 
                        itertools.product(
                            *out
                            )]
            else:
                new_chars = [''.join(out)]

            chars += [new_chars]
        else:
            chars += ['-']

    # correct for all combinations
    D[line[0]] = chars

# now iterate and create the matrix
matrices = []
patterns = []
chars = []
for k in sorted(D):
    uchars = []
    for pt in D[k]:
        for p in pt:
            uchars += [p]
    uchars = sorted(set(uchars + allchars))
    matrix = [[0 for i in uchars] for j in uchars]
    fill_matrix(uchars, matrix)

    matrices += [matrix]
    patterns += [D[k]]
    chars += [uchars]
    concepts += [k]

protos = [C[c] for c in concepts]

a,b = best_tree_brute_force(
        patterns,
        taxa,
        matrices,
        chars,
        verbose=True,
        proto_forms = protos
        )

print(Tree(a[0][0]).asciiArt())
print(b)
print(len(a))
        

for t,s in a:

    tree = LingPyTree(t)
    # get the "proto-form"
    match = 0
    smatch = 0
    mmatch = 0
    
    # check parsimony
    for p,m,c,k in zip(patterns, matrices, chars, concepts):

        W = sankoff_parsimony_up(p,taxa,tree, m, c)
        
        smin = min(W[tree.root].values())
        proto = [a for a,b in W[tree.root].items() if b == smin]

        rprot = C[k]

        if C[k] in proto:
            if len(proto) == 1:
                match += 1
                #print('[*] '+rprot)
            else:
                smatch += 1
                #print('[+] '+' '.join(proto)+' / '+rprot)
        else:
            mmatch += 1
            #print('[!] '+ ' '.join(proto) + ' / ' + rprot)
    
    print(Tree(t).asciiArt())
    print(t, match, smatch, mmatch)
    

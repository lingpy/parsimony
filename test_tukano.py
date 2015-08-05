# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2015-08-05 21:02
# modified : 2015-08-05 21:02
"""
<++>
"""

__author__="Johann-Mattis List"
__date__="2015-08-05"

import json

out = json.loads(open('tukano.simple.json').read())

from lingpy import *
from newick import *
from parsimony import *

tree = LingPyTree(open('tukano.tre').read().strip())

match = 0
smatch = 0
mmatch = 0
for p,m,c,pr in zip(out['patterns'], out['matrix'], out['chars'],
        out['protos']):

    w = sankoff_parsimony_up(
            p,
            tree.leaves,
            tree,
            m,
            c
            )

    # get best match
    minw = min(w[tree.root].values())
    minc = [a for a,b in w[tree.root].items() if b == minw]


    if [x for x in pr[0].split(' / ') if x in minc]:
        if len(minc) == 1:
            match += 1
            print('[*] '+pr[0], minw)
        else:
            smatch += 1
            print('[+] '+pr[0]+' / '+','.join(minc), minw)
    else:
        mmatch += 1
        print('[!] '+pr[0] + ' / '+','.join(minc), minw)

print(match, smatch, mmatch)



sankoff_parsimony(
        out['patterns'][12],
        tree['leaves'],
        tree,
        out['matrix'][12],
        out['chars'][12],
        pprint = True
        )

outs = best_tree_brute_force(
        out['patterns'],
        tree.leaves,
        out['matrix'],
        out['chars'],
        proto_forms = False,
        verbose=True
        )
print(outs[0])
print('')
print(outs[1])

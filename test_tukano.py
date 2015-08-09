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



#sankoff_parsimony(
#        out['patterns'][12],
#        out['taxa'],
#        tree,
#        out['matrix'][12],
#        out['chars'][12],
#        pprint = True
#        )

#outs = best_tree_brute_force(
#        out['patterns'],
#        tree.leaves,
#        out['matrix'],
#        out['chars'],
#        proto_forms = False,
#        verbose=True
#        )
#print(outs[0])
#print('')
#print(outs[1])



#trees = branch_and_bound(
trees = heuristic_parsimony(        
        out['taxa'],
        out['patterns'],
        out['matrix'],
        out['chars'],
        lower_bound = 200,
        #guide_tree = "((((((((((Bar,Tat),(Kub,Tan)),Kar),Pis),Tuy,Yur),(Pir,(Tuk,Wan))),((Des,Yup),Sir)),((Kor,(Sek,Sio)),Kue)),Mak),Bas);", #tree.newick,
        iterations = 50000,
        sample_steps = 100,
        )


with open('tukano.trees','w') as f:
    for t in trees[0]:
        f.write(t+'\n')
    f.write(str(trees[1]))
        


import networkx as nx
taxa = out['taxa']
matrices = out['matrix']
chars = out['chars']
patterns = out['patterns']

G = nx.Graph()
for i,tA in enumerate(taxa):
    for j,tB in enumerate(taxa):
        if i < j:

            # get score between two taxa
            all_scores = []
            for pattern,matrix,charset in zip(patterns, matrices, chars):

                patternA = pattern[i]
                patternB = pattern[j]

                # calculate minimal weight
                scores = []
                for pA in patternA:
                    pAidx = charset.index(pA)
                    for pB in patternB:
                        pBidx = charset.index(pB)

                        score = matrix[pAidx][pBidx]
                    scores += [score]
                all_scores += [min(scores)]

            G.add_edge(tA, tB, weight=sum(all_scores))

g = nx.minimum_spanning_tree(G)
print(sum([w[2]['weight'] for w in g.edges(data=True)]))

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
import networkx as nx


# create tukano data
sc = csv2list('tukano.network', strip_lines=False)

D = {}
for line in sc[1:]:
    no = line[0]
    p = line[1].strip()
    c = line[2].strip()
    s,t = line[-2].strip(), line[-1].strip()
    if not s or not t:
        pass
    else:
        if (p,c) not in D:
            D[p,c] = [(s,t)]
        else:
            D[p,c] += [(s,t)]

for a,b in out['protos']:
    if (a,b) not in D:
        print(a,b)

for idx,(a,b) in enumerate(out['protos']):

    if (a,b) in D:
        # make the graph
        g = nx.DiGraph()
        g.add_edges_from(D[a,b])
        
        nlen = len(g)
        matrix = [[0 for x in range(nlen)] for y in range(nlen)]
        chars = sorted(g.nodes())

        for i,nA in enumerate(chars):
            for j,nB in enumerate(chars):
                try:
                    d = nx.shortest_path_length(g,nA,nB)
                except nx.NetworkXNoPath:
                    d = 10 * len(chars)
                
                matrix[i][j] = d

                if d < 10:
                    print(nA, '->', nB, ':', d)
        out['matrix'][idx] = matrix
        out['chars'][idx] = chars
input('pause')
    
    

tree = LingPyTree(open('tukano.trees').read().split('\n')[3].strip())

match = 0
smatch = 0
mmatch = 0
for idx,(p,m,c,pr) in enumerate(zip(out['patterns'], out['matrix'], out['chars'],
        out['protos'])):

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
            print('[*] '+pr[0], minw, idx)
        else:
            smatch += 1
            print('[+] '+pr[0]+' / '+','.join(minc), minw, idx)
    else:
        mmatch += 1
        print('[!] '+pr[0] + ' / '+','.join(minc), minw, idx)

print(match, smatch, mmatch)
input('next run')

txt = '<ul>'
for idx,(p,m,c,pr) in enumerate(zip(out['patterns'], out['matrix'], out['chars'],
        out['protos'])):

    w,p,r = sankoff_parsimony(
            p,
            out['taxa'],
            tree,
            m,
            c
            )

    labels = {}
    for t,c in p[0]:
        labels[tree[t]['label']] = '<b>['+c+']</b> <sup>'+tree[t]['label']+'</sup>'

    tree.output('html', filename='patterns/pattern_'+str(idx+1), labels=labels)
    txt += '<li><a href="patterns/pattern_{0}.html">Proto-Tukano *{1} ({1})</a></li>\n'.format(idx+1,
            pr[0])

with open('navi.html', 'w') as f:
    f.write('<html><head><meta charset="utf-8"</meta></head><body>'+txt+'</body></html>')

input('next run')
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
        #lower_bound = 200,
        #guide_tree = tree.newick, #"((((((((((Bar,Tat),(Kub,Tan)),Kar),Pis),Tuy,Yur),(Pir,(Tuk,Wan))),((Des,Yup),Sir)),((Kor,(Sek,Sio)),Kue)),Mak),Bas);", #tree.newick,
        iterations = 50000,
        sample_steps = 100,
        )


with open('tukano.trees','w') as f:
    for t in trees[0]:
        f.write(t+'\n')
    f.write(str(trees[1]))
        

#import networkx as nx
#taxa = out['taxa']
#matrices = out['matrix']
#chars = out['chars']
#patterns = out['patterns']
#
#G = nx.Graph()
#for i,tA in enumerate(taxa):
#    for j,tB in enumerate(taxa):
#        if i < j:
#
#            # get score between two taxa
#            all_scores = []
#            for pattern,matrix,charset in zip(patterns, matrices, chars):
#
#                patternA = pattern[i]
#                patternB = pattern[j]
#
#                # calculate minimal weight
#                scores = []
#                for pA in patternA:
#                    pAidx = charset.index(pA)
#                    for pB in patternB:
#                        pBidx = charset.index(pB)
#
#                        score = matrix[pAidx][pBidx]
#                    scores += [score]
#                all_scores += [min(scores)]
#
#            G.add_edge(tA, tB, weight=sum(all_scores))
#
#g = nx.minimum_spanning_tree(G)
#print(sum([w[2]['weight'] for w in g.edges(data=True)]))

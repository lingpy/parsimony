# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2015-07-14 09:38
# modified : 2015-07-14 09:38
"""
Compute parsimony analyses.

Notes
-----

Code currently computes weighted parsimony. Required as input data are:

    * the patterns (the states in the leaves, passed as a list, multiple values
      are allowed and are interpreted as potential states of which the best
      states are then chosen)
    * taxa (the name of the languages, taxonomic units)
    * tree (the tree as a lingpy.cogent.object)
    * transitions (the matrix defining the transitions among the characters
    * characters: all characters that occur in the process


Ideas
-----

We can use this version of parsimony now in order to test all possible trees,
provided the number of taxonomic units is not too big. For this, we would just
test all patterns.
"""

__author__="Johann-Mattis List"
__date__="2015-07-14"

from lingpy import *
import itertools



def sankoff_parsimony_up(
        patterns, # the patterns in each taxonomic unit
        taxa, # the taxonomic units corresponding to the patterns
        tree, # the reference tree
        transitions, # the transition matrix,
        characters # the characters as they are provided in the transition matrix
        ):

    W = {}

    # get all characters

    # start iteration
    for node in tree.postorder():
        
        # name of node for convenience
        nname = node.Name

        if node.isTip():

            W[nname] = {}
            for char in characters:
                if char in patterns[taxa.index(nname)]:
                    W[nname][char] = 0
                else:
                    W[nname][char] = 1000000
        else:
            
            W[nname] = {}

            # iterate over the states
            for nchar in characters:
                
                nscores = []
                
                # iterate over the children
                for child in node.Children:
                    cname = child.Name

                    scores = []
                    for cchar in characters:
                        
                        # get the weight in child
                        wchild = W[cname][cchar]

                        # get the new weight due to transition process
                        wnew = wchild + transitions[
                                characters.index(nchar)
                                ][characters.index(cchar)]

                        # append to scores
                        scores += [wnew]

                    # get the minimal score for the char
                    smin = min(scores)

                    nscores += [smin]

                W[nname][nchar] = sum(nscores)
    
    return W
                    
def sankoff_parsimony_down(
        weights,
        patterns,
        taxa,
        tree,
        transitions,
        characters
        ):
    
    # get the root
    root = tree.root().Name

    # get the root chars
    smin = min(weights[root].values())

    # get the starting chars
    rchars = [a for a,b in weights[root].items() if b == smin]

    # initialize the queue by appending the children of the root
    queue = []
    scenario = {}
    root_scens = [(root,char,[]) for char in rchars]
    for child in tree.Children:
        for char in rchars:
            scenario[root,char] = set()
            queue += [[child, child.Name, root, char]]
    
    while queue:
        
        # get the children of the root
        node, name, parent, pchar = queue.pop(0)
        
        # get the smalles transition
        scores = []
        for char in characters:
            i = characters.index(pchar)
            j = characters.index(char)
            score = transitions[i][j]
            
            if weights[parent][pchar] - score == weights[name][char]:
                scores += [0]
            else:
                scores += [1]

        # get the min score
        smin = min(scores)

        # get the min-chars
        mchars = [a for a,b in zip(characters,scores) if b == smin]
        #print(pchar,name,'mchars',mchars,list(zip(scores,characters)))

        # append the children
        for child in node.Children:
            for char in mchars:
                scenario[parent,pchar].add((name,char))
                scenario[name,char] = set()
                queue += [
                        (child, child.Name,name, char)
                        ]

        if node.isTip():
            nchars = patterns[taxa.index(name)]
            for nchar in nchars:
                scenario[parent,pchar].add((name,nchar))
                scenario[name,nchar] = set()
    
     
    # now go back and get all scenarios
    queue = root_scens
    outs = []
    while queue:
        node,char,scen = queue.pop()
        scen += [(node,char)]

        # get the next value
        children = scenario[node,char]
        for cnode,cchar in children:
            queue += [(cnode,cchar,scen)]
        
        if len(scen) == len(tree.getNodeNames()):
            outs += [scen]
    
    # reduce duplicates (no idea how they would come up anyway) XXX
    outs = sorted(set([tuple(sorted(scen)) for scen in outs]))
    
    return outs

def sankoff_parsimony(
        patterns,
        taxa,
        tree,
        transitions,
        characters,
        pprint=False
        ):

    W = sankoff_parsimony_up(
            patterns,
            taxa,
            tree,
            transitions,
            characters
            )
    
    # get minimal weight
    smin = min(W[tree.root().Name].values())
    weights = [b for a,b in W[tree.root().Name].items() if b == smin]
    scenarios = sankoff_parsimony_down(
            W,
            patterns,
            taxa,
            tree,
            transitions,
            characters
            )

    if pprint:
        for i,out in enumerate(scenarios):
            tr = tree.asciiArt()
            for k,v in out:
                target = v+len(k) * '-'
                tr = tr.replace(k,target[:len(k)])
            print(tr)
            print(weights[i])
            print('')

    
    return weights, scenarios, W

def ftree(taxlist):
    """
    Function returns all possible trees for a given set of taxa.
    """
    
    # initialize queue
    queue = [[t for t in taxlist]]

    # store trees in this one
    trees = []
    

    while queue:
        
        taxa = queue.pop(0)
        
        if len(taxa) == 1:
            yield taxa
        else:
            # iterate over all pairwise combinations of taxa
            for a,b in itertools.combinations(taxa,2):
                
                # identify the combined elements
                new_taxa = [t for t in taxa if t not in (a,b)]+[(a,b)]

                queue += [new_taxa]
    
    return trees


def all_nodes_of_binary_tree(newick):
    """
    Function returns all nodes of a tree passed as Newick string.
    
    Notes
    -----

    This function employs a simple search algorithm and splits a tree in binary
    manner in pieces right until all nodes are extracted.

    """
    
    if newick.endswith(';'):
        newick = newick[:-1]
        
    out = [newick]
    queue = [newick[1:-1]]

    while queue:
                
        # find un-bracketed part inbetween
        nwk = queue.pop(0)

        brackets = ''
        idxs = [-1]
        for i,k in enumerate(nwk):

            if k == '(':
                brackets += '('
            elif k == ')':
                brackets = brackets[:-1]
            
            if not brackets and k == ',':
                idxs += [i]

        idxs += [i+1]

        for i,idx in enumerate(idxs[1:]):

            nwk_tmp = nwk[idxs[i]+1:idx]
            out += [nwk_tmp]

            if nwk_tmp.startswith('('):
                nwk_tmp = nwk_tmp[1:-1]
            if ',' in nwk_tmp:
                queue += [nwk_tmp]

                
            #    nwkA = nwk[:i]
            #    nwkB = nwk[i+1:]
            #    out += [nwkA, nwkB]

            #    if nwkA.startswith('('):
            #        nwkA = nwkA[1:-1]
            #    if nwkB.startswith('('):
            #        nwkB = nwkB[1:-1]

            #    if ',' in nwkA:
            #        queue += [nwkA]
            #    if ',' in nwkB:
            #        queue += [nwkB]
            #    break

    return out

def all_rooted_binary_trees(*taxa):
    """
    Compute all rooted trees.

    Notes
    -----

    This procedure yields all rooted binary trees for a given set of taxa, as
    described in :bib:`Felsenstein1978`.
    """

    if len(taxa) < 2:
        yield '('+','.join(taxa)+');'

    # construct the queue
    queue = [('('+','.join(taxa[:2])+')',taxa[:2]+tuple(['('+','.join(taxa[:2])+')'],), taxa[2:])]

    while queue:
        
        # get the three items
        tree, visited, tovisit = queue.pop(0)
        
        # get the new taxon and the new list of nodes to be visited
        new_taxon = tovisit[0]
        new_tovisit = tovisit[1:]

        # add the new taxon to all internal points in the tree
        for v in visited:
            new_tree = tree.replace(v, '('+v+','+new_taxon+')')
            
            # get all nodes of a binary tree to attach new taxa
            # XXX note that we may save time here by not computing all internal
            # nodes in each run, but rather by adding new nodes in case they
            # have not yet been added. The problem is, however, that there may
            # be nodes which are no longer nodes in the tree, since introducing
            # new nodes may split them up, I have not yet an idea how to catch
            # these cases on the run
            new_visited = all_nodes_of_binary_tree(new_tree)

            if not new_tovisit:
                yield new_tree+';'
            else:
                queue += [(new_tree, new_visited, new_tovisit)]




patterns = [['b'],['c'],['a'],['a'],['b']]
taxa = ['A','B','C','D','E']
tree = Tree('(((A,B),(C,D)),E);')
tree2 = Tree('((A,E),(B,C),D);')
transitions = [
         # a b c ø
        [0, 2, 4, 1], # a
        [2, 0, 6, 0], # b
        [4, 6, 0, 10], # c
        [1, 0, 10, 0], # ø
        ]
characters = ['a','b','c',"ø"]

weights, scenarios, W = sankoff_parsimony(patterns, taxa, tree, transitions,
        characters, pprint=True)
weights, scenarios, W = sankoff_parsimony(patterns, taxa, tree2, transitions,
        characters, pprint=True)


#wmins = {}
#for t in ftree(taxa):
#    
#    t = str(t[0])+';'
#    w = sankoff_parsimony_up(
#            patterns,
#            taxa,
#            Tree(t),
#            transitions,
#            characters
#            )
#
#    wmin = min(w['root'].values())
#    
#    wmins[t] = wmin
#
#mint = min(wmins.values())
#gt = [(a,b) for a,b in wmins.items() if b == mint]



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

#wmins = {}
#for t in all_rooted_binary_trees(*taxa):
#    
#    w = sankoff_parsimony_up(
#            patterns,
#            taxa,
#            Tree(t),
#            matrix,
#            chars
#            )
#
#    wmin = min(w['root'].values())
#    
#    wmins[t] = wmin
#
#mint = min(wmins.values())
#gt = [(a,b) for a,b in wmins.items() if b == mint]





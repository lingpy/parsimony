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
import newick as nwk

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
    for node in tree.postorder:
        
        # name of node for convenience
        nname = node

        if tree[node]['leave']:

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
                for child in tree[node]['children']:
                    cname = child

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
    root = tree.root

    # get the root chars
    smin = min(weights[root].values())

    # get the starting chars
    rchars = [a for a,b in weights[root].items() if b == smin]
    
    # prepare the queue
    queue = []
    for char in rchars:
        nodes = []
        for child in tree[tree.root]['children']:
            nodes += [child]
        queue += [([(nodes, tree.root, char)], [(tree.root, char)])]
    
    # prepare the scenarios which are written to output
    outs = []
    
    # start the loop
    while queue:

        nodes, scenario = queue.pop(0)

        if not nodes:
            outs += [scenario]
        else:
            # get children and parent
            children, parent, pchar = nodes.pop()
            pidx = characters.index(pchar)

            # get the best scoring combination for scenario and children
            pscore = weights[parent][pchar]

            combs = itertools.product(*len(children) * [characters])
            
            for comb in combs:
                score = 0
                for i,char in enumerate(comb):
                    cidx = characters.index(char)
                    score += transitions[pidx][cidx]
                    score += weights[children[i]][char]
                
                if score == pscore:
                    new_nodes = [n for n in nodes]
                    new_scenario = [s for s in scenario]
                    
                    for child,char in zip(children,comb):
                        new_nodes += [(tree[child]['children'], child, char)]
                        new_scenario += [(child, char)]

                    queue += [(new_nodes, new_scenario)]
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
    smin = min(W[tree.root].values())
    weights = [b for a,b in W[tree.root].items() if b == smin]
    scenarios = sankoff_parsimony_down(
            W,
            patterns,
            taxa,
            tree,
            transitions,
            characters
            )

    if pprint:
        tmp_tree = Tree(tree.newick)
        C = {}
        for k,v in tmp_tree.getNodesDict().items():
            C[str(v)[:-1]] = k

        for i,out in enumerate(scenarios):
            tr = tmp_tree.asciiArt()
            for k,v in out:
                target = v+len(C[k]) * '-'
                
                # get the nodes dict
                tr = tr.replace(C[k], target[:len(C[k])])
            print(tr)
            print(smin)
            print('')

    
    return weights, scenarios, W

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
            new_visited = [x[0] for x in
                    nwk.all_nodes_of_newick_tree(new_tree)]+[new_tree]

            if not new_tovisit:
                yield new_tree+';'
            else:
                queue += [(new_tree, new_visited, new_tovisit)]

def best_tree_brute_force(
        patterns,
        taxa,
        transitions,
        characters,
        proto_forms=False,
        verbose=False
        ):
    """
    This is an experimental parsimony version that allows for ordered
    character states.
    """

    minScore = 1000000000
    bestTree = []

    for idx,tree in enumerate(all_rooted_binary_trees(*taxa)):
        t = nwk.LingPyTree(tree)
        if verbose:
            print('[{0}] {1}...'.format(idx+1, t.newick))

        score = 0
        for i,(p,m,c) in enumerate(zip(patterns, transitions, characters)):
            weights = sankoff_parsimony_up(
                    p,
                    taxa,
                    t,
                    m,
                    c
                    )
            if not proto_forms:
                minWeight = min(weights[t.root].values())
            else:
                minWeight = weights[t.root][proto_forms[i]]
                
            score += minWeight
            
            if score > minScore:
                break

        if score == minScore:
            bestTree += [(t.newick,score)]
        elif score < minScore:
            minScore = score
            bestTree = [(t.newick,score)]

    return bestTree, minScore

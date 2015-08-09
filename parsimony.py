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
import lingpy
import itertools
import newick as nwk
import networkx as nx
import random

def sankoff_parsimony_up(
        patterns, # the patterns in each taxonomic unit
        taxa, # the taxonomic units corresponding to the patterns
        tree, # the reference tree tree nodes in post-order
        transitions, # the transition matrix,
        characters, # the characters as they are provided in the transition matrix
        weight_only = False, # specify wether only weights should be returned
        weight_and_chars = False
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
    
    if weight_only:
        return min(W[tree.root].values())
    if weight_and_chars:
        minw = min(W[tree.root].values())
        minchars = [x for x,y in W[tree.root].items() if y == minw]
        return minw,minchars
    
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


def swap_tree(tree):
    
    # make safe tree
    if not '"' in tree:
        tree = nwk.safe_newick_string(tree)

    # swap two nodes of the tree
    nodes = list(nwk.nodes_in_tree(tree))[1:]
    random.shuffle(nodes)
    
    # choose two nodes to be swapped
    nodeA = nodes.pop(0)

    # get another node that can be interchanged
    while nodes:
        nodeB = nodes.pop(0)
        if nodeB in nodeA or nodeA in nodeB:
            pass
        else:
            break

    tree = tree.replace(nodeA+',', '#dummyA#,')
    tree = tree.replace(nodeA+')', '#dummyA#)')
    tree = tree.replace(nodeB+',', '#dummyB#,')
    tree = tree.replace(nodeB+')', '#dummyB#)')

    tree = tree.replace('#dummyA#', nodeB)
    tree = tree.replace('#dummyB#', nodeA)

    return nwk.sort_tree(tree)

def mst_weight(
        taxa,
        patterns,
        matrices,
        characters
        ):
    """
    Calculate minimal weight of unsorted trees.
    """

    G = nx.Graph()
    for i,tA in enumerate(taxa):
        for j,tB in enumerate(taxa):
            if i < j:
                all_scores = []
                for pt,mt,cs in zip(patterns, matrices, characters):
                    ptA = pt[i]
                    ptB = pt[j]
                    scores = []
                    for pA in ptA:
                        idxA = cs.index(pA)
                        for pB in ptB:
                            idxB = cs.index(pB)
                            score = mt[idxA][idxB]
                        scores += [score]
                    all_scores += [min(scores)]
                G.add_edge(tA, tB, weight=sum(all_scores))
    g = nx.minimum_spanning_tree(G)
    
    return sum([w[2]['weight'] for w in g.edges(data=True)]) / 2

def heuristic_parsimony(
        taxa, 
        patterns,
        transitions,
        characters,
        guide_tree = False,
        verbose = True,
        lower_bound = False,
        iterations = 300,
        sample_steps = 100,
        log = False,
        ):
    """
    Try to make a branch and bound parsimony calculation.
    """
    
    if not guide_tree:
        guide_tree = lingpy.basic.tree.random_tree(taxa)


    lower_bound = 0
    ltree = nwk.LingPyTree(guide_tree)
    for p,t,c in zip(patterns, transitions, characters):
        lower_bound += sankoff_parsimony_up(
                p,
                taxa,
                ltree,
                t,
                c,
                weight_only=True
            )
    print("[i] Lower Bound (in guide tree):",lower_bound)

    # we start doing the same as in the case of the calculation of all rooted
    # trees below
    if len(taxa) <= 2:
        return '('+','.join(taxa)+');'

    # make queue with taxa included and taxa to be visited
    tree = nwk.sort_tree(guide_tree)
    queue = [(tree, lower_bound)]
    visited = [tree]
    trees = [tree]

    # create a generator for all rooted binary trees
    gen = all_rooted_binary_trees(*taxa)
    previous = 0

    while queue:       

        # modify queue
        queue = sorted(queue, key=lambda x: x[1])

        # check whether tree is in data or not
        #if tree in visited:
            
        # try creating a new tree in three steps:
        # a) swap the tree
        # b) make a random tree
        # c) take a generated tree

        # add next taxon
        tree, bound = queue.pop(0)

        if tree.endswith(';'):
            tree = tree[:-1]
        
        forest = []

        # try and get the derivations from the best trees
        for i in range(4 * len(taxa)):
            new_tree = swap_tree(random.choice(trees))
            if new_tree not in visited:
                forest += [new_tree]
                visited += [new_tree]
            if previous < len(visited) and len(visited) % sample_steps == 0:
                print("[i] Investigated {0} trees so far, currently holding {1} trees with best score of {2}.".format(len(visited), len(trees), lower_bound)) 
                previous = len(visited)
        
        for i in range(len(taxa)):
            new_tree = swap_tree(tree)
            if new_tree not in visited:
                forest += [new_tree]
                visited += [new_tree]
            if previous < len(visited) and len(visited) % sample_steps == 0:
                print("[j] Investigated {0} trees so far, currently holding {1} trees with best score of {2}.".format(len(visited), len(trees), lower_bound)) 
                previous = len(visited)    

        # go on with b
        for i in range(len(taxa) // 2):
            new_tree = nwk.sort_tree(lingpy.basic.tree.random_tree(taxa))
            if new_tree not in visited:
                forest += [new_tree]
                visited += [new_tree]
            if previous < len(visited) and len(visited) % sample_steps == 0:
                print("[k] Investigated {0} trees so far, currently holding {1} trees with best score of {2}.".format(len(visited), len(trees), lower_bound)) 
                previous = len(visited)

        for i in range(1 * len(taxa) // 4):
            new_tree = nwk.sort_tree(next(gen))
            if new_tree not in visited:
                forest += [new_tree]
                visited += [new_tree]
            if previous < len(visited) and len(visited) % sample_steps == 0:
                print("[l] Investigated {0} trees so far, currently holding {1} trees with best score of {2}.".format(len(visited), len(trees), lower_bound)) 
                previous = len(visited)

        best_scores = []
        for tree in forest:
            score = 0
            lp_tree = nwk.LingPyTree(tree)
            for p,t,c in zip(patterns, transitions, characters):
                weight,chars  = sankoff_parsimony_up(
                        p,
                        taxa,
                        lp_tree,
                        t,
                        c,
                        weight_and_chars =True
                        )
                score += weight

            # append stuff to queue
            best_scores += [(tree, score)]

        for k in sorted(best_scores, key=lambda x:
                x[1])[:len(best_scores) // 4]:
            queue += [(k[0], k[1])]
            
            if k[1] < lower_bound:
                trees = [tree]
                lower_bound = k[1]
            elif k[1] == lower_bound:
                trees += [tree]

        if len(visited) > iterations:
            answer = input("[?] Number of chosen iterations is reached, do you want to go on with the analysis? y/n ")
            if answer == 'y':
                while True:
                    number = input("[?] How many iterations? ")
                    try:
                        number = int(number)
                        iterations += number
                        break
                    except:
                        pass
            else:
                break

    return trees, lower_bound

def branch_and_bound(
        taxa, 
        patterns,
        transitions,
        characters,
        guide_tree = False,
        verbose = True,
        lower_bound = False,
        sample_steps = 100,
        ):
    """
    Try to make a branch and bound parsimony calculation.
    """
    
    # calculate the lower bound
    if lower_bound:
        pass
    elif not guide_tree:
        lower_bound = mst_weight(taxa, patterns, transitions, characters) * 2

        print("[i] Lower Bound (estimated):", lower_bound)
    else:
        lower_bound = 0
        ltree = nwk.LingPyTree(guide_tree)
        for p,t,c in zip(patterns, transitions, characters):
            lower_bound += sankoff_parsimony_up(
                    p,
                    taxa,
                    ltree,
                    t,
                    c,
                    weight_only=True
                )
        print("[i] Lower Bound (in guide tree):",lower_bound)

    trees = []

    # we start doing the same as in the case of the calculation of all rooted
    # trees below
    if len(taxa) <= 2:
        return '('+','.join(taxa)+');'

    # make queue with taxa included and taxa to be visited
    queue = [('('+','.join(taxa[:2])+')', taxa[2:], lower_bound)]

    visited = 0
    all_trees = []
    previous = 0

    M = {}

    while queue:
        queue = sorted(queue, key = lambda x: (len(x[1]), x[2]))
        
        # add next taxon
        tree, rest, bound = queue.pop(0)
        
        if rest:
            next_taxon = rest.pop()
            nodes = list(nwk.nodes_in_tree(tree))
            random.shuffle(nodes)
            for node in nwk.nodes_in_tree(tree):
                new_tree = tree.replace(node, '('+next_taxon+','+node+')')
                visited += 1

                # parsimony evaluation and lower bound comes here
                score = 0
                ltree = nwk.LingPyTree(new_tree)
                for p,t,c in zip(patterns, transitions, characters):

                    weight,chars  = sankoff_parsimony_up(
                            p,
                            taxa,
                            ltree,
                            t,
                            c,
                            weight_and_chars =True
                            )
                    score += weight

                if rest:
                    if tuple(rest) in M:
                        mst = M[tuple(rest)]
                    else:
                        mst = mst_weight(rest, patterns, transitions, characters)
                        M[tuple(rest)] = mst
                    score += mst
                
                else:
                    all_trees += [new_tree]

                if score <= lower_bound:
                    r = [x for x in rest]
                    random.shuffle(r)
                    queue += [(new_tree, r, score)]
                    
                    if not rest:
                        if score < lower_bound:
                            lower_bound = score
                            trees = [nwk.sort_tree(new_tree)]
                            queue = [q for q in queue if q[2] <= lower_bound]
                        else:
                            trees += [nwk.sort_tree(new_tree)]
                if len(all_trees) % sample_steps == 0 and len(all_trees) > previous:
                    print("[i] Checked {0} trees so far, current lower bound is {1} with {2} best trees.".format(
                        len(all_trees),
                        lower_bound,
                        len(trees)
                        ))
                    previous = len(all_trees)

    return trees, lower_bound

def all_rooted_binary_trees(*taxa):
    """
    Compute all rooted trees.

    Notes
    -----

    This procedure yields all rooted binary trees for a given set of taxa, as
    described in :bib:`Felsenstein1978`. It implements a depth-first search.
    """
    if len(taxa) <= 2:
        yield '('+','.join(taxa)+');'

    # make queue with taxa included and taxa to be visited
    queue = [('('+','.join(taxa[:2])+')', list(taxa[2:]))]

    out = []

    while queue:
        
        # add next taxon
        tree, rest = queue.pop()

        if rest:
            next_taxon = rest.pop()
            
            nodes = list(nwk.nodes_in_tree(tree))
            random.shuffle(nodes)
            for node in nodes: 
                new_tree = tree.replace(node, '('+next_taxon+','+node+')')
                
                r = [x for x in rest]
                random.shuffle(r)
                queue += [(new_tree, r)]
                if not rest:
                    yield new_tree

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



#def nodes_in_tree(tree):
#    """
#    This methods yields all nodes in a tree in newick format without labels for
#    nodes and branch lengths being assigned.
#    """
#
#    stack = [tree[1:-1]]
#    out = [tree]
#
#    yield tree
#    
#    while stack:
#        tmp = stack.pop()
#        brackets = 0
#
#        idx = 0
#        for i,c in enumerate(tmp):
#            if c == '(':
#                brackets += 1
#            elif c == ')':
#                brackets -= 1
#
#            if not brackets and c == ',':
#                tree = tmp[idx:i]
#                idx = i+1
#
#                yield tree
#
#                if tree.startswith('('):
#                    tree = tree[1:-1]
#
#                if ',' in tree:
#                    stack += [tree]
#        
#        tree = tmp[idx:]
#        if tree.startswith('('):
#            tree = tree[1:-1]
#        
#        yield tree
#
#        if ',' in tree:
#            stack += [tree]
#
#    return out


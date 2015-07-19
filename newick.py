# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2015-07-19 12:00
# modified : 2015-07-19 12:00
"""
Newick Module of LingPy for parsing of Newick strings.
"""

__author__="Johann-Mattis List"
__date__="2015-07-19"


import json

def clean_newick_string(newick):
    """
    Helper function to reduce all branch-lengths from a newick string.
    """
    start = newick
    out = ''
    
    while start:
        idxA = start.find(':') 
        idxB = start.find(')')

        colon = False

        if idxA != -1 and idxB != -1:
            if idxA < idxB:
                idx = idxA
                colon = True
            else:
                idx = idxB
        elif idxA != -1:
            idx = idxA
            colon = True
        elif idxB != -1:
            idx = idxB
        else:
            return out + start
        
        if colon:
            out += start[:idx]
            start = start[idx+1:]
            while start and (start[0].isdigit() or start[0] == '.'):
                start = start[1:]
        else:
            out += start[:idx+1]
            start = start[idx+1:]
            while start and start[0] not in ',;)':
                start = start[1:]
    
    out += start

    return out

def parse_newick(newick):
    """
    Function parses a newick tree to json format.

    Notes
    -----
    The format is a dictionary with sufficient information to further parse the
    tree, and also to use it as input for d3 and other javascript libraries.
    """
    
    D = {}

    if newick.endswith(';'):
        newick = newick[:-1]
    
    nwk, label, blen = label_and_blen(newick)

    root = '('+clean_newick_string(nwk)+')'
    
    D['root'] = root
    D['nodes'] = [root]
    D['leaves'] = []

    D[root] = dict(
            children=[], 
            branch_length = blen, 
            root=True,
            leave=False,
            label = label
            )

    for node, label, blen, parent in all_nodes_of_newick_tree(newick):
        D['nodes'] += [node]

        D[node] = dict(parent=parent, children=[], branch_length=blen, root=False)

        if ',' in node:
            D[node]['leave'] = False
        else:
            D[node]['leave'] = True
            D['leaves'] += [node]
        
        D[parent]['children'] += [node]

    return D

def label_and_blen(nwk):
    """
    Helper function parses a Newick string and returns the highest-order label
    and branch-length.
    """
    idx = nwk[::-1].find(')')

    # no brackets means we are dealing with a leave node
    if idx == -1:
        nwk_base = nwk
        idx = nwk[::-1].find(':')
        if idx == -1:
            return nwk, nwk, '0'
        else:
            nwk_base = nwk[:-idx-1]
            return nwk_base, nwk_base, nwk[-idx:]
    
    # if index is 0, there's no label and no blen
    elif idx == 0:
        return nwk[1:-1], '', '0'

    # else, we carry on
    nwk_base = nwk[:-idx]
    label = nwk[-idx:]

    idx = label[::-1].find(':')
    if idx == -1:
        label = label
        blen = 0
    else:
        blen = label[-idx:]
        label = label[:-idx-1]
    
    return nwk_base[1:-1], label, blen

def all_nodes_of_newick_tree(newick, parent_nodes=False):
    """
    Function returns all nodes of a tree passed as Newick string.
    
    Notes
    -----

    This function employs a simple search algorithm and splits a tree in binary
    manner in pieces right until all nodes are extracted.

    """
    # look for bracket already, don't assume they are around the tree! 
    if newick.endswith(';'):
        newick = newick[:-1]

    # fill the queue
    nwk, label, blen = label_and_blen(newick)
    queue = [(nwk, label, blen)]
    
    # raise error if the number of brackets doesn't fit
    nr_opening_brackets = newick.count('(')
    nr_closing_brackets = newick.count(')')
    if nr_opening_brackets != nr_closing_brackets:
        raise ValueError("The number of brackets is wrong!")
    
    while queue:
                
        # find un-bracketed part inbetween
        nwk, label, blen = queue.pop(0)

        brackets = 0
        idxs = [-1]
        for i,k in enumerate(nwk):

            if k == '(':
                brackets += 1
            elif k == ')':
                brackets -= 1
            
            if not brackets and k == ',':
                idxs += [i]

        idxs += [i+1]

        for i,idx in enumerate(idxs[1:]):

            nwk_tmp = nwk[idxs[i]+1:idx]
            nwk_tmp, label, blen = label_and_blen(nwk_tmp)

            nnwk = clean_newick_string(nwk_tmp)
            npar = clean_newick_string(nwk)

            nnwk = '('+nnwk+')' if ',' in nnwk else nnwk
            npar = '('+npar+')' if ',' in npar else npar
            
            yield nnwk, label, blen, npar

            if ',' in nwk_tmp:
                
                queue += [(nwk_tmp, label, blen)]


def postorder(tree):
    """
    Carry out a post-order traversal of the tree.
    """
    
    # make the stack
    stack = [tree['root']]

    # make the output
    out = []

    # make copy of tree
    ctree = dict([(k,tree[k]['children']) for k in tree['nodes']])

    # climb down the tree

    while stack:
        
        node = stack[-1]
        children = ctree[node]
        
        # if we are at a leave-node, we remove the item from the stack 
        if not children:
            stack.pop()
            out += [node]
            if stack:
                ctree[stack[-1]].pop(0)

        else:
            stack += [children[0]]
    
    return out

class LingpyTree(object):

    def __init__(self, newick):

        self._dict = parse_newick(newick)
        self.root = self._dict['root']
        self.nodes = self._dict['nodes']
        self.leaves = self._dict['leaves']
        self.preorder = self.nodes
        self.postorder = postorder(self._dict)

    def output(self, dtype, filename=None):
        
        if dtype == 'json':
            if filename:
                with open(filename+'.'+dtype, 'w') as f:
                    f.write(json.dumps(self._dict, indent=2))
            else:
                return json.dumps(self._dict, indent=2)

    def __getitem__(self, idx):

        try:
            return self._dict[idx]
        except:
            raise KeyError(idx)

a = LingpyTree('((a:1,b:2)ab:3,(c:1,d:2)cd:4)abcd;')

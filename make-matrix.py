# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2015-08-11 11:09
# modified : 2015-08-11 11:09
"""
<++>
"""

__author__="Johann-Mattis List"
__date__="2015-08-11"

import json

with open('tukano.simple.json') as f:
    data = json.loads(f.read())

out = ''
for i,(p,refs) in enumerate(zip(data['protos'], data['chars'])):

    out += str(i+1)+'\t'+p[0]+'\t' + p[1]+'\t'+', '.join(refs)+'\n'

with open('tukano.data.tsv', 'w') as f:
    f.write(out)


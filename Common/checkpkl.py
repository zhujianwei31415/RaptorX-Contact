#！/usr/bin/env python3
import numpy as np
import sys
import pickle

if len(sys.argv) != 3:
    sys.exit('Usage: %s <py2pkl> <py3pkl>' % sys.argv[0])
py2pkl, py3pkl = sys.argv[1:3]

hhm1 = pickle.load( open(py2pkl, 'rb'), encoding='latin' )
hhm2 = pickle.load( open(py3pkl, 'rb') )

assert hhm1.keys() == hhm2.keys()
for k, v in hhm1.items():
    if isinstance( v, (str, list, np.int32, np.float32) ):
        if v != hhm2[k]:
            print(hhm1['name'], hhm2['name'], False)
            sys.exit()
    else:
        if not np.array_equal(v, hhm2[k]):
            print(hhm1['name'], hhm2['name'], False)
            sys.exit()
print(hhm1['name'], hhm2['name'], True)

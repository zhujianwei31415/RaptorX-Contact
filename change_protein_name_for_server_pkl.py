#!/usr/bin/env python2
import sys

try:
    import cPickle as pickle
except:
    import pickle

if len(sys.argv) != 4:
    sys.exit('Usage: %s <target_name> <server_pkl> <output_pkl>' % sys.argv[0])
tarname, serverpkl, outpkl = sys.argv[1:4]

with open(serverpkl, 'rb') as fin:
    samples = pickle.load(fin)
feats = []
for s in samples:
    s['name'] = tarname
    feats.append(s)
with open(outpkl, 'wb') as fout:
    pickle.dump(feats, fout)
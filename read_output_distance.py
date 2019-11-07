from __future__ import absolute_import, division, print_function

import numpy as np
import sys

try:
    import cPickle as pickle
except:
    import pickle


def convert_25_to_12_bin(prob_bin25):
    prob_bin12 = np.array([0. for _ in range(12)])
    for i in range(24):
        prob_bin12[i // 2] += prob_bin25[i]
    prob_bin12[-1] += prob_bin25[-1]
    return prob_bin12


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit('Usage: %s <XXX.pkl> <XXX.epad_prob> <XXX.epad_prob_25> <XXXXX.con_mat>' % sys.argv[0])
    pkl_file, epad_file, epad_25_file, con_file = sys.argv[1:5]
  
    # read in pkl file
    sample = None
    with open(pkl_file, 'r') as fin:
        sample = pickle.load(fin)
    assert isinstance(sample, tuple) and 6 == len(sample), (
        'ERROR: wrong file format %s' % pkl_file)
    assert isinstance(sample[2], dict) and 'CbCb_Discrete25C' in sample[2], (
        'ERROR: wrong dict in file %s' % pkl_file)
    seq, dist, con = sample[1], sample[2]['CbCb_Discrete25C'], sample[3]['CbCb']
    assert (len(seq), len(seq), 25) == dist.shape and \
        (len(seq), len(seq)) == con.shape, \
        ('ERROR: wrong sequence length %s' % pkl_file)

    # write contact matrix
    np.savetxt(con_file, con)

    # output epad for 25 bin and 12 bin
    f12h = open(epad_file, 'w')
    f25h = open(epad_25_file, 'w')
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            # write 25 bin
            assert np.isclose(con[i, j], np.sum(dist[i,j,:8]))
            print(i, j, ' '.join(['%.6f'%(_*100.) for _ in dist[i,j]]), file=f25h)
            # convert 25 bins to 12 bins
            tmp12bin = convert_25_to_12_bin(dist[i, j])
            # write 12 bin
            assert np.isclose(con[i, j], np.sum(tmp12bin[:4]))
            print(i, j, ' '.join(['%.6f'%(_*100.) for _ in tmp12bin]), file=f12h)
    f12h.close()
    f25h.close()
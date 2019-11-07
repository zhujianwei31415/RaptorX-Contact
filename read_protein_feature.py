from __future__ import absolute_import, division, print_function

import numpy as np
import os
import sys

try:
    import cPickle as pickle
except:
    import pickle


##in the future we may add some code to check correctness
def Load1DFeatureFromFile(file, seqName=None, seq=None):
    # check input feature file
    if not os.path.isfile(file):
        sys.exit('ERROR: cannot find the 1D feature file %s' % file)
    # read in data from file
    allprobs, AAs = [], []
    with open(file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0 or line[0] == '#': continue
            array = line.split()
            assert len(array) == 4 or len(array) == 6 or len(array) == 11, (
                'ERROR: wrong file format for predicted 1D feature %s' % file)
            probs = [ np.float32(_) for _ in array[3:] ]
            allprobs.append(probs)
            AAs.append(array[1])
    AAs = ''.join(AAs)
    # check length
    if seq is not None:
        assert len(seq) == len(allprobs), 'ERROR: wrong length for %s' % file
        assert seq == AAs, 'ERROR: wrong sequence for %s' % file
    # convert probs to numpy array and check
    probs = np.array(allprobs)
    if np.isnan( np.sum(probs) ):
        sys.exit('ERROR: there are NaNs in file %s' % file)

    return probs


def LoadProfile(file, seqName=None, seq=None):
    # check input profile
    if not os.path.isfile(file):
        sys.exit('ERROR: cannot find temporary profile file %s' % file)
    # read in data from file
    fh = open(file, 'r')
    content = list(fh)
    fh.close()
    # check data format
    seqLen = int(content[1].strip())
    assert len(content) == (3 + 2*seqLen), 'ERROR: wrong length for %s' % file
    if seqName is not None:
        assert content[0].strip() == seqName, 'ERROR: wrong name for %s' % file
    if seq is not None:
        assert content[2].strip() == seq, 'ERROR: wrong sequence for %s' % file
    # convert PSFM to float
    allprobs = []
    for line in content[3 : 3+seqLen]:
        probs = [ np.float32(_) for _ in line.split(',') ]
        allprobs.append(probs)
    # convert PSSM to float
    allscores = []   
    for line in content[3+seqLen : 3+2*seqLen ]:
        scores = [ np.float32(_) for _ in line.split(',') ]
        allscores.append(scores)
    # convert list ot numpy array
    probs, scores = np.array(allprobs), np.array(allscores)
    if np.isnan( np.sum(probs) ) or np.isnan( np.sum(scores) ):
        sys.exit('ERROR: there are NaNs in file %s' % file)

    return probs, scores


def LoadECMatrix(file, seqName=None, seq=None):
    # check input profile
    if not os.path.isfile(file):
        sys.exit('ERROR: cannot find EC matrix file %s' % file)
    # read in data from file
    allECs = []
    with open(file, 'r') as fin:
        for line in fin:
            ECs = [ np.float16(_) for _ in line.split() ]
            if seq is not None:
                assert len(seq) == len(ECs), 'ERROR: wrong columns %s' % file
            allECs.append(ECs)
    if seq is not None:
        assert len(seq) == len(allECs), 'ERROR: wrong rows in %s' % file
    # convert list to numpy array
    ECMatrix = np.array(allECs)
    if np.isnan( np.sum(ECMatrix.astype(np.float32)) ):
        sys.exit('ERROR:, there is at least one NaN in %s' % file)

    return ECMatrix


def LoadOtherPairFeatures(file, seqName=None, seq=None):
    # check input potential file
    if not os.path.isfile(file):
        sys.exit('ERROR: cannot find contact potential file %s' % file)
    if seq is None:
        sys.exit('Please provide sequence for %s' % seqName)
    # read in data from file
    indexList = []
    valueList = []
    with open(file, 'r') as fin:
        for line in fin:
            fields = line.split()
            assert len(fields) == 5, 'ERROR: wrong data format for %s' % file
            indexList.append( [ np.int16(_)-1 for _ in fields[:2] ] )
            valueList.append( [ np.float16(_) for _ in fields[2:] ] )
    # check matrix index
    indexArr = np.transpose( np.array(indexList) )
    if np.amin(indexArr) < 0 or np.amax(indexArr) >= len(seq):
        sys.exit('ERROR: index out of sequence length for %s' % file)
    # convert data to numpy array
    allPairs = np.zeros((len(seq),len(seq),len(valueList[0])), dtype=np.float16)
    allPairs[ indexArr[0], indexArr[1] ] = valueList # copy data to half matrix
    allPairs[ indexArr[1], indexArr[0] ] = valueList # make the matrix symmetric
    if np.isnan( np.sum( allPairs.astype(np.float32) ) ):
        sys.exit('ERROR: there are NaNs in file %s' % file)

    return allPairs


def ReadFeatures(p=None, DataSourceDir=None):
    if p is None:
        sys.exit('Please specify a valid target name!')
    if DataSourceDir is None:
        sys.exit('Please specify a folder containing all features for target!')
    if not os.path.isdir(DataSourceDir):
        sys.exit('The feature directory does not exist: ', DataSourceDir)

    OneProtein = dict()
    OneProtein['name'] = p

    # read in sequence from .seq file
    fastafh=open(DataSourceDir + p + ".seq", "r")
    fastacontent = [ line.strip() for line in list(fastafh) ]
    if not fastacontent[0].startswith('>'):
        OneProtein['sequence'] = ''.join(fastacontent)
    else:
        OneProtein['sequence'] = ''.join(fastacontent[1: ])
    fastafh.close()

    # read in 3-state secondary structure from .ss3 file
    OneProtein['SS3'] = Load1DFeatureFromFile(file=DataSourceDir+p+".ss3",
                                              seqName=p,
                                              seq=OneProtein['sequence'])

    # read in 8-state secondary structure from .ss8 file
    OneProtein['SS8'] = Load1DFeatureFromFile(file=DataSourceDir+p+".ss8",
                                              seqName=p,
                                              seq=OneProtein['sequence'])

    # read in solvent accessibility from .acc file
    OneProtein['ACC'] = Load1DFeatureFromFile(file=DataSourceDir+p+".acc",
                                              seqName=p,
                                              seq=OneProtein['sequence'])

    # read in disorder prediction from .diso file
    OneProtein['DISO'] = Load1DFeatureFromFile(file=DataSourceDir+p+".diso",
                                               seqName=p,
                                               seq=OneProtein['sequence'])

    # read in PSFM and PSSM from .tgt file
    import subprocess
    import tempfile
    with tempfile.NamedTemporaryFile(delete=True) as tmp:
        # check input file
        tgtf = DataSourceDir + p + ".tgt"
        if not os.path.isfile(tgtf):
            sys.exit('ERROR: cannot find the following file %s' % tgtf)
        # check helper program
        PrintTGT = os.path.abspath(os.path.dirname(__file__)) + "/PrintTGT"
        if not os.path.isfile(PrintTGT):
            sys.exit('ERROR: cannot find the helper program: %s' % PrintTGT)
        # write temporary profile file
        with open(tmp.name, 'w') as fout:
            cmdStr = [PrintTGT, p, DataSourceDir]
            proc = subprocess.Popen(cmdStr, stdout=fout, stderr=fout)
            proc.wait()
            fout.seek(0)
        # read in PSFM and PSSM from temporary profile
        OneProtein['PSFM'], OneProtein['PSSM'] = LoadProfile(
            file=tmp.name,
            seqName=p,
            seq=OneProtein['sequence'])
    
    # read in evolutionary couplings from .ccmZ file
    OneProtein['ccmpredZ'] = LoadECMatrix(file=DataSourceDir+p+'.ccmZ',
                                          seqName=p,
                                          seq=OneProtein['sequence'])

    # read in contact potential from .pot file
    OneProtein['OtherPairs'] = LoadOtherPairFeatures(file=DataSourceDir+p+'.pot',
                                                     seqName=p,
                                                     seq=OneProtein['sequence'])

    return OneProtein


def main(listFile, featureMetaDir, outPickle):
    print('protein list file=', listFile)
    print('feature directory=', featureMetaDir)
    print('output picklefile=', outPickle)

    # check file and directory
    if not os.path.isdir(featureMetaDir):
        sys.exit('the provided feature folder is invalid: %s' % featureMetaDir)
    if not os.path.isfile(listFile):
        sys.exit('the provided protein list file is invalid: %s' % listFile)
    
    # read in proteins
    with open(listFile, 'r') as fin:
        proteins = [_.strip() for _ in fin]

    # read in features
    pFeatures = []
    for p in proteins:
        thisFeatureDir = featureMetaDir + '/'
        pFeature = ReadFeatures( p=p, DataSourceDir=thisFeatureDir )
        pFeatures.append(pFeature)

        if len(pFeatures) % 500 == 0:
            print('finished loading features for ', len(pFeatures), ' proteins')

    # writing the result to pickle file
    with open(outPickle, 'wb') as fout:
        pickle.dump( pFeatures, fout )


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit( """Usage: {} <listFile> <featureMetaDir> <outPickle>
    listFile:        the file containing a list of proteins, each protein
                     name in one line
    featureMetaDir:  specify a folder containing all the features, under
                     which each protein has an independent feature folder
                     named after feat_proteinName_contact
    outPickle:       the output feature file in pickle format, in which
                     echo protein is one dict, the keys are ['name',
                     'sequence', 'SS3', 'SS8', 'ACC', 'DISO', 'PSFM',
                     'PSSM', 'ccmpredZ', 'OtherPairs']

    This script only reads protein features for distance prediciton.""".format(
        sys.argv[0]) )
    listFile, featureMetaDir, outPickle = sys.argv[1:4]
    
    ## read a list of protein features into a single PKL file
    main(listFile, featureMetaDir, outPickle)

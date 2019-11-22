#!/bin/bash

if [ $# != 3 ]; then
    echo "Usage: $0 <target_name> <feature_pickle_file> <output_directory>"
    exit 0
fi

# get python
PYTHON="$HOME/anaconda3/envs/raptorx/bin/python"

# get parameters from command line
tarname=$1
featpkl=$2
outdir=$3

# get absolute directory of current script
HOMEDIR=`dirname $(readlink -f "${BASH_SOURCE[0]}")`

# check output directory
outdir=$(readlink -f "$outdir")
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

# check features
echo "Checking feature ..."
if [ ! -f $featpkl ]; then
    echo "Usage: $featpkl does not exist"
    exit 1
fi

# change target protein name for server pickle file
tmpfile=$(mktemp) || exit 1
trap "rm -rf $tmpfile" EXIT
echo "Change server pickle name to $tmpfile"
$PYTHON $HOMEDIR/change_protein_name_for_server_pkl.py $tarname $featpkl $tmpfile

# predict contact matrix
echo "Predicting contact ..."
model="$HOMEDIR/Models/RXContact-DeepModel11410.pkl"
model="$HOMEDIR/Models/RXContact-DeepModel10820.pkl"
distance_file="$outdir/$tarname.predictedDistMatrix.pkl"
$PYTHON $HOMEDIR/DL4DistancePrediction2/RunDistancePredictor2.py -m $model -p $tmpfile -d $outdir
if [ ! -f $distance_file ]; then
    echo "ERROR: distance prediction error $distance_file"
    exit 2
fi

# convert output pickle to text format
echo "Convert pickle file to epad and matrix ..."
epad_file="$outdir/$tarname.epad_prob"
epad_25_file="$outdir/$tarname.epad_prob_25"
mat_file="$outdir/$tarname.contactMatrix.txt"
$PYTHON $HOMEDIR/read_output_distance.py $distance_file $epad_file $epad_25_file $mat_file
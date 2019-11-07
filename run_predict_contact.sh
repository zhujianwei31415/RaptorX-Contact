#!/bin/bash

if [ $# != 3 ]; then
    echo "Usage: $0 <target_name> <sequence_file> <output_directory>"
    exit 0
fi

# get python
PYTHON="$HOME/anaconda3/envs/raptorx/bin/python"

# get parameters from command line
tarname=$1
seqfile=$2
outdir=$3

# get absolute directory of current script
HOMEDIR=`dirname $(readlink -f "${BASH_SOURCE[0]}")`

# check output directory
outdir=$(readlink -f "$outdir")
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

# generate feature
echo "Predicting feature ..."
#$HOMEDIR/Generate_Feature/Generate_Feature.sh $tarname $seqfile $outdir

# read in feature to pkl
echo "Converting feature to pkl file ..."
pkl_file="$outdir/$tarname.contactFeatures.pkl"
tmpfile=$(mktemp)
echo $tarname > $tmpfile
$PYTHON $HOMEDIR/read_protein_feature.py $tmpfile $outdir $pkl_file
rm $tmpfile

# predict contact matrix
echo "Predicting contact ..."
model="$HOMEDIR/Models/RXContact-DeepModel11410.pkl"
model="$HOMEDIR/Models/RXContact-DeepModel10820.pkl"
distance_file="$outdir/$tarname.predictedDistMatrix.pkl"
$PYTHON $HOMEDIR/DL4DistancePrediction2/RunDistancePredictor2.py -m $model -p $pkl_file -d $outdir
if [ ! -f $distance_file ]; then
    echo "ERROR: distance prediction error $distance_file"
    exit 1
fi

# convert output pickle to text format
echo "Convert pickle file to epad and matrix ..."
epad_file="$outdir/$tarname.epad_prob"
epad_25_file="$outdir/$tarname.epad_prob_25"
mat_file="$outdir/$tarname.contactMatrix.txt"
$PYTHON $HOMEDIR/read_output_distance.py $distance_file $epad_file $epad_25_file $mat_file
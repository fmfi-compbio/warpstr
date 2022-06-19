#!/bin/bash
#script calling guppy basecaller for sample path with .fast5 files located in root folder of a disease

if [[ ${#} -le 5 ]]; then
   invoking_script="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")";
   echo "missing arguments when calling "$invoking_script;
   exit 1
fi

SAMPLEPATH=$1;
GUPPY=$2;
FLOWCELL=$3;
KIT=$4;
OUT_SUFFIX=$5;
THREADS=$6;

SAVEPATH=$SAMPLEPATH/"aux"
mkdir -p $SAVEPATH


echo "calling script $GUPPY --input_path $SAMPLEPATH --save_path $SAVEPATH --flowcell $FLOWCELL --kit $KIT --fast5_out --cpu_threads_per_caller $THREADS"
( $GUPPY --input_path $SAMPLEPATH --save_path $SAVEPATH --flowcell $FLOWCELL --kit $KIT --fast5_out --cpu_threads_per_caller $THREADS)

WORKSPACEPATH=$SAVEPATH/"workspace"
ANNOTPATH=$SAMPLEPATH/$OUT_SUFFIX
mkdir -p $ANNOTPATH

echo "find $WORKSPACEPATH -maxdepth 1 -name "*.fast5" -print0 | xargs  -0 cp -t $ANNOTPATH;"
find $WORKSPACEPATH -maxdepth 1 -name "*.fast5" -print0 | xargs  -0 cp -t $ANNOTPATH;

rm -r $SAVEPATH
echo "Finished Guppy basecalling for "$SAMPLEPATH

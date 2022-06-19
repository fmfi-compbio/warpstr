#!/bin/bash
#wrapper script traversing rootpath and calling guppy_annotate.sh for each sample_run

function usage {
        echo "Usage: $(basename $0) [-rgfkat]" 2>&1
        echo '   -r   rootpath for reads of various samples'
        echo '   -g   path to guppy basecaller'
        echo '   -f   flowcell name'
        echo '   -k   kit name'
        echo '   -a   subdir name for output'
        echo '   -t   number of threads'
        exit -1
}

if [[ ${#} -eq 0 || ${#} -le 9 ]]; then
   usage
fi

while getopts r:g:f:k:a:t: arg; do
  case "${arg}" in
    r) ROOTPATH=${OPTARG};;
    g) GUPPY=${OPTARG};;
    f) FLOWCELL=${OPTARG};;
    k) KIT=${OPTARG};;
    a) OUT_SUFFIX=${OPTARG};;
    t) THREADS=${OPTARG};;
    
    ?)
      echo "Invalid option: -${OPTARG}."
      echo
      usage
      ;;
  esac
done

#get path to the second script calling guppy for each fast5 file
ANNOT=$(dirname "$(readlink -f "$BASH_SOURCE")")"/guppy_annotate.sh"

#check if guppy params (flowcell and kit) are correct
test_guppy=`$GUPPY --print_workflows | grep "^$FLOWCELL " | grep $KIT" " | wc -l`
if [[ $test_guppy -eq 0 ]]; then
   echo "incorrect guppy parameters: flowcell "$FLOWCELL "and kit "$KIT
   exit 1
fi

#call guppy annotate for each sample
find $ROOTPATH -mindepth 1 -maxdepth 2 -type d -print0 | xargs -0 -I{} $ANNOT {} $GUPPY $FLOWCELL $KIT $OUT_SUFFIX $THREADS

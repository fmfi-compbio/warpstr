#~/bin/bash

echo Running WarpSTR test run case
echo Please provide path to the GRCh38 reference, i.e. /path/GRCh38.fa
read REF
echo Please provide path to the guppy executable, i.e. /path/guppy_basecaller
read GUPPY

cp 'test/config_template.yaml' 'test/config_test.yaml'
sed -i "s|PLACEHOLDER_REFERENCE|$REF|" 'test/config_test.yaml'
sed -i "s|PLACEHOLDER_GUPPY|$GUPPY|" 'test/config_test.yaml'

echo "Running WarpSTR with config: 'test/config_test.yaml'"
echo "... This could take approx. 3-5 minutes ..."
python WarpSTR.py test/config_test.yaml 2> test/log.err

hmm=1

if [ ! -e "bin" ]; then mkdir bin; fi

echo "Starting SEGMENTS.pl"
pp -f Bleach -o bin/dmrs src/SEGMENTS.pl
echo "Starting pipeline.pl"
pp -f Bleach -o bin/MethylScore src/pipeline.pl
echo "Starting generate_genome_matrix.pl"
pp -f Bleach -o bin/generate_genome_matrix src/generate_genome_matrix.pl
echo "Starting split_MRfile.pl"
pp -f Bleach -o bin/split_MRfile src/split_MRfile.pl
echo "Starting merge_DMRs.pl"
pp -f Bleach -o bin/merge_DMRs src/merge_DMRs.pl

echo "Starting context-dep bins"
pp -f Bleach -o bin/dmrs-contexts src/SEGMENTS-contexts.pl
pp -f Bleach -o bin/MethylScore-contexts src/pipeline-contexts.pl
pp -f Bleach -o bin/merge_DMRs-contexts src/merge_DMRs-contexts.pl

dir=`pwd`

# compile hmm_mrs
if [ $hmm -eq 1 ]; then
    cd src/hmm_mrs/src/smithlab_cpp
    make clean && make
    cd ../analysis
    export SMITHLAB_CPP=../smithlab_cpp
    make
    cd $dir
    cp src/hmm_mrs/src/analysis/hmm_mrs bin/
    chmod +x bin/hmm_mrs
else
    cp src/hmm_mrs/src/analysis/hmm_mrs bin/
fi
###

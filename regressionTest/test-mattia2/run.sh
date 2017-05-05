#!/bin/bash

OLDDIR=`pwd`

JSON=full.json
GTF=annotated-isoforms.gtf
GENOME=genomic.txt
ESTS=ests.txt


BINDIR=/home/dellavedova/Devel/Splicing/PIntron/dist/pintron-latest/bin
cd $BINDIR/../../..
export USE_PAR=no
make dist
cd $OLDDIR

rm -f $GTF *.log
rm -f $JSON \
TEMP_COMPOSITION_TRANS1_* TRANSCRIPTS1_* VariantGTF.txt  pintr* raw-multifasta-out.txt \
build-ests.txt meg-edges.txt  out-after-intron-agree.txt  predicted-introns.txt  processed-megs-info.txt  \
CCDS_transcripts.txt  config-dump.ini \
genomic-exonforCCDS.txt  megs.txt  out-agree.txt  processed-ests.txt \
processed-megs.txt ests-alignments.sam


$BINDIR/pintron -c -k -a --bin-dir=$BINDIR -g $GENOME -s $ESTS \
--output=$JSON  \
--gtf=$GTF \
--extended-gtf=pintron-all-isoforms.gtf \
--logfile=pintron-pipeline-log.txt \
--general-logfile=pintron-log.txt \
--organism="human"  --gene=AAMP  \
2>&1 > pintron.log

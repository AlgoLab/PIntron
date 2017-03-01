The PIntron pipeline consists of the follwing steps

|**Program**            |Output                         |Explanation                               |  
|-----------------------|-------------------------------|------------------------------------------|  
|est-fact               |raw-multifasta-out.txt, processed-ests.txt         |Alignments of EST sequences to the genome |  
|min-factorization      |out-agree.txt                  |Minimize the EST alignments               |   
|intron-agreement       |out-after-intron-agree.txt, predicted-introns.txt              |Compute the agreement among the introns   |  
|compact-compositions   |build-ests.txt, genomic-exonforCCDS.txt                 |Compact the EST alignments                |  
|maximal-transcripts    |isoforms.txt                   |Compute all full-length transcripts       |  
|cds-annotation         |CCDS_transcripts.txt, variantGTF.txt           |Annotate the full-length transcripts      |  

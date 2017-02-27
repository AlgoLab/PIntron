The PIntron pipeline consists of the follwing steps

|**Program**            |Output                         |Explanation                               |  
|-----------------------|-------------------------------|------------------------------------------|  
|est-fact               |raw-multifasta-out.txt         |Alignments of EST sequences to the genome |  
|min-factorization      |out-agree.txt                  |Minimize the EST alignments               |   
|intron-agreement       |intron-agreement               |Compute the agreement among the introns   |  
|compact-compositions   |build-ests.txt                 |Compact the EST alignments                |  
|maximal-transcripts    |isoforms.txt                   |Compute all full-length transcripts       |  
|cds-annotation         |CCDS_transcripts.txt           |Annotate the full-length transcripts      |  

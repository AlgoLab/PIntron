The PIntron pipeline consists of the follwing steps

|**Program**            |Output                         |Explanation                            |  
|-----------------------|-------------------------------|---------------------------------------|  
|est-fact               |raw-multifasta-out.txt         |Factorizations                         |  
|min-factorization      |out-agree.txt                  |Minimize the factorization             |   
|intron-agreement       |intron-agreement               |Compute the agreement                  |  
|compact-compositions   |build-ests.txt                 |Compact the compositions               |  
|maximal-transcripts    |isoforms.txt                   |Compute all full-length transcripts    |  
|cds-annotation         |CCDS_transcripts.txt           |Select the annotated transcripts       |  

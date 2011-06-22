### JSON Output

The output file specified with the `--output` option is a JSON file, containing informations
  on the gene structure, the isoforms and the exons. Since JSON is mainly a
  (key, value) dictionary, we have chosen keys that are as self-explaining as
  possible. We have exploited the nesting nature of JSON files to
encode a set of gene where each gene has a set of isoforms and a set of
  introns.

Each gene has keys  "length_genomic_sequence",
  "version",  "number_isoforms", "isoforms" (a hash where a progressive ID is
  the only key and each value is a single isoform), "introns" (a list).

Each isoform has keys "canonical start codon?",
  "canonical end codon?", "NMD flag", "polyA?", "annotated CDS?", "start", "PAS?"
  "end", "number exons", "Type", "length", "coding length", "reference?", "from RefSeq?",
  "exons" (a list of exons).

Each exon has keys  "sequence", "relative start", "relative end",
  "chromosome start", "chromosome end" (absolute coordinates), "3utr length",
  "5utr length".

Each intron has keys "relative start", "relative end", "chromosome start",
"chromosome end", "length", "number supporting EST", "EST_list" (a list of EST
IDs), "donor alignment average error", "acceptor alignment average error", "donor score",
"acceptor score", "BPS score", "BPS position", "type", "pattern",
"repeat sequence", "donor suffix", "prefix", "suffix",
"acceptor prefix"

#!/usr/bin/env python3
####
#
#
#                              PIntron
#
# A novel pipeline for computational gene-structure prediction based on
# spliced alignment of expressed sequences (ESTs and mRNAs).
#
# Copyright (C) 2010  Gianluca Della Vedova, Yuri Pirola
#
# Distributed under the terms of the GNU Affero General Public License (AGPL)
#
#
# This file is part of PIntron.
#
# PIntron is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIntron is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with PIntron.  If not, see <http://www.gnu.org/licenses/>.
#
####


import random
import sys
import re
import os
import subprocess
import time
import logging
import json
import pprint
import traceback
import csv

# from Bio import SeqIO
from optparse import OptionParser


class PIntronError(Exception):
    """Base class for exceptions of the PIntron pipeline."""
    pass


class PIntronIOError(PIntronError):
    """Exception raised for errors related to I/O operations.

    Attributes:
        e_file -- file on which the error occurred
        msg    -- explanation of the error
    """

    def __init__(self, e_file, msg):
        self.e_file = e_file
        self.msg = msg

    def __str__(self):
        return '{0} (offending file: "{1}")'.format(self.msg, self.e_file)


def parse_command_line():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    # parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
    #                   default=False, help="print status messages to stdout")
    parser.add_option("-g", "--genomic",
                      dest="genome_filename",
                      default="genomic.txt",
                      help="FILE containing the genomic sequence",
                      metavar="GENOMIC_FILE")
    parser.add_option("-s", "--est",
                      dest="est_filename", default="ests.txt",
                      help="FILE containing the ESTs", metavar="ESTs_FILE")
    # parser.add_option("-1", "--no-full-lengths", action="store_true",
    #                   dest="step1", default=False,
    #                   help="do not compute the full-lengths")
    parser.add_option("-o", "--output",
                      dest="output_filename",
                      default="pintron-full-output.json",
                      help="full output file (default = 'pintron-full-output.json')",
                      metavar="FILE")
    parser.add_option("-z", "--compress", action="store_true",
                      dest="compress", default=False,
                      help="compress output (default = False)")
    parser.add_option("-l", "--logfile",
                      dest="logfile", default="pintron-log.txt",
                      help="log filename (default = 'pintron-log.txt')")
    parser.add_option("-c", "--continue", action="store_true",
                      dest="from_scratch", default=False,
                      help="resume a previosly interrupted computation (default = False)")
    parser.add_option("-b", "--bin-dir",
                      dest="bindir", default="",
                      help="DIRECTORY containing the programs (default = system PATH)")
    parser.add_option("-n", "--organism",
                      dest="organism", default="unknown",
                      help="Organism originating the ESTs (default = 'unknown')")
    parser.add_option("-e", "--gene",
                      dest="gene", default="unknown",
                      help="Gene symbol (or ID) of the locus which the ESTs refer to (default = 'unknown')")
    parser.add_option("-k", "--keep-intermediate-files", action="store_true",
                      dest="no_clean", default=False,
                      help="keep all intermediate or temporary files (default = False)")
    parser.add_option("-t", "--gtf",
                      dest="gtf_filename",
                      default="pintron-cds-annotated-isoforms.gtf",
                      help="output GTF FILE with the isoforms that have a CDS annotation "
                      "(default = 'pintron-cds-annotated-isoforms.gtf')",
                      metavar="GTF_FILE")
    parser.add_option("--extended-gtf",
                      dest="extended_gtf_filename",
                      default="pintron-all-isoforms.gtf",
                      help="output GTF FILE with all the predicted isoforms "
                      "(default = 'pintron-all-isoforms.gtf')",
                      metavar="GTF_FILE")

    # parser.add_option("--strand",
    #                   dest="strand", type="int", default=1,
    #                   help="[Expert use only] print status messages to stdout")
    # parser.add_option("--chromosome",
    #                   dest="chromosome", default=False,
    #                   help="[Expert use only] print status messages to stdout")
    # parser.add_option("--est-cluster",
    #                   dest="est_cluster", default=False,
    #                   help="[Expert use only] print status messages to stdout")
    # parser.add_option("--min-factor-length",
    #                   dest="min_factor_length", type="int", default=15,
    #                   help="[Expert use only] minimum factor length")
    # parser.add_option("--min-intron-length",
    #                   dest="min_intron_length", type="int", default=60,
    #                   help="[Expert use only] minimum intron length")
    # parser.add_option("--max-intron-length",
    #                   dest="max_intron_length", type="int", default=0,
    #                   help="[Expert use only] max intron length")
    # parser.add_option("--min-string-depth-rate",
    #                   dest="min_string_depth_rate", type="float", default=0.2,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-prefix-discarded-rate",
    #                   dest="max_prefix_discarded_rate", type="float", default=0.6,
    #                   help="[Expert use only] largest prefix of an EST that might be discarded (expressed as a fraction of the EST length)")
    # parser.add_option("--max-suffix-discarded-rate",
    #                   dest="max_suffix_discarded_rate", type="float", default=0.6,
    #                   help="[Expert use only] largest suffix of an EST that might be discarded (expressed as a fraction of the EST length)")
    # parser.add_option("--max-prefix-discarded",
    #                   dest="max_prefix_discarded", type="int", default=50,
    #                   help="[Expert use only] largest prefix of an EST that might be discarded during factorization construction (expressed as number of bases)")
    # parser.add_option("--max-suffix-discarded",
    #                   dest="max_suffix_discarded", type="int", default=50,
    #                   help="[Expert use only] largest suffix of an EST that might be discarded during factorization construction (expressed as number of bases)")
    # parser.add_option("--min-distance-of-splice-sites",
    #                   dest="min_distance_splice_sites", type="int", default=50,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-no-of-factorizations",
    #                   dest="max_factorizations", type="int", default=0,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-difference-of-coverage",
    #                   dest="max_difference_coverage", type="float", default=0.05,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-difference-of-no-of-exons",
    #                   dest="max_difference_no_exons", type="int", default=5,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-difference-of-gap-length",
    #                   dest="max_difference_gap_length", type="int", default=20,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--retain-externals",
    #                   dest="retain_externals", default=True, action="store_true",
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-pairings-in-CMEG",
    #                   dest="max_pairings_CMEG", type="int", default=80,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--max-shortest-pairing-frequence",
    #                   dest="max_shortest_pairing_frequence", type="float", default=0.4,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--suffix-prefix-length-intron",
    #                   dest="suffix_prefix_length_intron", type="int", default=70,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--suffix-prefix-length-est",
    #                   dest="suffix_prefix_length_est", type="int", default=30,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--suffix-prefix-length-genomic",
    #                   dest="suffix_prefix_length_genomic", type="int", default=30,
    #                   help="[Expert use only] TODO")
    # parser.add_option("--no-transitive-reduction",
    #                   dest="no_transitive_reduction", default=False, action="store_true",
    #                   help="[Expert use only] TODO")
    # parser.add_option("--no-short-edge",
    #                   dest="no_short_edge", default=False, action="store_true",
    #                   help="[Expert use only] TODO")
    parser.add_option("--set-max-factorization-time",
                      dest="max_factorization_time", type="int", default=60,
                      help="[Expert use only] Set a time limit (in mins) for the factorization step")
    parser.add_option("--set-max-factorization-memory",
                      dest="max_factorization_memory", type="int", default=3000,
                      help="[Expert use only] Set a limit (in MiB) for the memory used by the factorization step"
                      " (default = 3000 MiB, approx. 3GB)")
    parser.add_option("--set-max-exon-agreement-time",
                      dest="max_exon_agreement_time", type="int", default=15,
                      help="[Expert use only] Set a time limit (in mins) for the exon agreement step")
    parser.add_option("--set-max-intron-agreement-time",
                      dest="max_intron_agreement_time", type="int", default=30,
                      help="[Expert use only] Set a time limit (in mins) for the intron agreement step")
    parser.add_option("--pas-tolerance",
                      dest="pas_tolerance", type="int", default=30,
                      help="[Expert use only] Maximum allowed difference on the exon final coordinate to identify a PAS")

    (options, args) = parser.parse_args()
    if options.bindir:
        options.bindir = os.path.normpath(options.bindir)

    return(options)


# Transform a JSON file into a GTF
def json2gtf(infile, outfile, genomic_seq, gene_name, all_isoforms):
    def write_gtf_line(file, seqname, feature, start, end, score, strand, frame, gene, transcript):
        if end >= start:
            f.write("\t".join([seqname, "PIntron", feature, str(start), str(end), score, strand, str(frame),
                           "gene_id \"{0}\"; transcript_id \"{0}.{1}\";\n".format(gene, str(transcript))]))

    logging.debug(str(time.localtime()))
    logging.debug(json2gtf)
    logging.debug(infile + ':' + outfile + ':' + genomic_seq)
    # Find the sequence ID
    # It is stored in the first line of the genomic sequence
    with open(genomic_seq, 'r', encoding='utf-8') as f:
        line=f.readline().rstrip("\r\n")
        m=re.search('^>[^\d]*(\d+):.*:([+-]?\d+)$',line)
        sequence_id = line[1:]
        strand=m.group(2)
        if strand == '-1' or strand == '-':
            strand = '-'
        else:
            strand = '+'

    with open(infile, 'r', encoding='utf-8') as f:
        entry = json.load(f)

    # The genomic sequence length stored in the JSON file
    # cannot be trusted.
    # seq_record=next(SeqIO.parse(genomic_seq, "fasta"))
    # entry['length_genomic_sequence']=len(seq_record)
    # print(entry['length_genomic_sequence'])

    with open(outfile, 'w', encoding='utf-8') as f:
        data_strand = {'first': {"label": "5UTR",
                                 "delimiter": "start_codon",
                                 "codons": ["ATG"],
                                 "offset": 0,
                             },
                       'last' : {"label" : "3UTR",
                                 "delimiter" : "stop_codon",
                                 "codons"        : ["TGA", "TAG", "TAA"],
                                 "offset"        : 3,
                             }
                   }
        # if strand == '-':
        #     data_strand['temp']=data_strand['first']
        #     data_strand['first']=data_strand['last']
        #     data_strand['last']=data_strand['temp']
        #     del data_strand['temp']

        #pprint.pprint(data_strand)
        for isoform_id, isoform in entry["isoforms"].items():
            if not all_isoforms or not isoform["annotated CDS?"]:
                continue
            whole_cds_len = 0
            total_cds_length = isoform['CDS length'] -3 # Because the stop codon is outside the CDS
            for p in ['first', 'last']:
                data_strand[p]['codon']    = ''
                data_strand[p]['codon_ok'] = False
            if strand == '-':
                isoform["exons"].reverse()

            for exon in isoform["exons"]:
                # if the transcript is taken with a direction opposite to
                # that of the genome, then chromosome end < chromosome start
                rel_start   = exon['relative start']
                rel_end     = exon['relative end']
                exon_length = len(exon['sequence'])
                bad_prefix  = exon.get("5utr length", 0)
                bad_suffix  = exon.get("3utr length", 0)
                cds_start   = rel_start + bad_prefix
                cds_end     = rel_end - bad_suffix
                frame       = (3-(whole_cds_len % 3)) % 3 #http://mblab.wustl.edu/GTF22.html

#                pprint.pprint(data_strand)
                if cds_end > cds_start:
                    whole_cds_len += cds_end - cds_start +1

                for p in ['first', 'last']:
                    data_strand[p]['print'] = False
                    data_strand[p]['new'] = ''
                    if not data_strand[p]['codon_ok']:
                        if p == 'first':
                            data_strand['first']['new'] = exon['sequence'][bad_prefix:bad_prefix+3]
                        else:
                            if whole_cds_len > total_cds_length:
                                data_strand['last']['new'] = exon['sequence'][-(bad_suffix+3):][:3]
                cds_start += min(len(data_strand['first']['new']),data_strand['first']['offset'])
                cds_end   -= min(len(data_strand['last']['new']),data_strand['last']['offset'])

                for p in ['first', 'last']:
                    if len(data_strand[p]['new']) > 0 and not data_strand[p]['codon_ok']:
                        data_strand[p]['codon'] += data_strand[p]['new']
                        data_strand[p]['print'] = True
                        if not data_strand[p]['codon_ok'] and isoform["annotated CDS?"] and not (data_strand[p]['codon'].upper() in data_strand[p]['codons']):
                            print("Warning: wrong " + data_strand[p]['delimiter'] +
                                  ". Found " + data_strand[p]['codon'] + " instead of " +
                                  "/".join(data_strand[p]['codons']))
                            pprint.pprint(exon)
                            print (exon_length)
                        if len(data_strand[p]['codon']) == 3:
                            data_strand[p]['codon_ok'] = True

                write_gtf_line(f, sequence_id, "exon", rel_start, rel_end, "0", "+", ".", gene_name, isoform_id)
                write_gtf_line(f, sequence_id, data_strand['first']['label'], rel_start, cds_start-1, "0", "+",  ".", gene_name, isoform_id)
                if data_strand['first']['print']:
                    write_gtf_line(f, sequence_id, data_strand['first']['delimiter'], cds_start, cds_start+2,
                                   "0", "+", len(data_strand['first']['codon']) % 3, gene_name, isoform_id)
                write_gtf_line(f, sequence_id, "CDS", cds_start, cds_end, "0", "+", frame, gene_name, isoform_id)
                if data_strand['last']['print']:
                    write_gtf_line(f, sequence_id, data_strand['last']['delimiter'], cds_end-2+3, cds_end+3,
                                   "0", "+", len(data_strand['last']['codon']) % 3, gene_name, isoform_id)
                write_gtf_line(f, sequence_id, data_strand['last']['label'], cds_end+1+len(data_strand['last']['new']), rel_end, "0", "+",  ".", gene_name, isoform_id)


def compute_json(ccds_file, variant_file, logfile, output_file, from_scratch, pas_tolerance):
    gene={
        'version': 2, # Hardcoding version number
        'isoforms': {},
        'introns': {},
        'exons': {},
    }

    for file in [ccds_file, variant_file]:
        if not os.access(file, os.R_OK):
            # throw exception and die
            logging.exception("*** Fatal error: Could not read " + file + "\n")

    factorizations = {}

    with open('out-after-intron-agree.txt', mode='r', encoding='utf-8') as fd:
        current = ''
        for line in fd:
            l = line.rstrip()
            if l[0] == '>':
                current = l
                factorizations[current] = {
                    'polyA?' : False,
                    'PAS' : False,
                }
            elif re.match('#polya=1', l):
                factorizations[current]['polyA?'] = True
            elif re.match('#polyad(\S*)=1', l):
                factorizations[current]['PAS'] = True
            elif factorizations[current]['PAS'] and re.match('(\d+) (\d+) (\d+) (\d+)( \S+)? \S+$', l):
                # Since we use the factorizations only for detecting PAS, there is no need
                # for storing unused information
                new = re.match('(\d+) (\d+) (\d+) (\d+)( \S+)? \S+$', l).groups()
                exon = {
                    'relative start'   : int(new[2]),
                    'relative end'     : int(new[3]),
                }
                factorizations[current]['exon']=exon
    # At the end, remove all factorizations without exons, since they are not useful
    PAS_factorizations = {k:v for k,v in factorizations.items() if factorizations[k]['PAS'] }
#    pprint.pprint(PAS_factorizations)
#    del factorizations

    with open(variant_file, mode='r', encoding='utf-8') as fd:
        for line in fd:
            row = re.split(' /', line.rstrip())
            index=int(re.sub('^.*\#', '', row.pop(0)))
            isoform = {
                'exons' : [],
                'polyA?' : False,
                'PAS?' : False,
                'annotated CDS?' : False,
                'Reference frame?' : False,
            }
            for t in row:
                (k, v) = re.split('=', t, 2)
                logging.debug("Reading VariantGTF: " + k + "=>" + v + "!\n")
                print("Reading VariantGTF: " + k + "=>" + v + "!\n")

                if k == "nex":
                    isoform['number exons'] = int(v)
                elif k == "L":
                    isoform["length"] = int(v)
                elif k == "CDS":
                    if v != '..':
                        isoform["annotated CDS?"] = True
                        isoform['CDS length'] = 0
                        m = re.match('^(<?)(\d+)\.\.(\d+)(>?)$', v)
                        (a, isoform["CDS start"], isoform["CDS end"], b) = (m.group(1), int(m.group(2)),
                                                                            int(m.group(3)), m.group(4))
                        isoform['canonical start codon?'] = False if a == '<' else True
                        isoform['canonical end codon?']   = False if b == '>' else True
                elif k == "RefSeq":
                    m = re.match('^(.*?)(\(?([NY])([NY])\)?)?$', v, flags=re.IGNORECASE)
                    if m:
                        (r, a, b) = (m.group(1), m.group(3), m.group(4))
                        isoform['reference CDS start codon?'] = False if a == 'N' else True
                        isoform['reference CDS end codon?']   = False if b == 'N' else True
                        if r != None and r:
                            isoform['RefSeq'] = r
                elif k == "ProtL":
                    if v != '..' and isoform["annotated CDS?"]:
                        m = re.match('^(>?)(\d+)$', v, flags=re.IGNORECASE)
                        (a, isoform['protein length']) = (m.group(1), int(m.group(2)))
                        isoform['protein missing codon?'] = False if a != '>' else True
                elif k == "Frame":
                    m = re.match('^y', v, flags=re.IGNORECASE)
                    if m != None and isoform["annotated CDS?"]:
                        isoform['Reference frame?'] = True
                elif k == "Type":
                    ref = True if v == 'Ref' else False
                    # if ref != isoform['Reference frame?']:
                    #     raise ValueError(format("Wrong reference for isoform n. {}\n{}",
                    #                             str(index), line))
                    if isoform['Reference frame?']:
                        if 'RefSeq' in isoform:
                            isoform['Type'] = isoform['RefSeq'] + " (Reference TR)"
                        else:
                            isoform['Type'] = "(Reference TR)"
                    else:
                        isoform['Type'] = re.sub('\s+$', '', v)
                elif not re.match('^\s*\#', line):
                    raise ValueError("Could not parse GTF file " + variant_file + "(" + k + "=>" + v + ")\n" + line + "\n")
            gene['isoforms'][index] = isoform

    with open(ccds_file, mode='r', encoding='utf-8') as fd:
        gene['number_isoforms'] = int(fd.readline().rstrip())
        gene['length_genomic_sequence'] = int(fd.readline().rstrip())
        for line in fd:
            l = line.rstrip()
            l = re.sub('\s+', '', l)
            l = re.sub('#.*', '', l)

            if re.match('^>', l):
                # New isoform
                l = l[1:]
                # print(l)
                fields = [int(x) for x in  re.split(':', l)]

                # pprint.pprint(fields)
                index = fields[0]
                if not index in gene['isoforms']:
                    raise ValueError("CCDS file " + ccds_file + "contains isoform with index " + index + " not in variant_file\n")
                if fields[1] > gene['isoforms'][index]['number exons']:
                    import pdb; pdb.set_trace()
                    raise ValueError("Wrong number of exons: " + str(index) + "\n " + str(fields[1]) + "!= " +
                                     str(isoform['number exons']) +"\n")

                gene['isoforms'][index]['reference?'] = False if fields[2] == 0 else True
                gene['isoforms'][index]['from RefSeq?'] = False if fields[3] == 0 else True
                gene['isoforms'][index]['NMD flag'] = fields[4]

            elif re.match('^(\d+:){5}(-?\d+:)(-?\d+)$', l):
                # Row contains exon metadata
                exon = {}
                (exon["chromosome start"], exon["chromosome end"], exon["relative start"], exon["relative end"],
                 polyA, exon["5utr length"], exon["3utr length"]) = [max(0, int(x)) for x in  re.split(':', l)]
                if (polyA == 1):
                    gene['isoforms'][index]['polyA?'] = True
                # pprint.pprint(exon)
                logging.debug("Reading CCDS_transcripts: Row contains exon metadata\n")
                logging.debug(line)
                logging.debug(max(exon["relative end"], exon["relative start"]))
                logging.debug(min(exon["relative end"], exon["relative start"]))
                logging.debug(exon["5utr length"])
                logging.debug(exon["3utr length"])
                if gene['isoforms'][index]['annotated CDS?']:
                    gene['isoforms'][index]["CDS length"] += (max(exon["relative end"], exon["relative start"]) -
                                                             min(exon["relative end"], exon["relative start"]) + 1 -
                                                             exon["5utr length"] - exon["3utr length"])

                if int(re.split(':', l)[4]) < 0:
                    del(exon["5utr length"])
                if int(re.split(':', l)[5]) < 0:
                    del(exon["3utr length"])
                gene['isoforms'][index]['exons'].append(exon)

            elif re.match('^[acgtACGT]+$', l):
                last_exon = gene['isoforms'][index]['exons'][-1]['sequence'] = l
            elif not re.match('^\s*\#', line):
                raise ValueError("Could not parse CCDS file " + ccds_file + " at line:\n" + line + "\n")

    with open('predicted-introns.txt', mode='r', encoding='utf-8') as fd:
        index=1
        for line in fd:
            intron = {}
            (intron['relative start'], intron['relative end'],
             intron['chromosome start'], intron['chromosome end'], intron['length'], intron['number supporting EST'], EST_list,
             intron['donor alignment average error'], intron['acceptor alignment average error'], intron['donor score'],
             intron['acceptor score'], intron['BPS score'], intron['BPS position'], intron['type'], intron['pattern'],
             intron['repeat sequence'], intron['donor suffix'], intron['prefix'], intron['suffix'],
             intron['acceptor prefix']) = re.split("\t", line.rstrip())
            intron['EST list'] = [i for i in re.split(',', EST_list) if i != '']

            for field in ('relative start', 'relative end', 'chromosome start', 'chromosome end', 'length', 'number supporting EST',
                          'BPS position'):
                intron[field] = int(intron[field])
            for field in ('donor alignment average error', 'acceptor alignment average error', 'donor score',
                          'acceptor score', 'BPS score'):
                intron[field] = float(intron[field])

            if intron['BPS position'] < 0:
                del intron['BPS position']

            gene['introns'][index]=intron
            index += 1

    # add introns to each isoform
    for isoform in gene['isoforms'].values():
        isoform['exons'].sort(key=lambda x: x['relative end'])
        isoform['introns'] = []
        pairs = zip(isoform['exons'][1:], isoform['exons'][:-1])
#        import pdb; pdb.set_trace()
        for pair in pairs:
            list_extremes = [pair[0]['chromosome end'], pair[0]['chromosome start'],
                             pair[1]['chromosome end'], pair[1]['chromosome start']]
            list_extremes.sort()
            left_border = list_extremes[1]+1
            right_border = list_extremes[2]-1
            for index in gene['introns'].keys():
                intron=gene['introns'][index]
                if intron['chromosome start'] == left_border and intron['chromosome end'] == right_border or intron['chromosome end'] == left_border and intron['chromosome start'] == right_border:

                    isoform['introns'].append(index)

    def same_coordinates(a, b):
        return True if (a['relative start'] == b['relative start'] and
                        30 >=  a['relative end'] - b['relative end'] >= -30) else False


    # pprint.pprint(gene)
    for isoform in gene['isoforms'].keys():
        # Check if we have to add PAS
        if not gene['isoforms'][isoform]['polyA?']:
            continue
        exon = gene['isoforms'][isoform]['exons'][-1]
        # If PAS_factorizations has an exon with the same coordinates,
        # we have a PAS
        if any(x for x in PAS_factorizations.values() if same_coordinates(x['exon'], exon)):
            gene['isoforms'][isoform]['PAS?']=True

    # import pdb; pdb.set_trace()
    with open(output_file, mode='w', encoding='utf-8') as fd:
        fd.write(json.dumps(gene, sort_keys=True, indent=4))


def exec_system_command(command, error_comment, logfile, output_file="",
                        from_scratch=True):
    if not from_scratch or (not output_file == "" and not os.access(output_file, os.R_OK)):
        logging.debug(str(time.localtime()))
        logging.debug(command)

        try:
            retcode = subprocess.call(command + " 2>> " + logfile, shell=True)
            if retcode < 0:
                print(error_comment, -retcode, file=sys.stderr)
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)
            raise PIntronError


def check_executables(bindir, exes):
    """Check if the executables are in the path or in the specified directory.
    """

    full_exes = {}
    if bindir:
        if bindir[0] == '~':
            bindir = os.environ["HOME"] + bindir[1:]
        paths = [bindir] + os.environ["PATH"].split(os.pathsep)
    else:
        paths = os.environ["PATH"].split(os.pathsep)
    for exe in exes:
        full_exes[exe] = None
        for path in paths:
            if os.access(os.path.join(path, exe), os.X_OK):
                real_path = os.path.realpath(os.path.abspath(os.path.join(path, exe)))
                logging.debug("Using program '{}' in dir '{}'".format(exe, real_path))
                full_exes[exe] = real_path
                break
        if full_exes[exe] == None:
            raise PIntronIOError(exe, "Could not find program '{}'!\n"
                                 "Search path: {}".format(exe,
                                                          ":".join(paths)))
    return full_exes


def pintron_pipeline(options):
    """Executes the whole pipeline, using the input options.
    """

    pintron_version='_____%PINTRON_VERSION%_____'
    if (pintron_version[0] == '_' and
        pintron_version[1:] == '____%PINTRON_VERSION%_____'):
        pintron_version= ''

    logging.info("PIntron%s", pintron_version)
    logging.info("Copyright (C) 2010  Paola Bonizzoni, Gianluca Della Vedova, Yuri Pirola, Raffaella Rizzi.")
    logging.info("This program is distributed under the terms of the GNU Affero General Public License (AGPL), either version 3 of the License, or (at your option) any later version.")
    logging.info("This program comes with ABSOLUTELY NO WARRANTY. See the GNU Affero General Public License for more details.")
    logging.info("This is free software, and you are welcome to redistribute it under the conditions specified by the license.")

    logging.info("Running: " + " ".join(sys.argv))


    # Check and copy input data
    logging.info("STEP  1:  Checking executables and preparing input data...")

    exes= check_executables(options.bindir, ["est-fact",
                                             "min-factorization",
                                             "intron-agreement",
                                             "gene-structure",
                                             "compact-compositions",
                                             "maximal-transcripts",
                                             "cds-annotation"
                                             ])

    if not os.path.isfile(options.genome_filename) or not os.access(options.genome_filename, os.R_OK):
        raise PIntronIOError(options.genome_filename,
                             'Could not read file "' + options.genome_filename + '"!')
    if not os.path.isfile(options.est_filename) or not os.access(options.est_filename, os.R_OK):
        raise PIntronIOError(options.est_filename,
                             'Could not read file "' + options.est_filename + '"!')
    if ( os.access('genomic.txt', os.F_OK) and
         not os.path.samefile('genomic.txt', options.genome_filename) and
         not os.access('genomic.txt', os.W_OK) ):
        raise PIntronIOError('genomic.txt',
                             'Could not write file "genomic.txt"!')
    if ( os.access('ests.txt', os.F_OK) and
         not os.path.samefile('ests.txt', options.est_filename) and
         not os.access('ests.txt', os.W_OK) ):
        raise PIntronIOError('ests.txt',
                             'Could not write file "ests.txt"!')

    if os.path.isfile('genomic.txt') and os.path.samefile('genomic.txt', options.genome_filename):
        logging.debug('Files "%s" and "genomic.txt" refer to the same file: skip copy.',
                      options.genome_filename)
    else:
        exec_system_command(
            command="cp " + options.genome_filename + " genomic.txt ",
            error_comment="Could not prepare genomic input file",
            logfile=options.logfile,
            output_file='raw-multifasta-out.txt',
            from_scratch=options.from_scratch)

    if os.path.isfile('ests.txt') and os.path.samefile('ests.txt', options.est_filename):
        logging.debug('Files "%s" and "ests.txt" refer to the same file: skip copy.',
                      options.est_filename)
    else:
        exec_system_command(
            command="cp " + options.est_filename + " ests.txt ",
            error_comment="Could not prepare ESTs input file",
            logfile=options.logfile,
            output_file='raw-multifasta-out.txt',
            from_scratch=options.from_scratch)


    # Compute factorizations
    logging.info("STEP  2:  Pre-aligning transcript data...")

    exec_system_command(
        command="ulimit -t "+ str(options.max_factorization_time * 60) + " && ulimit -v " +
        str(options.max_factorization_memory * 1024) + " && " + exes["est-fact"],
        error_comment="Could not compute the factorizations",
        logfile=options.logfile,
        output_file='raw-multifasta-out.txt',
        from_scratch=options.from_scratch)


    # Min factorization agreement
    logging.info("STEP  3:  Computing a raw consensus gene structure...")

    exec_system_command(
        command="ulimit -t "+ str(options.max_exon_agreement_time * 60) +" && " +
        exes["min-factorization"] + " < raw-multifasta-out.txt >out-agree.txt",
        error_comment="Could not minimize the factorizations",
        logfile=options.logfile,
        output_file='out-agree.txt',
        from_scratch=options.from_scratch)


    # Intron prediction
    logging.info("STEP  4:  Predicting introns...")

    exec_system_command(
        command="ulimit -t "+ str(options.max_intron_agreement_time*60) +" && " +
        exes["intron-agreement"],
        error_comment="Could not compute the factorizations",
        logfile=options.logfile,
        output_file='out-after-intron-agree.txt',
        from_scratch=options.from_scratch)


    # Compute the Gene structure
    logging.info("STEP  5:  Computing the final gene structure...")

    exec_system_command(
        command=exes["gene-structure"] + " ./",
        error_comment="Could not compute maximal transcripts",
        logfile=options.logfile,
        output_file='out-after-intron-agree.txt',
        from_scratch=options.from_scratch)

    # TODO: Gene structure browser

    # The computation of the full-length isoforms should not be avoided
    # if options.step1:
    #     sys.exit(0)


    # Transform compositions into exons
    logging.info("STEP  6:  Computing the final transcript alignments...")

    exec_system_command(
        command=exes["compact-compositions"] + " < out-after-intron-agree.txt > build-ests.txt",
        error_comment="Could not transform factorizations into exons",
        logfile=options.logfile,
        output_file='build-ests.txt',
        from_scratch=options.from_scratch)


    # Compute maximal transcripts
    logging.info("STEP  7:  Computing the final full-length isoforms...")

    exec_system_command(
        command=exes["maximal-transcripts"] + " < build-ests.txt",
        error_comment="Could not compute maximal transcripts",
        logfile=options.logfile,
        output_file='CCDS_transcripts.txt',
        from_scratch=options.from_scratch)
    exec_system_command(
        command="cp -f TRANSCRIPTS1_1.txt isoforms.txt",
        error_comment="Could not link isoforms",
        logfile=options.logfile,
        output_file='CCDS_transcripts.txt',
        from_scratch=options.from_scratch)


    # Annotate CDS
    logging.info("STEP  8:  Annotating CDS...")

    exec_system_command(
        command=exes["cds-annotation"] + " ./ ./ " + options.gene + " " + options.organism,
        error_comment="Could not annotate the CDSs",
        logfile=options.logfile,
        output_file='CCDS_transcripts.txt',
        from_scratch=options.from_scratch)

    # TODO: Transcripts browser


    # Output the desired file
    logging.info("STEP  9:  Saving outputs...")

    json_output=compute_json(ccds_file="CCDS_transcripts.txt",
                             variant_file="VariantGTF.txt",
                             logfile=options.logfile,
                             output_file=options.output_filename,
                             from_scratch=options.from_scratch,
                             pas_tolerance=options.pas_tolerance)

    if options.gtf_filename:
        json2gtf(options.output_filename, options.gtf_filename, options.genome_filename,
                 options.gene, False)
    if options.extended_gtf_filename:
        logging.debug("""WARNING: you are creating a file that is not consistent with the GTF specifications.
        See http://mblab.wustl.edu/GTF22.html""")
        json2gtf(options.output_filename, options.extended_gtf_filename, options.genome_filename,
                 options.gene, True)


    # Clean mess
    logging.info("STEP 10:  Finalizing...")

    if options.compress:
        exec_system_command("gzip -q9 " + options.output_filename + " "+ options.logfile,
                            error_comment="Could not compress final files",
                            logfile="/dev/null",
        output_file=options.output_filename + '.gz')

    if not options.no_clean:
        tempfiles=("TEMP_COMPOSITION_TRANS1_1.txt", "TEMP_COMPOSITION_TRANS1_2.txt",
                   "TEMP_COMPOSITION_TRANS1_3.txt", "TEMP_COMPOSITION_TRANS1_4.txt",
                   "TRANSCRIPTS1_1.txt", "TRANSCRIPTS1_2.txt", "TRANSCRIPTS1_3.txt", "TRANSCRIPTS1_4.txt",
                   "VariantGTF.txt", "build-ests.txt", "CCDS_transcripts.txt", "config-dump.ini",
                   "genomic-exonforCCDS.txt", "info-pid-*.log", "isoforms.txt", "gene-struct.txt",
                   "meg-edges.txt", "megs.txt", "out-after-intron-agree.txt", "out-agree.txt", "out-fatt.txt",
                   "predicted-introns.txt", "processed-ests.txt", "processed-megs-info.txt",
                   "processed-megs.txt", "raw-multifasta-out.txt", "time-limits")
        subprocess.call("rm -f " + " ".join(tempfiles), shell=True)


def prepare_loggers(options):
    """Prepare loggers.

    Save DEBUG and higher messages to options.logfile, and INFO and higher messages to stdout
    Code adapted from
    http://docs.python.org/py3k/library/logging.html?highlight=logging#logging-to-multiple-destinations
    """
    logging.basicConfig(filename=options.logfile,
                        filemode='w',
                        format='%(levelname)s:%(name)s:%(asctime)s%(msecs)d:%(message)s',
                        datefmt='%Y%m%d-%H%M%S',
                        level=logging.DEBUG)
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(levelname)-8s] %(asctime)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


if __name__ == '__main__':


    try:
        options = parse_command_line()
        prepare_loggers(options)
        pintron_pipeline(options)
    except PIntronError as err:
        logging.exception("*** Fatal error caught during the execution of the pipeline! ***\n"
                          "%s", err)

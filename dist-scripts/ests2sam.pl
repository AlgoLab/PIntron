#!/usr/bin/env perl
####
#
#
#                              PIntron
#
# A novel pipeline for computational gene-structure prediction based on
# spliced alignment of expressed sequences (ESTs and mRNAs).
#
# Copyright (C) 2012 Gianluca Della Vedova and  Raffaella Rizzi
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

use strict;
use warnings;
use Data::Dumper;
use List::MoreUtils qw(each_array);
use Getopt::Long;

my $dir='.';
my $genome='genomic.txt';

my $result = GetOptions ("directory=s" => \$dir,
                         "genome=s" => \$genome);
$genome='genomic.txt' if ($genome eq '');
$dir='.' if ($dir eq '');

print STDERR "Analyzing directory $dir\n";
open INC, "<" ,$dir."/out-after-intron-agree.txt" or die "Can't open file. $!";
open SAM, ">" ,$dir."/est-alignments.sam" or die "Can't open file. $!";
die "Can't open file. $!" unless (-e $genome);
my $rname=`head -n 1 $genome`;
chomp $rname;
$rname=~s/^>//;
$rname=~s/[^!-()+-<>-~]//g;

my %ests=();
my $current_label='';
while (<INC>) {
    next if (/^\s*#/);
    
    if (/^>/) {
	# Header line
	my $current_line=$_;
	chomp $current_line;
	write_est($current_label, @{$ests{$current_label}}) unless ($current_label eq '' );
	$current_label=$current_line;
	# Sanitize name according to SAM specs
	$current_label=~s/^>\///;
	$current_label=~s/\/.*//;
	$current_label=~s/[^!-?A-~]//g;

	die "Gene $current_label appears more then once in $!\n" if exists $ests{$current_label};
	$ests{$current_label}=[];
	next;
    } 
    # EST alignment line
    my $factorization={};
    ($factorization->{relative_l}, $factorization->{relative_r},
     $factorization->{absolute_l}, $factorization->{absolute_r},
     $factorization->{est_factor}, $factorization->{genomic_factor})=
             split(/ /, $_, 6);
    push @{$ests{$current_label}}, $factorization;
}
write_est($current_label, @{$ests{$current_label}});
close INC;
close SAM;

sub write_est {
    my $est_label=shift;
    my @factors=@_;
    my $num_factors=scalar @factors;

    # SAM formats requires the next fragment for each read/fragment
    # @rotated array keeps a reference to the next fragment
    my @rotated=@factors;
    my $elem=shift @rotated;
    push @rotated, $elem;

    my @indices=(1 .. $num_factors);
    my $iter = each_array(@factors, @rotated, @indices);
    while ( my ($f, $next, $index) = $iter->() ) {
#	print Dumper($f);
	my $flag=2;
	$flag++ if (scalar @factors > 1);
	$flag+=128 if $index==1;
	$flag+=256 if $index==$num_factors;
	print SAM join "\t", ($est_label,       # QNAME
                          $flag, # FLAG
                          $rname, # RNAME
                          $f->{absolute_l}, # POS
                          '*', # MAPQ
                          '*', # CIGAR
                          '*', # RNEXT
                          $next->{absolute_l}, # PNEXT
                          0, # TLEN
                          $f->{est_factor}, # SEQ
                          '*', # QUAL
                             );
	print SAM "\n";
    }
}

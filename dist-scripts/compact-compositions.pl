#!/usr/bin/env perl
####
#
#
#                              PIntron
#
# A novel pipeline for computational gene-structure prediction based on
# spliced alignment of expressed sequences (ESTs and mRNAs).
#
# Copyright (C) 2010  Raffaella Rizzi
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

####
#
# Script che produce il file build-ests.txt (STDOUT) dal file
# out-after-intron-agree.txt (STDIN) (see PIntron pipeline).
# Questo script legge anche genomic.txt (genomic file)
# ATTENZIONE: Bisogna essere sicuri che ad ogni EST corrisponda una sola
# composizione.
# L'argomento della linea di comando e' un parametro per simulare la gestione degli
# esterni di ASPIC (dovrebbe essere 2 volte il fattore minimo usato in fattorizzazione
# Per ora e' hard-coded
# Lo script produce anche il file genomic-exonforCCDS.txt degli esoni da RefSeq
#
####

use Time::HiRes;
use strict;
use warnings;

my ($seconds, $microseconds) = Time::HiRes::gettimeofday();
print STDERR "Start ", $seconds, " seconds ", $microseconds, " microseconds\n";
my $t0=[$seconds, $microseconds];

my $min_ext=30;

open(GEN, "<genomic.txt") or die "Can't open the genomic file! $!\n";

my $gen_header=<GEN>;
close GEN;

chomp $gen_header;
$gen_header =~ m/^>chr([xXyY\d]+):(\d+):(\d+):([-+]{0,1}1)/i or die "Header in genomic.txt uncorrect!\n";

my $abs_left=($2 < $3)?($2):($3);
my $abs_right=($2 < $3)?($3):($2);
my $strand=$4;
my $boundary=0;

print $abs_left, "\n";
print $abs_right, "\n";
print $strand, "\n";
print $boundary, "\n";

my $gen_length=$abs_right-$abs_left+1;

undef $/;
my $file=<>;

#print $file;

my @file_str=split /^>/m, $file;

my %composition_hash=();    #key=GB id
my %left_exon_hash=();      #key=left
my %right_exon_hash=();     #key=right
my %polya_hash=();
my %polya_exon_hash=();

my $counter=1;

open(OUT, ">./genomic-exonforCCDS.txt") or die "Can't create RefSeq exon file! $!\n";

my %compact_composition=();
foreach my $composition(@file_str){
   if($composition ne ""){
        chomp $composition;
        my @compo_str=split /\n/, $composition;
        my $header=shift @compo_str;
        #print $header, "\n";
        
        my $polya=0;
        my $gb="";
        
        my $is_RefSeq=0;
        
        if($header =~ m/\/gb=(\w+)/){
            $gb=$1;
            if($gb =~ m/^NM_/){
                $is_RefSeq=1;
            }
         }
        else{
            die "No GB ID found for $header\n";
        }
        #print $gb, "\n";

        my @exon_list=();
        foreach my $row(@compo_str){
            if($row ne ""){
                #print $row, "\n";
                if($row =~ m/^\#/){
                    if($row =~ m/^\#polya=(\d+)/){
                        $polya=$1;
                    }
                }
                else{
                    if($row =~ m/\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\w+)/){
                        my $gen_left=$3;
                        my $gen_right=$4;
                        my $est_seq=$5;
                        my $gen_seq=$6;
                        
                        #$key_str.=$gen_left."-".$gen_right."-";
  
                        my @coord_list=($gen_left, $gen_right, $est_seq, $gen_seq);
                        push @exon_list, \@coord_list;
                    }
                    else{
                        die "Wrong format file!\n";
                    }
                }
            }
        }
      if(substr($gb, 0, 2) eq "NM"){
        foreach my $ref(@exon_list){
            my @coord_list=@{$ref};
            print OUT $coord_list[0], " ", $coord_list[1], " ", $coord_list[2], "\n";
        }
      }

        my $key_str="";
        my $key_must_not_exist=0;

        if(scalar(@exon_list) > 1){
            my @first_coord_list=@{$exon_list[0]};
            $key_str.=$first_coord_list[1]."-";
            for(my $i=1; $i < $#exon_list; $i++){
                 my @coord_list=@{$exon_list[$i]};
                 $key_str.=$coord_list[0]."-".$coord_list[1]."-";
            }
            my @last_coord_list=@{$exon_list[$#exon_list]};
            $key_str.=$last_coord_list[0]."-";
        
            if(substr($gb, 0, 2) eq "NM"){
                $key_str.=$gb;
                $key_must_not_exist=1;
            }
        }
        
        if($key_str ne "" && exists $compact_composition{$key_str}){
            if($key_must_not_exist == 1){
                die "Failure 1!\n";
            }
            
            my @gb_id=@{$compact_composition{$key_str}};
            my $j=0;
            my $stop=0;
            while($j <= $#gb_id && $stop == 0){
                my $id=$gb_id[$j];
                if(!exists $composition_hash{$id}){
                    die "Failure 3!\n";
                }
                my @temp_list=@{$composition_hash{$id}};
                my $ESTs=shift @temp_list;
                if(scalar(@exon_list) != scalar(@temp_list)){
                    die "Failure 4!\n";
                }
                my @first_coords=@{$temp_list[0]};
                my @last_coords=@{$temp_list[$#temp_list]};
                my @add_first_coords=@{$exon_list[0]};
                my @add_last_coords=@{$exon_list[$#exon_list]};
                if($first_coords[1] != $add_first_coords[1] || $last_coords[0] != $add_last_coords[0]){
                      die "Failure 5!\n";
                }
                my $ok=0;
                my $new_last_right=-1;
                my $new_last_est_seq;
                my $new_last_gen_seq;
                if($polya == 1){
                    if($polya_hash{$id} == 1){
                        if($last_coords[1] == $add_last_coords[1]){
                            $new_last_right=$last_coords[1];
                            $new_last_est_seq=$last_coords[2];
                            $new_last_gen_seq=$last_coords[3];
                            $ok=1;
                        }
                    }
                    else{
                        if($last_coords[1] <= $add_last_coords[1]){
                            $new_last_right=$add_last_coords[1];
                            $new_last_est_seq=$add_last_coords[2];
                            $new_last_gen_seq=$add_last_coords[3];
                            $ok=1;
                        }
                    }
                }
                else{
                    if($polya_hash{$id} == 1){
                        if($last_coords[1] >= $add_last_coords[1]){
                            $new_last_right=$last_coords[1];
                            $new_last_est_seq=$last_coords[2];
                            $new_last_gen_seq=$last_coords[3];
                            $ok=1;
                        }
                    }
                    else{
                        $new_last_right=($last_coords[1] >= $add_last_coords[1])?($last_coords[1]):($add_last_coords[1]);
                        $new_last_est_seq=($last_coords[1] >= $add_last_coords[1])?($last_coords[2]):($add_last_coords[2]);
                        $new_last_gen_seq=($last_coords[1] >= $add_last_coords[1])?($last_coords[3]):($add_last_coords[3]);
                        $ok=1;
                    }
                }
                
                if($ok){
                    my $new_first_left=($first_coords[0] <= $add_first_coords[0])?($first_coords[0]):($add_first_coords[0]);
                    my $new_first_est_seq=($first_coords[0] <= $add_first_coords[0])?($first_coords[2]):($add_first_coords[2]);
                    my $new_first_gen_seq=($first_coords[0] <= $add_first_coords[0])?($first_coords[3]):($add_first_coords[3]);
  
                    #print "Updated composition!\n";
                    #print "\tGB: $id; polya: ".$polya_hash{$id}."; ESTs: ", $ESTs, "\n";
                    #for(my $k=0; $k <= $#temp_list; $k++){
                    #    my @coord_list=@{$temp_list[$k]};
                    #    print "\t", $coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
                    #}
                    
                    #print "\tWith composition!\n";
                    #print "\t\tGB: $gb; polya: ".$polya."\n";
                    #for(my $k=0; $k <= $#exon_list; $k++){
                    #    my @coord_list=@{$exon_list[$k]};
                    #    print "\t\t", $coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
                    #}                 
                    
                     if($polya_hash{$id} == 0){
                        $polya_hash{$id}=$polya;
                    }
                    
                    $first_coords[0]=$new_first_left;
                    $first_coords[2]=$new_first_est_seq;
                    $first_coords[3]=$new_first_gen_seq;
            if($first_coords[1]-$first_coords[0]+1 != length($first_coords[3])){
                die "Failure 13!\n";
            }
                  $temp_list[0]=\@first_coords;
                    
                    $last_coords[1]=$new_last_right;
                    $last_coords[2]=$new_last_est_seq;
                    $last_coords[3]=$new_last_gen_seq;
            if($last_coords[1]-$last_coords[0]+1 != length($last_coords[3])){
                die "Failure 13!\n";
            }


                    $temp_list[$#temp_list]=\@last_coords;
                    
                    $ESTs=$ESTs+1;                    
                    unshift @temp_list, $ESTs;
                    
                    $composition_hash{$id}=\@temp_list;
  
                    #print "\tTo composition!\n";
                    #print "\tGB: $id; polya: ".$polya_hash{$id}."; ESTs: ", $temp_list[0], "\n";
                    #for(my $k=1; $k <= $#temp_list; $k++){
                    #    my @coord_list=@{$temp_list[$k]};
                    #    print "\t\t", $coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
                    #}
                    #print "\n";
  
                    $stop=1;
                }
                $j++;
             }
             
             if($stop == 0){
                 push @gb_id, $gb;
                 $compact_composition{$key_str}=\@gb_id;
                 unshift @exon_list, 1;
                 $composition_hash{$gb}=\@exon_list;
                 $polya_hash{$gb}=$polya;
                 
                 #print "Added new composition!\n";
                 #print "\tGB: $gb; polya: ".$polya_hash{$gb}."; ESTs: ", $exon_list[0], "\n";
                 #for(my $k=1; $k <= $#exon_list; $k++){
                 #   my @coord_list=@{$exon_list[$k]};
                 #   print "\t", $coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
                 #}
                 #print "\n";
             }
        }
        else{
            if($key_str ne ""){
                my @gb_id=($gb);
                $compact_composition{$key_str}=\@gb_id;
                if(exists $composition_hash{$gb}){
                     die "Failure 2!\n";
                }
            }
            unshift @exon_list, 1;
            $composition_hash{$gb}=\@exon_list;
            $polya_hash{$gb}=$polya;
            
            #print "Added new composition!\n";
            #print "\tGB: $gb; polya: ".$polya_hash{$gb}."; ESTs: ", $exon_list[0], "\n";
            #for(my $k=1; $k <= $#exon_list; $k++){
            #    my @coord_list=@{$exon_list[$k]};
            #    print "\t", $coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
            #}
            #print "\n";
         }
    }
}

close(OUT);

#print scalar(keys %composition_hash), "\n";

my @keys;

#@keys=keys %composition_hash;
#foreach my $key(@keys){
#    print $key, "\n";
#    my @temp_list=@{$composition_hash{$key}};
#    my $ESTs=shift @temp_list;
#    print "\t".$ESTs."(".$polya_hash{$key}.")", "\n";
#    foreach my $ref(@temp_list){
#        my @coord_list=@{$ref};
#        if($coord_list[1]-$coord_list[0]+1 != length($coord_list[3])){
#            die "Failure 6!\n";
#        }
#        print "\t\t".$coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
#    }
#}
#exit;

@keys=keys %composition_hash;
foreach my $key(@keys){
    #print $key, "\n";
    my @temp_list=@{$composition_hash{$key}};
    my $ESTs=shift @temp_list;
    #print "\t".$ESTs."(".$polya_hash{$key}.")", "\n";
    my $i=0;
    #foreach my $ref(@temp_list){
    while($i <= $#temp_list){
        my @coord_list=@{$temp_list[$i]};
        #print "\t\t".$coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
        my $exon_key=$coord_list[0]."-".$coord_list[1];

        if($i == $#temp_list){
            if(!exists $polya_exon_hash{$exon_key} || $polya_exon_hash{$exon_key} == 0){
                $polya_exon_hash{$exon_key}=$polya_hash{$key};
            }
        }
        else{
            $polya_exon_hash{$exon_key}=0;
        }
        
        my %left_temp_hash;
        if(exists $left_exon_hash{$coord_list[0]}){
            %left_temp_hash=%{$left_exon_hash{$coord_list[0]}};
        }
        else{
           %left_temp_hash=();
        }
        my @left_id_list;
        if(exists $left_temp_hash{$coord_list[1]}){
           @left_id_list=@{$left_temp_hash{$coord_list[1]}};
        }
        else{
           @left_id_list=();
        }
        push @left_id_list, $key;
        $left_temp_hash{$coord_list[1]}=\@left_id_list;
        $left_exon_hash{$coord_list[0]}=\%left_temp_hash;

        my %right_temp_hash;
        if(exists $right_exon_hash{$coord_list[1]}){
            %right_temp_hash=%{$right_exon_hash{$coord_list[1]}};
        }
        else{
           %right_temp_hash=();
        }
        my @right_id_list;
        if(exists $right_temp_hash{$coord_list[0]}){
           @right_id_list=@{$right_temp_hash{$coord_list[0]}};
        }
        else{
           @right_id_list=();
        }
        push @right_id_list, $key;
        $right_temp_hash{$coord_list[0]}=\@right_id_list;
        $right_exon_hash{$coord_list[1]}=\%right_temp_hash;
        
        $i++;
    }
}

#@keys=keys %composition_hash;
#foreach my $key(@keys){
#    print $key, "\n";
#    my @temp_list=@{$composition_hash{$key}};
#    my $ESTs=shift @temp_list;
#    print "\t".$ESTs."(".$polya_hash{$key}.")", "\n";
#    foreach my $ref(@temp_list){
#        my @coord_list=@{$ref};
#        if($coord_list[1]-$coord_list[0]+1 != length($coord_list[3])){
#            die "Failure 6!\n";
#        }
#        print "\t\t".$coord_list[0]."-".$coord_list[1]." ".$coord_list[2]." ".$coord_list[3]."\n";
#    }
#}
#exit;

#my @keysp=keys %left_exon_hash;
#foreach my $key(@keysp){
#    print $key, "\n";
#    my %temp_hash=%{$left_exon_hash{$key}};
#    my @keysi=keys %temp_hash;
#    foreach my $keyi(@keysi){
#        my @temp_list=@{$temp_hash{$keyi}};
#        print "\t", $keyi, " ", join(",", @temp_list), "\n";
#    }
#}

#my @keysp=keys %right_exon_hash;
#foreach my $key(@keysp){
#    print $key, "\n";
#    my %temp_hash=%{$right_exon_hash{$key}};
#    my @keysi=keys %temp_hash;
#    foreach my $keyi(@keysi){
#        my @temp_list=@{$temp_hash{$keyi}};
#        print "\t", $keyi, " ", join(",", @temp_list), "\n";
#    }
#}

@keys=keys %composition_hash;
foreach my $key(@keys){
    #print $key, "\n";
    my @temp_list=@{$composition_hash{$key}};
    
    if(substr($key, 0, 2) ne "NM" && scalar(@temp_list) > 2){
    
        my @first_coords=@{$temp_list[1]};
        my @next_first_coord=@{$temp_list[2]};

        #print "Reduce ".$first_coords[0]."-".$first_coords[1]."\n";
        #print "\tnext ".$next_first_coord[0]."-".$next_first_coord[1]."\n";
    
        my %right_temp_hash;
        if(!exists $right_exon_hash{$first_coords[1]}){
            die "Failure 7!\n";
        }
        %right_temp_hash=%{$right_exon_hash{$first_coords[1]}};
        if(!exists $right_temp_hash{$first_coords[0]}){
            die "Failure 8!\n";
        }
        my @left_ordered_keys= sort {$a <=> $b} keys %right_temp_hash;
        
        #print "\t check ", join(",",@left_ordered_keys), "\n";
        
        my $i=0;
        my $stop=0;
        while($i <= $#left_ordered_keys && $stop == 0){
            if($left_ordered_keys[$i] == $first_coords[0]){
                $stop=1;
            }
            else{
                my @id_list=@{$right_temp_hash{$left_ordered_keys[$i]}};
                my $j=0;
                my $stop2=0;
                while($j <= $#id_list && $stop2 == 0){
                    my @exon_list=@{$composition_hash{$id_list[$j]}};
                    
                    my $k=1;
                                        
                    my $stop3=0;
                    #while($k <= $#exon_list && $stop3 == 0){
                    while($k < $#exon_list && $stop3 == 0){
                        my @coord_list=@{$exon_list[$k]};
                        if($coord_list[0] == $left_ordered_keys[$i] && $coord_list[1] == $first_coords[1]){
                            $stop3=1;
                        }
                        else{                    
                            $k++;
                        }
                    }
                    #if($stop3 != 1){
                    #    die "Failure 11!\n";
                    #}
                
                    #if($k < $#exon_list){
                    if($stop3 == 1){
                        #my @next_coord_list=@{$exon_list[$k+1]};
                        #if($next_coord_list[0] == $next_first_coord[0]){
                            $stop2=1;
                            my @coord_list=@{$exon_list[$k]};
                            $first_coords[0]=$coord_list[0];
                            $first_coords[1]=$coord_list[1];
                            $first_coords[2]=$coord_list[2];
                            $first_coords[3]=$coord_list[3];
                            if($first_coords[1]-$first_coords[0]+1 != length($first_coords[3])){
                                die "Failure 13!\n";
                            }
                            #print "\treduced to ".$first_coords[0]."-".$first_coords[1]."\n";
                            $temp_list[1]=\@first_coords;
                        #}
                    }
                    
                    $j++;
                }
                
                $stop=$stop2;
            
                $i++;
            }
        }
  
        if($polya_hash{$key} == 0){  
            my @last_coords=@{$temp_list[$#temp_list]};
            my @prev_last_coord=@{$temp_list[$#temp_list-1]};
            
            #print "Reduce ".$last_coords[0]."-".$last_coords[1]."\n";
            #print "\tprevious ".$prev_last_coord[0]."-".$prev_last_coord[1]."\n";
            
            my %left_temp_hash;
            if(!exists $left_exon_hash{$last_coords[0]}){
                die "Failure 9!\n";
            }
            %left_temp_hash=%{$left_exon_hash{$last_coords[0]}};
            if(!exists $left_temp_hash{$last_coords[1]}){
                die "Failure 10!\n";
            }
            my @right_ordered_keys= sort {$b <=> $a} keys %left_temp_hash;
            
            #print "\t", join(",",@right_ordered_keys), "\n";
            
            $i=0;
            $stop=0;
            while($i <= $#right_ordered_keys && $stop == 0){
                if($right_ordered_keys[$i] == $last_coords[1]){
                    $stop=1;
                }
                else{
                    #print "\tCheck ", $right_ordered_keys[$i], "\n";
                    my @id_list=@{$left_temp_hash{$right_ordered_keys[$i]}};
                    my $j=0;
                    my $stop2=0;
                    while($j <= $#id_list && $stop2 == 0){
                        #print "\t\tID: ", $id_list[$j], "\n";
                        my @exon_list=@{$composition_hash{$id_list[$j]}};
                        
                        #my $k=1;
                        my $k=2;
                        
                        my $stop3=0;
                        
                        while($k <= $#exon_list && $stop3 == 0){
                            my @coord_list=@{$exon_list[$k]};
                            #print "\t\t", $coord_list[0], "-", $coord_list[1], "\n";
                            if($coord_list[0] == $last_coords[0] && $coord_list[1] == $right_ordered_keys[$i]){
                                $stop3=1;
                            }
                            else{                    
                                $k++;
                            }
                        }
                        #if($stop3 != 1){
                        #    die "Failure 12(2)!\n";
                        #}
 
                        #if($k > 1){
                        if($stop3 == 1){
                            #my @prev_coord_list=@{$exon_list[$k-1]};
                            #print "\tTo ".${$exon_list[$k]}[0]."-".${$exon_list[$k]}[1]."\n";
                            #print "\t\tpreviuos ".$prev_coord_list[0]."-".$prev_coord_list[1]."\n";
 
                            #if($prev_coord_list[1] == $prev_last_coord[1]){
                                 $stop2=1;
                                 my @coord_list=@{$exon_list[$k]};
                                 $last_coords[0]=$coord_list[0];
                                 $last_coords[1]=$coord_list[1];
                                 $last_coords[2]=$coord_list[2];
                                 $last_coords[3]=$coord_list[3];
                                 if($last_coords[1]-$last_coords[0]+1 != length($last_coords[3])){
                                    die "Failure 14!\n";
                                 }

                                 $temp_list[$#temp_list]=\@last_coords;
                                 $polya_hash{$key}=$polya_hash{$id_list[$j]};
                             #}
                        }
            
                        $j++;
                    }
                
                    $stop=$stop2;
          
                    $i++;
                }
            }
        }
    }
    $composition_hash{$key}=\@temp_list;    
}

my %print_compositions=();

my @print_exon_list=();
my %ordered_print_exon_hash=();

my %print_exon_hash=();
my @print_exon_seq_list=();
 
my $exon_index=0;

@keys=keys %composition_hash;
my $min_left=$gen_length+1;
my $max_right=0;

foreach my $key(@keys){
    #print $key, "\n";
    my @temp_list=@{$composition_hash{$key}};
    my $ESTs=shift @temp_list;
    
    my $composition_str="";
  
    my $isRefSeq=0;
    if(substr($key, 0, 2) eq "NM"){
        $isRefSeq=1;
    }
     
    foreach my $coord_ref(@temp_list){
        my @coord_list=@{$coord_ref};
        
        if($max_right < $coord_list[1]){
            $max_right=$coord_list[1];
        }
 
         if($min_left > $coord_list[0]){
            $min_left=$coord_list[0];
        }
        
        my $exon_key=$coord_list[0]."-".$coord_list[1];
        my $polya=$polya_exon_hash{$exon_key};
        
        if($isRefSeq == 1){
            $exon_key.=":".$key;
        } 
        
        if(!exists $print_exon_hash{$exon_key}){
             $print_exon_hash{$exon_key}=$exon_index;
             $exon_index++;
             
             my $print_str=$coord_list[0].":".$coord_list[1].":".$polya;

         #print $print_str, " ", $exon_key, "\n";

             push @print_exon_list, $print_str;
             
             my $add_seq=($isRefSeq == 1)?($coord_list[2]):($coord_list[3]);
             push @print_exon_seq_list, $add_seq;

          if(exists $ordered_print_exon_hash{$coord_list[0]}){
            #print "\tExist key ", $coord_list[0], "\n";
            my %temp_hash=%{$ordered_print_exon_hash{$coord_list[0]}};
            if(exists $temp_hash{$coord_list[1]}){
                #print "\t\tExist key ", $coord_list[1], "\n";
                my @temp_list=($polya, $print_exon_hash{$exon_key}, $add_seq);
                my @add_list=@{$temp_hash{$coord_list[1]}};
                push @add_list, \@temp_list;
                $temp_hash{$coord_list[1]}=\@add_list;
                $ordered_print_exon_hash{$coord_list[0]}=\%temp_hash;
            }
            else{
                my @temp_list=($polya, $print_exon_hash{$exon_key}, $add_seq);
                my @add_list=();
                push @add_list, \@temp_list;
                $temp_hash{$coord_list[1]}=\@add_list;
                $ordered_print_exon_hash{$coord_list[0]}=\%temp_hash;
            }
          }
          else{
            #print "\tNot Exist key ", $coord_list[0], "\n";
            my %temp_hash=();
            my @temp_list=($polya, $print_exon_hash{$exon_key}, $add_seq);
            my @add_list=();
            push @add_list, \@temp_list;
            $temp_hash{$coord_list[1]}=\@add_list;
            $ordered_print_exon_hash{$coord_list[0]}=\%temp_hash;
        }
        }
        
        my $index=$print_exon_hash{$exon_key};
        $composition_str.=$index.".";
    }
 
    chop $composition_str;
    
    if(exists $print_compositions{$composition_str}){
        if($isRefSeq == 1){
            die "Failure 12(1)!\n";
        }
        my @temp_list=@{$print_compositions{$composition_str}};
        my $counter=$temp_list[0];
        $counter+=$ESTs;
        $temp_list[0]=$counter;
        $print_compositions{$composition_str}=\@temp_list;
    }
    else{
        my @temp_list=($ESTs);
        if($isRefSeq == 1){
            push @temp_list, $key;
        }
        $print_compositions{$composition_str}=\@temp_list;
    }
}

#@keys=sort {$a <= $b} keys %ordered_print_exon_hash;
#foreach my $key(@keys){
#   print $key, "\n";
#   my %temp_hash=%{$ordered_print_exon_hash{$key}};    
#   my @keys2=sort {$a <= $b} keys %temp_hash;
#   foreach my $key2(@keys2){
#       print "\t", $key2, "\n";
#       my @temp_list=@{$temp_hash{$key2}};
#       foreach my $ref(@temp_list){
#           my @temp_list2=@{$ref};
#           print "\t\t", $temp_list2[1], "\n";
#       }
#   }
#}
#exit;

my @keys_compositions=keys %print_compositions;
print scalar(@keys_compositions), "\n"; #number of compositions
print scalar(@print_exon_list), "\n";   #number of exons
my $coverage_length=$max_right;
print $coverage_length, "\n";

#print join("\n", @print_exon_list);

@keys=sort {$a <=> $b} keys %ordered_print_exon_hash;
my $ordered_index=0;
my %hash_map=();
foreach my $key(@keys){
    #print $key, "\n";
    my %temp_hash=%{$ordered_print_exon_hash{$key}};    
    my @keys2=sort {$a <=> $b} keys %temp_hash;
    foreach my $key2(@keys2){
        #print "\t", $key2, "\n";
        my @temp_list=@{$temp_hash{$key2}};

        foreach my $ref(@temp_list){
            my @temp_list2=@{$ref};
            if(exists $hash_map{$temp_list2[1]}){
                die "Failure 14!\n";
            }
            print $key.":".$key2.":".$temp_list2[0]."\n";
            #$hash_map{$temp_list2[1]}=[$ordered_index, $temp_list2[2]];
            $hash_map{$temp_list2[1]}=$ordered_index;
            $ordered_index++;
        }
    }
}
if($ordered_index != scalar(@print_exon_list)){
    die "Failure 15!\n";
}

foreach my $key(@keys_compositions){
    my @temp_list=@{$print_compositions{$key}};
    my $header="";
    foreach my $str(@temp_list){
        $header.=".$str";
    }
    print $header."\n";
    
    my @index_list=split /\./, $key;

    #print $key."\n";
    my $ordered_key="";
    foreach my $index(@index_list){
        $ordered_key.=$hash_map{$index}.".";
    }
    chop $ordered_key;
    print $ordered_key."\n";

    foreach my $index(@index_list){
        print $print_exon_seq_list[$index]."\n";
    }
}

print "\#\n*\n";

my $elapsed = Time::HiRes::tv_interval($t0);
print STDERR "Time elapsed: ", ($elapsed*1000000), " microsec\n";

exit;

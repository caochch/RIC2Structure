#!/usr/bin/perl
die "perl $0 virus.final.CON domains_bedpe \n" if(@ARGV != 2);
my $phase1_final_constrain_file=shift;
my $domains_bedpe=shift;

#public parameters
my $ref_seq_fasta_file="WuHan_Hu_1.fasta";      #the file name of target sequence, in fasta format
my $whole_RNA_length=29903;              #the length of target sequence
my $RICseq_count_matrix_for_correlation="../1.Loop_by_HiCCUPS/virus.onlyPart.resolution5.format.matrix";          #matrix used for correlation analysis
my $ruler_used=5;
my $window_len=1000;
my $max_basepair_length_for_whole=750;
#public parameters over

##test parameters
my $print_paired_pixels=0;
my $print_structure_pvalue=1;
##test parameter over

#loop from enriched pixels
my %all_loops;
open(CPB,$domains_bedpe) || die;
while(my $line=<CPB>){
        if($line=~/^#/){
                next;
        }
        chomp $line;
        my @sub=split/\s+/,$line;
        my $x_start=$sub[1];
        my $x_end=$sub[2];
        my $y_start=$sub[4];
        my $y_end=$sub[5];
        my $name=$sub[6];

        my @four=($x_start,$x_end,$y_start,$y_end);
        @four=sort {$a<=>$b} @four;

        $raw_loops{$name}=$four[0]."\t".$four[3];
        $four[0]=$four[0] < $extended_slop ? 0 : $four[0]-$extended_slop;
        $four[3]=$four[3] > $whole_RNA_length-$extended_slop ? $whole_RNA_length : $four[3]+$extended_slop;
        $all_loops{$name}=$four[0]."\t".$four[3];
}
#finished loop reading

#read RIC-seq count matrix
my %RICseq_matrix;
open(MT,$RICseq_count_matrix_for_correlation) || die;
my $head=<MT>;
chomp $head;
my @arr_head=split/\s+/,$head;
while(my $line=<MT>){
        chomp $line;
        my @sub=split/\s+/,$line;
        foreach (1..$#sub){
                if($arr_head[$_] <= $sub[0]){
                        next;
                }
                $RICseq_matrix{$sub[0]}{$arr_head[$_]}=$sub[$_];
        }
}

#read genome sequence
my $genome_seq;
open(FA,$ref_seq_fasta_file) || die;
while(my $line=<FA>){
        chomp $line;
        if($line=~/>/){
                next;
        }
        else{
                $genome_seq.=$line;
        }
}
my @genome_letter=split//,$genome_seq;
#finished genome sequence reading

#initialize dynamic solved
my %dynamic_constrains;
my %dynamic_solved;
my %dynamic_basepairs;
my %dynamic_unsolved;

foreach (1..$whole_RNA_length){
	$dynamic_unsolved{$_}=1;
}

open(PFCF,$phase1_final_constrain_file) || die;
while(my $line=<PFCF>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($line=~/-1/){
		next;
	}
	elsif($line=~/(\d+)\s+(\d+)/){
		my @pair=($1,$2);
		@pair=sort {$a<=>$b} @pair;
		$dynamic_constrains{$pair[0]."\t".$pair[1]}=$pair[0];
		$dynamic_basepairs{$pair[0]}=$pair[1];	#1-based
		$dynamic_basepairs{$pair[1]}=$pair[0];	#1-based
		$dynamic_solved{$pair[0]}=1;		#1-based
		$dynamic_solved{$pair[1]}=1;		#1-based
		foreach ($pair[0]..$pair[1]){
			$dynamic_unsolved{$_}=0;
		}
	}
}
#finished initialize dynamic solved 

#sort domains
my %loop_solved_ratio;
foreach my $id (keys %all_loops){
	my ($start,$end)=split/\s+/,$all_loops{$id};
	my $resolved_ratio=count_solved_in_domain($start,$end);
	$loop_solved_ratio{$id}=$resolved_ratio;
}

my @running_loops;
foreach (sort {$loop_solved_ratio{$b} <=> $loop_solved_ratio{$a}} keys %loop_solved_ratio){
	push (@running_loops,$_);
}
#sort domains over;

#
my %dynamic_single;
my $run_cycle;
#start to iteration
foreach my $id (@running_loops){
	my ($loop_start,$loop_end)=split/\s+/,$all_loops{$id};

	$run_cycle++;
	print "-----------running cycle $run_cycle----------------------\n";
	print $loop_start,"\t",$loop_end,"\t",$loop_solved_ratio{$id},"\tselected loop\n";

	my ($extended_loop_start,$extended_loop_end)=extended_to_farest_bp($loop_start,$loop_end);
	$extended_loop_end=$extended_loop_end+1;
	my $extended_loop_len=$extended_loop_end-$extended_loop_start;
	my $extended_loop_seq=substr($genome_seq,$extended_loop_start,$extended_loop_len);

        open(LOOPS,">tmp.cycle$run_cycle.fa") || die;
        print LOOPS ">cycle$run_cycle\n$extended_loop_seq\n";
        close LOOPS;
        my $fa_file="tmp.cycle$run_cycle.fa";
        my $ct_file="tmp.cycle$run_cycle.ct";

	#add constrains
        my $have_cons=collect_constrain($extended_loop_start,$extended_loop_end);

	print $extended_loop_start,"\t",$extended_loop_end,"\t",$have_cons,"\tloop basic\n";

        #run prediction
        if($have_cons){
		$have_cons = $have_cons < $max_basepair_length_for_whole ? $max_basepair_length_for_whole : $have_cons;
		print "Command: /mnt/pub/work/caochch/project_My/project030_RIC_Covid19Virus/software/RNAstructure/exe/Fold $fa_file $ct_file --constraint tmp.CON --maxdistance $have_cons\n";
                `/mnt/pub/work/caochch/project_My/project030_RIC_Covid19Virus/software/RNAstructure/exe/Fold $fa_file $ct_file --constraint tmp.CON --maxdistance $have_cons`;
        }
        else{
                #`/mnt/pub/work/caochch/project_My/project030_RIC_Covid19Virus/software/RNAstructure/exe/Fold $fa_file $ct_file --maxdistance 750`;
		warn "this will not happen\n";
		die;
        }
        if(-e $ct_file){
        }
        else{
                print "No available structure\n";
		die;
        }
        my ($index,$energy,$pvalue,$ctinfo)=select_by_correlation($extended_loop_start,$extended_loop_end,$ct_file,$extended_loop_start,$extended_loop_end);
        print $index,"\t",$energy,"\t",$pvalue,"\tselected\n";
        #`rm -rf $fa_file $ct_file`;
        if($pvalue > 0.05){
		warn "should check this region carefully\n";
		die;
        }
        if($pvalue eq "NA"){
		warn "should check this region carefully\n";
		die;
        }

        my @ct_lines=split/\n/,$ctinfo;
        shift @ct_lines;
	my $left_most_paired_in_domain;
	my $right_most_paired_in_domain;
        foreach my $l (@ct_lines){
                $l=~s/^\s+//;
                my @sub=split/\s+/,$l;
		if($#sub != 5){
			next;
		}
               	my $real_loci_left=$extended_loop_start+$sub[0];
             	my $real_loci_right=$extended_loop_start+$sub[4];
		if($sub[4]){	#paired
			if($real_loci_left >= $extended_loop_start and $real_loci_left <= $extended_loop_end and $real_loci_right >= $extended_loop_start and $real_loci_right <= $extended_loop_end){
	                 	my @pair=($real_loci_left,$real_loci_right);
	                        @pair=sort {$a<=>$b} @pair;
				if(exists $dynamic_basepairs{$pair[0]} and $dynamic_basepairs{$pair[0]} ne $pair[1]){
					warn "conflict basepair\n";
					die;
				}
				if(exists $dynamic_basepairs{$pair[1]} and $dynamic_basepairs{$pair[1]} ne $pair[0]){
					warn "conflict basepair\n";
					die;
				}
				$dynamic_constrains{$pair[0]."\t".$pair[1]}=$pair[0];
	                       	$dynamic_basepairs{$pair[0]}=$pair[1];
	                        $dynamic_basepairs{$pair[1]}=$pair[0];
				if(defined($left_most_paired_in_domain)){
					if($pair[0] < $left_most_paired_in_domain){
						$left_most_paired_in_domain=$pair[0];
					}
					if($pair[1] > $right_most_paired_in_domain){	
						$right_most_paired_in_domain=$pair[1];
					}
				}
				else{
					$left_most_paired_in_domain=$pair[0];
					$right_most_paired_in_domain=$pair[1];
				}
				$dynamic_solved{$pair[0]}=1;
				$dynamic_solved{$pair[1]}=1;
			}
		}
	}
	foreach my $l (@ct_lines){
		$l=~s/^\s+//;
		my @sub=split/\s+/,$l;
		if($#sub != 5){
			next;
		}
		my $real_loci_left=$extended_loop_start+$sub[0];
		if($sub[4]){
			next;
		}
		else{	#single
			if($real_loci_left >= $extended_loop_start and $real_loci_left <= $extended_loop_end and $real_loci_left >= $left_most_paired_in_domain and $real_loci_left <= $right_most_paired_in_domain){

				$dynamic_single{$real_loci_left}=1;
				$dynamic_solved{$real_loci_left}=1;
                	}
		}
        }

	print $left_most_paired_in_domain,"\t",$right_most_paired_in_domain,"\tNewSolvedRegion\n";

	foreach ($left_most_paired_in_domain..$right_most_paired_in_domain){
		$dynamic_unsolved{$_}=0;
	}

	my $tmp_prefix="phase4_cycle$run_cycle";
	visualize_by_IGV($tmp_prefix);
	visualize_by_ct($tmp_prefix);
	write_constrains($tmp_prefix);

        #update resolved blocks
        #open(RESOLVE,">log.resolved_bases.cycle$run_cycle.bed") || die;
        #my $resolved_bases_blocks=print_resolved_region();
        #print RESOLVE $resolved_bases_blocks;
        #close RESOLVE;
        #update resolved regions and print over
}

write_constrains("phase4_final");
visualize_by_dot("phase4_final");
visualize_by_IGV("phase4_final");
visualize_by_ct("phase4_final");

sub write_constrains{
        my $prefix=shift;
        open(WCON,">virus.$prefix.CON") || die;
        print WCON "DS:\n-1\n";
        print WCON "SS:\n";
        print WCON "-1\nMod:\n-1\nPairs:\n";
        foreach my $p (sort {$dynamic_constrains{$a} <=> $dynamic_constrains{$b}} keys %dynamic_constrains){
                my ($loci_left,$right_loci)=split/\s+/,$p;
                print WCON $loci_left," ",$right_loci,"\n";
        }
        print WCON "-1 -1\n";
        print WCON "FMN:\n-1\nForbids:\n-1 -1\n";
        close WCON;
}

sub print_resolved_region{
        my @all_bases=sort {$a<=>$b} keys %dynamic_solved;
        my $bedfiles;
        foreach my $i (0..$#all_bases){
                if($i == 0){
                        my $zerobased=$all_bases[$i]-1;
                        $bedfiles="NC_045512.2\t$zerobased\t";
                }
                elsif($i == $#all_bases){
                        my $tmp_loci=$all_bases[$i];
                        $bedfiles.="$tmp_loci\n";
                }
                else{
                        if($all_bases[$i] == $all_bases[$i-1]+1){
                                next;
                        }
                        else{
                                my $tmp_loci=$all_bases[$i-1];
                                $bedfiles.="$tmp_loci\n";
                                my $zerobased=$all_bases[$i]-1;
                                $bedfiles.="NC_045512.2\t$zerobased\t";
                        }
                }
        }
        return $bedfiles;
}


sub visualize_by_IGV{
        #visualize by IGV
        my $prefix=shift;
        open(IGV,">virus.$prefix.igv.bed") || die;
        print IGV "track graphType=arc\n";
        foreach (keys %dynamic_constrains){
                my ($loci_left,$loci_right)=split/\s+/,$_;
                $loci_left=$loci_left-1;
                $loci_right=$loci_right-1;
                print IGV "NC_045512.2\t$loci_left\t$loci_right\n";
        }
        close IGV;
}

sub visualize_by_ct{
        #visualize by ct
        my $prefix=shift;
        open(WCT,">virus.$prefix.whole.ct") || die;
        print WCT "$whole_RNA_length\tCovid19\n";
        my $offset=0;   #for the purpose of start at a given loci
        foreach my $i (1..$whole_RNA_length){
                print WCT $i-$offset,"\t";
                print WCT $genome_letter[$i-1],"\t";
                print WCT $i-1-$offset,"\t",$i+1-$offset,"\t";
                if(exists $dynamic_basepairs{$i}){
                        print WCT $dynamic_basepairs{$i}-$offset,"\t";
                }
                else{
                        print WCT "0\t";
                }
                print WCT $i-$offset,"\n";
        }
        close WCT;
}

sub visualize_by_ct{
        #visualize by ct
        my $prefix=shift;
        open(WCT,">virus.$prefix.whole.ct") || die;
        print WCT "$whole_RNA_length\tCovid19\n";
        my $offset=0;   #for the purpose of start at a given loci
        foreach my $i (1..$whole_RNA_length){
                print WCT $i-$offset,"\t";
                print WCT $genome_letter[$i-1],"\t";
                print WCT $i-1-$offset,"\t",$i+1-$offset,"\t";
                if(exists $dynamic_basepairs{$i}){
                        print WCT $dynamic_basepairs{$i}-$offset,"\t";
                }
                else{
                        print WCT "0\t";
                }
                print WCT $i-$offset,"\n";
        }
        close WCT;
}

sub visualize_by_dot{
        #visualize by dot
        my $prefix=shift;
        open(DOT,">virus.$prefix.dot") || die;
        my $string="." x $whole_RNA_length;
        foreach my $i (1..$whole_RNA_length){
                if(substr($string,$i-1,1) eq "."){
                        if(exists $dynamic_basepairs{$i}){
                                if($dynamic_basepairs{$i} < $i){
                                        die "$i: should be revised as ) corresponding to loci $dynamic_basepairs{$i}\n";
                                }
                                else{
                                        substr($string,$i-1,1)="(";
                                        substr($string,$dynamic_basepairs{$i}-1,1)=")";
                                }
                        }
                        else{
                                #single strand; nothing todo
                        }
                }
                else{
                        next;
                }

        }
        print DOT $genome_seq,"\n";
        print DOT $string,"\n";

        #my $test_len=$test_region_end-$test_region_start;
        #print DOT substr($genome_seq,$test_region_start,$test_len),"\n";       #test
        #print DOT substr($string,$test_region_start,$test_len),"\n";           #test
}

sub extended_to_farest_bp{
	my $start=shift;
	my $end=shift;
	$start=$start+1;
	my $left_most=$start;
	my $right_most=$end;

	print "From $left_most to $right_most before extending\n";

	#extend left boundary
	if(!$dynamic_unsolved{$start}){	#resolved, extend to left unreolved region at left duplex
		my $next_to_left=$left_most;
		while(1){
			my $next=$next_to_left-1;
			if($dynamic_unsolved{$next}){
				last;
			}
			if($next < 1){
				last;
			}
			if(exists $dynamic_basepairs{$next} and $dynamic_basepairs{$next} >= $start){
				$left_most=$next;
			}
			$next_to_left=$next_to_left-1;
		}
		print "Left is solved, from $start to $left_most after extend solved\n";
		while(1){
			my $next=$left_most-1;
			if(!$dynamic_unsolved{$next}){
				last;
			}
			if($next < 1){
				last;
			}
			$left_most=$next;
		}
		print "Left is solved, from $start to $left_most after extend unsolved left at duplex\n";
		if($left_most < $start - 1000){
			$left_most=$start;
			while(1){
				my $next=$left_most+1;
				if($dynamic_unsolved{$next}){
					last;
				}
				if($next > $end){
					last;
				}
				$left_most=$next;
			}
			print "Left is solved, but duplex is too long, back to shrink to unsolved, from $start to $left_most\n";
		}
		if(!$dynamic_unsolved{$left_most} and $left_most > 1){
			$left_most=$left_most-1;
		}
	}
	else{	#un resolved, extended unsolved
		while(1){
			my $next=$left_most-1;
			if(!$dynamic_unsolved{$next}){
				last;
			}
			if($next < 1){
				last;
			}
			$left_most=$next;
		}	
		print "Left is unsolved, from $start to $left_most after extending unsolved\n";
	}
	if(!$dynamic_unsolved{$end}){	#extend to right unreolved region at right duplex
		my $next_to_right=$right_most;
		while(1){
			my $next=$next_to_right+1;
			if($dynamic_unsolved{$next}){
				last;
			}
			if($next > $whole_RNA_length){
				last;
			}
			if(exists $dynamic_basepairs{$next} and $dynamic_basepairs{$next} <= $end){
				$right_most=$next;
			}
			$next_to_right=$next_to_right+1;
		}
		print "Right is solved, from $end to $right_most after extend solved\n";
		while(1){
			my $next=$right_most+1;
			if(!$dynamic_unsolved{$next}){
				last;
			}
			if($next > $whole_RNA_length){
				last;
			}
			$right_most=$next;
		}
		print "Right is solved, from $end to $right_most after extend unsolved right ad duplex\n";
		if($right_most > $end + 1000){#extended too long, we select to give up
			$right_most=$end;	#trace back, to shrink
			while(1){
				my $next=$right_most-1;
				if($dynamic_unsolved{$next}){
					last;
				}
				if($next < $start){
					last;
				}
				$right_most = $next;
			}
			print "Right is solved, but duplex is too long, back to shrink to unsolved, from $end to $right_most\n";
		}
		if(!$dynamic_unsolved{$right_most} and $right_most < $whole_RNA_length){
			$right_most=$right_most+1;
		}
			
	}
	else{	#unresolved, extended
		while(1){
			my $next=$right_most+1;
			if(!$dynamic_unsolved{$next}){
				last;
			}
			if($next > $whole_RNA_length){
				last;	
			}
			$right_most=$next;
		}
		print "Right is unsolved, Right from $end to $right_most after extending unsolved\n";
	}
	
	print "From $left_most to $right_most after extending\n";
        $left_most=$left_most-1;
        $right_most=$right_most-1;
	return ($left_most,$right_most);
}

sub count_solved_in_domain{
	my $start=shift;
	my $end=shift;
	my $already_solved;
	foreach my $j ($start+1..$end){
		if(exists $dynamic_solved{$j}){
			$already_solved++;
		}
	}
	my $len=$end-$start;
	my $ratio=$already_solved/$len;
	return $ratio;
}


sub collect_constrain{
	my $start=shift;
	my $end=shift;
	my %tmp_constrain_pairs;
	my %tmp_constrain_singlestrand;
	print $start,"\t",$end,"\tloop region to call constrain\n";
	foreach (keys %dynamic_constrains){
		my ($loci_left,$loci_right)=split/\s+/,$_;
		if($loci_left >= $start+1 and $loci_left <= $end and $loci_right >= $start+1 and $loci_right <= $end){
			$tmp_constrain_pairs{$loci_left}=$loci_right;
		}
		elsif($loci_left >= $start+1 and $loci_left <= $end and ($loci_right < $start+1 or $loci_right > $end)){
			$tmp_constrain_singlestrand{$loci_left}=1;
		}
		elsif(($loci_left < $start+1 or $loci_left > $end) and ($loci_right >= $start+1 and $loci_right <= $end)){
			$tmp_constrain_singlestrand{$loci_right}=1;
		}
		else{
			next;
		}
	}

	foreach my $loci (keys %dynamic_single){
		if($loci >= $start+1 and $loci <= $end){
			$tmp_constrain_singlestrand{$loci}=1;
		}
	}

	if(%tmp_constrain_pairs or %tmp_constrain_singlestrand){	
		my $max_basepair_distance;
		open(CON,">tmp.CON") || die;
		my $len=$end-$start;
		my $testString="." x $len; 
		print CON "DS:\n-1\n";
		print CON "SS:\n";
		foreach my $loci (sort {$a<=>$b} keys %tmp_constrain_singlestrand){
			my $relative_loci=$loci-$start;
			print CON $relative_loci,"\n";
		}
		print CON "-1\nMod:\n-1\nPairs:\n";
		foreach my $loci_left (sort {$a<=>$b} keys %tmp_constrain_pairs){
			my $loci_right=$tmp_constrain_pairs{$loci_left};
			my $this_basepair_distance=$loci_right-$loci_left;
			if($this_basepair_distance > $max_basepair_distance){
				$max_basepair_distance=$this_basepair_distance;
			}
			my $relative_loci_left=$loci_left-$start;
			my $relative_right_loci=$loci_right-$start;
			print CON $relative_loci_left," ",$relative_right_loci,"\n";
			substr($testString,$relative_loci_left-1,1)="(";
			substr($testString,$relative_right_loci-1,1)=")";
			
		}
		#print $testString,"\n";
		print CON "-1 -1\n";
		print CON "FMN:\n-1\nForbids:\n-1 -1\n";
		close CON;
		$max_basepair_distance=$max_basepair_distance+50;
		return $max_basepair_distance;
	}
	else{
		return 0;
	}
}

sub select_by_correlation{
	my $start=shift;
	my $end=shift;
	my $ct_file=shift;
	my $this_loop_start=shift;
	my $this_loop_end=shift;
	my $slop_due_to_RICseq_resolution=5;
	my $least_distance_to_diagonal=5;
	my %structure;
	my $id_of_structure;
	$/="ENERGY";
	open(CT,$ct_file) || die;
	<CT>;
	while(my $block=<CT>){
		chomp $block;
		$id_of_structure++;
		my %paired;
		my @lines=split/\n/,$block;
		my $first=shift @lines;
		$first=~s/^\s+//;
		my $energy=(split/\s+/,$first)[1];
		foreach my $l (@lines){
			$l=~s/^\s+//;
			my @sub=split/\s+/,$l;
			if($sub[4]){#paired
				$paired{$sub[0]}{$sub[4]}=1;
			}
			else{#single strand
			}
		}
		
		#do t-test
		my %paired_RICseq_signal;
		foreach my $i ($start..$end){
			my $win_i=int($i/$ruler_used)*$ruler_used;
			foreach my $j ($start..$end){
				my $win_j=int($j/$ruler_used)*$ruler_used;
				if($win_j <= $win_i){
					next;
				}
				else{
					my $relative_loci_i=$i-$start+1;
					my $relative_loci_j=$j-$start+1;
					if($paired{$relative_loci_i}{$relative_loci_j}){
						#print $ct_file,"\t",$id_of_structure,"\t",$relative_loci_i,"\t",$relative_loci_j,"\t",$i,"\t",$j,"\tpaired\n";
						$paired_RICseq_signal{$win_i."\t".$win_j}=1;
					}
				}
			}
		}
	
		#print $start,"\t",$end,"\n";
		if($print_paired_pixels){
			print ">$ct_file"."_"."$id_of_structure\n";
		}
	
		my @paired;
		my @other;
		my %unique_windows;
		foreach my $i ($start..$end-1){
			my $win_i=int($i/$ruler_used)*$ruler_used;
			foreach my $j ($start..$end-1){
				my $win_j=int($j/$ruler_used)*$ruler_used;
				if($win_j <= $win_i+$least_distance_to_diagonal){
					next;
				}
				if($unique_windows{$win_i."\t".$win_j}){	#already assigned
					next;
				}
				if($resolved_bases{$i} and $resolved_bases{$j}){	#keep focus on newly added region
					next;
				}
				if($i < $this_loop_start or $i > $this_loop_end or $j < $this_loop_start or $j > $this_loop_end){	#keep focus on newly added by this loop
					next;
				}
				$unique_windows{$win_i."\t".$win_j}=1;
				my $signal=$RICseq_matrix{$win_i}{$win_j};
				if(defined($signal)){
				}
				else{
					warn $ct_file,"\t",$i,"\t",$j,"\t",$win_i,"\t",$win_j,"\tbug\n";
					die;
				}
				if(exists $paired_RICseq_signal{$win_i."\t".$win_j}){
					if($print_paired_pixels){
						print "hs2\t$win_i\t",$win_i+$ruler_used,"\ths2\t$win_j\t",$win_j+$ruler_used,"\t0,0,255\n";
					}
					push (@paired,$signal);
				}
				else{
					my $distance_to_paired=2*$slop_due_to_RICseq_resolution;
					foreach my $paired_pixels (keys %paired_RICseq_signal){
						my ($paired_win_i,$paired_win_j)=split/\s+/,$paired_pixels;
						my $tmp_dis;
						if(abs($win_i-$paired_win_i) > abs($win_j-$paired_win_j)){
							$tmp_dis=abs($win_i-$paired_win_i);
						}
						else{
							$tmp_dis=abs($win_j-$paired_win_j);
						}
						if($tmp_dis < $distance_to_paired){
							$distance_to_paired=$tmp_dis;
						}
					}
					if($distance_to_paired <= $slop_due_to_RICseq_resolution){
						if($print_paired_pixels){
							print "hs2\t$win_i\t",$win_i+$ruler_used,"\ths2\t$win_j\t",$win_j+$ruler_used,"\t0,255,255\n";
						}
						push (@paired,$signal);
					}
					else{
						if($print_paired_pixels){
							print "hs2\t$win_i\t",$win_i+$ruler_used,"\ths2\t$win_j\t",$win_j+$ruler_used,"\t0,255,0\n";
						}
						push (@other,$signal);
					}
				}
			}
		}
		
		open(TMPP,">tmp.paired.RICseq.signal.list") || die;
		open(TMPO,">tmp.other.RICseq.signal.list") || die;
		foreach (@paired){
			print TMPP $_,"\n";
		}
		foreach (@other){
			print TMPO $_,"\n";
		}
		close TMPP;
		close TMPO;

		if(@paired < 1 or @other < 1){
			next;
		}

		my $test_result=`Rscript 0.Rscript_ttest.r`;
		$test_result=~s/\s+$//;
		my ($pvalue,$mean_paired,$median_paired,$mean_other,$median_other)=split/\s+/,$test_result;

		my $paired_sum=$mean_paired*($#paired+1);
		my $other_sum=$mean_other*($#other+1);

		$structure{$id_of_structure}{"energy"}=$energy;
		$structure{$id_of_structure}{"sum"}=$paired_sum;
		$structure{$id_of_structure}{"pvalue"}=$pvalue;
		$structure{$id_of_structure}{"ct"}=$block;
		if($print_structure_pvalue){
			print ">",$id_of_structure,"\t",$#paired+1,"\t",$paired_sum,"\t",$mean_paired,"\t",$median_paired,"\t",$#other+1,"\t",$other_sum,"\t",$mean_other,"\t",$median_other,"\t",$pvalue,"\ttestPvalue\n";
		}
	}
	$/="\n";

        if(%structure){
                #my @sorted_strcuture=sort {$structure{$a}{"pvalue"} <=> $structure{$b}{"pvalue"}} keys %structure;     #min pvalue
                my @sorted_strcuture=sort {$structure{$b}{"sum"} <=> $structure{$a}{"sum"}} keys %structure;            #max sum
                #return ($sorted_strcuture[0],$structure{$sorted_strcuture[0]}{"energy"},$structure{$sorted_strcuture[0]}{"pvalue"},$structure{$sorted_strcuture[0]}{"ct"});

                my $optimal_value=$structure{$sorted_strcuture[0]}{"sum"};

                my @candidate_structure;
                foreach (@sorted_strcuture){
                        if($structure{$_}{"sum"} >= $optimal_value){
                                push (@candidate_structure,$_);
                        }
                }

                my @sorted_again_structure=sort {$structure{$a}{"energy"} <=> $structure{$b}{"energy"}} @candidate_structure;        #lowest energy
                return ($sorted_again_structure[0],$structure{$sorted_again_structure[0]}{"energy"},$structure{$sorted_again_structure[0]}{"pvalue"},$structure{$sorted_again_structure[0]}{"ct"});

        }
        else{
                return ("NA","NA","NA","NA");
        }
}

	


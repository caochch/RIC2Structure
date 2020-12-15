#!/usr/bin/perl
die "perl $0 enriched_pixels cluster_of_pixels.bedpe virus.final.CON\n" if(@ARGV != 3);
my $enriched_pixels_bedpe=shift;
my $cluster_of_pixels_bedpe=shift;
my $phase1_final_constrain_file=shift;

#public parameters
my $ref_seq_fasta_file="WuHan_Hu_1.fasta";	#the file name of target sequence, in fasta format
my $whole_RNA_length=29903;		#the length of target sequence
my $RICseq_count_matrix_for_correlation="../1.Loop_by_HiCCUPS/virus.onlyPart.resolution5.format.matrix";          #matrix used for correlation analysis
my $ruler_used=5; 
my $output_sequence_id="NC_045512.2";

my $longest_loop_length=3000;		#limit to save running time

my $extended_slop_ratio=0.05;		#slightly slop sloops 
my $min_extended_len=5;			#minimal extended slop length
my $max_extended_len=25;		#maximal extended slop length

my $boundary_effect_ratio=0.05;		#do not force structure at boundary
my $min_boundary_effect_len=10;		#minimal boundary effect length
my $max_boundary_effect_len=20;		#maximal boundart effect length 

my $cutoff_for_inner_distance=600;	#classify as long loop or not

my $cutoff_to_closest=10;		#required to remove basepairs not supported by RIC-seq

my $maxExtended_solved=750;		#length of extended resolved region
#public parameters over

##test parameters
my $print_paired_pixels=0;
my $print_pvalue_eachStructure=1;
my $test_region_start=0;
my $test_region_end=300000;
##test parameters over

#enriched pixels
my %enriched_pixels_loci;
open(EPB,$enriched_pixels_bedpe) || die;
<EPB>;
while(my $line=<EPB>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$enriched_pixels_loci{$sub[6]}=join"\t",@sub[0..5];
}

#loop from enriched pixels 
my %raw_loops;
my %all_loops;
my %loop_len;
my %inner_distance;
my %loops_bedpe;
my @running_loops;
my @running_info;
open(CPB,$cluster_of_pixels_bedpe) || die;
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

	$inner_distance{$name}=$y_start-$x_end;
	$loops_bedpe{$name}=$x_start."\t".$x_end."\t".$y_start."\t".$y_end;

	my @four=($x_start,$x_end,$y_start,$y_end);
	@four=sort {$a<=>$b} @four;

	if($four[3] - $four[0] > $longest_loop_length){	#limited region
		next;
	}

	my $extended_slop=int($extended_slop_ratio*($four[3]-$four[0]));
	$extended_slop = $extended_slop < $min_extended_len ? $min_extended_len : $extended_slop;
	$extended_slop = $extended_slop > $max_extended_len ? $max_extended_len : $extended_slop;

	$raw_loops{$name}=$four[0]."\t".$four[3];
	$four[0]=$four[0] < $extended_slop ? 0 : $four[0]-$extended_slop;
	$four[3]=$four[3] > $whole_RNA_length-$extended_slop ? $whole_RNA_length : $four[3]+$extended_slop;
	$all_loops{$name}=$four[0]."\t".$four[3];
	$loop_len{$name}=$four[3]-$four[0];
	push (@running_loops,$name);
	push (@running_info,$line);
}

#finished loop reading

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

#read RIC-seq count matrix
my %RICseq_matrix;
open(MT,$RICseq_count_matrix_for_correlation) || die;
my $head=<MT>;
chomp $head;
my @arr_head=split/\s+/,$head;
while(my $line=<MT>){
	chomp $line;
	my @sub=split/\s+/,$line;

	#test####################################
	if($sub[0] < $test_region_start){	#
		next;				#
	}					#
	if($sub[0] > $test_region_end){		#
		next;				#
	}					#
	#test####################################
	
	foreach (1..$#sub){
		if($arr_head[$_] <= $sub[0]){
			next;
		}
		$RICseq_matrix{$sub[0]}{$arr_head[$_]}=$sub[$_];
	}
}

#finished RIC-seq count matrix

#dynamic prediction results
my %dynamic_constrains;
my %dynamic_basepairs;

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
                $dynamic_basepairs{$pair[0]}=$pair[1];  #1-based
                $dynamic_basepairs{$pair[1]}=$pair[0];  #1-based
        }
}

#my %input_dynamic_basepairs=%dynamic_basepairs;
#finished initialize dynamic solved
#dynamic prediction results

#running
my %resolved_bases;
open(USEDLP,">log.used_loops.bed") || die;
#open(RESOLVE,">log.resolved_bases.bed") || die;

my $run_cycle;
foreach my $id (@running_loops){
	my ($loop_start,$loop_end)=split/\s+/,$all_loops{$id};
	my $loop_details=shift @running_info;

	#for test use
	if($loop_start < $test_region_start or $loop_end > $test_region_end){
		next;
	}
	#for test use over

	$run_cycle++;
	print "-----------------------running cycle $run_cycle-----------------------\n";
	print $loop_details,"\n";
	my ($raw_loop_start,$raw_loop_end)=split/\s+/,$raw_loops{$id};
	my $raw_loop_length=$raw_loop_end-$raw_loop_start;
	my ($extended_loop_start,$extended_loop_end)=extended_to_farest_bp($loop_start,$loop_end);
	if($extended_loop_start eq "All" and $extended_loop_end eq "All"){
		print "$id\tregions is all resolved\n";
		next;   #the loop covered regions is all resolved
	}

	my ($extended_loop_start,$extended_loop_end)=shrink_to_farest_unsolved($extended_loop_start,$extended_loop_end);
	$extended_loop_start=$extended_loop_start-1;
        $extended_loop_end=$extended_loop_end+1;
        $extended_loop_start=$extended_loop_start < 0 ? 0 : $extended_loop_start;
        $extended_loop_end=$extended_loop_end > $whole_RNA_length ? $whole_RNA_length : $extended_loop_end;

	my $extended_loop_len=$extended_loop_end-$extended_loop_start;
	my $extended_loop_seq=substr($genome_seq,$extended_loop_start,$extended_loop_len);

	if($extended_loop_len <= 0){
		print "$id\tregion is too short, maybe due to bug, please check\n";
		exit;
	}

	my $tmp_unresolved_blocks=get_resolved_blocks_in_a_loop($extended_loop_start,$extended_loop_end);

	if($inner_distance{$id} <= $cutoff_for_inner_distance){
		if($unique_unresolved_blocks{$tmp_unresolved_blocks}){
			print $tmp_unresolved_blocks,"The unresolved region in this extended loops is identical to previous extended loop\nPrevious loops are:\n";
			print $unique_unresolved_blocks{$tmp_unresolved_blocks};
			next;
		}
		else{
			print $tmp_unresolved_blocks,"The unresolved region in this extended loops is novel\n";
			$unique_unresolved_blocks{$tmp_unresolved_blocks}.=$id."\n";
		}
	}

	my $boundary_effect_len=int($boundary_effect_ratio*$extended_loop_len);
	$boundary_effect_len = $boundary_effect_len < $min_boundary_effect_len ? $min_boundary_effect_len : $boundary_effect_len;
	$boundary_effect_len = $boundary_effect_len > $max_boundary_effect_len ? $max_boundary_effect_len : $boundary_effect_len;

	open(LOOPS,">tmp.$id.fa") || die;
	print LOOPS ">$id\n$extended_loop_seq\n";
	close LOOPS;
	my $fa_file="tmp.$id.fa";
	my $ct_file="tmp.$id.ct";

	#add constrains
	my $have_cons=collect_constrain($extended_loop_start,$extended_loop_end);

	print $id,"\t",$have_cons,"\n";
	print $extended_loop_seq,"\n";

	#run prediction
	if($have_cons){
		`/mnt/pub/work/caochch/project_My/project030_RIC_Covid19Virus/software/RNAstructure/exe/Fold $fa_file $ct_file --constraint tmp.CON --maxdistance 2500`;
	}
	else{
		`/mnt/pub/work/caochch/project_My/project030_RIC_Covid19Virus/software/RNAstructure/exe/Fold $fa_file $ct_file --maxdistance 2500`;
	}
	if(-e $ct_file){
	}
	else{
		print "No available structure\n";
		next;
	}
	my ($index,$energy,$pvalue,$ctinfo)=select_by_correlation($extended_loop_start,$extended_loop_end,$ct_file,$loop_start,$loop_end);
	print $index,"\t",$energy,"\t",$pvalue,"\tselected\n";
	#`rm -rf $fa_file $ct_file`;
	if($pvalue eq "NA"){
		next;
	}
	if($pvalue > 0.05){
		next;
	}

	#update all structure
	my @ct_lines=split/\n/,$ctinfo;
	shift @ct_lines;
	foreach my $l (@ct_lines){
		$l=~s/^\s+//;
		my @sub=split/\s+/,$l;
		if($#sub != 5){
			next;
		}
		
		if($sub[4]){
			my $real_loci_left=$extended_loop_start+$sub[0];
			my $real_loci_right=$extended_loop_start+$sub[4];

                        if($inner_distance{$id} <= $cutoff_for_inner_distance){
                                if($real_loci_left <= $raw_loop_start + $boundary_effect_len or $real_loci_left >= $raw_loop_end - $boundary_effect_len){
                                        next;
                                }
                                if($real_loci_right <= $raw_loop_start + $boundary_effect_len or $real_loci_right >= $raw_loop_end - $boundary_effect_len){
                                        next;
                                }
                        }
                        else{
				my @loops_arms_loci=split/\s+/,$loops_bedpe{$id};
				if($real_loci_left > $loops_arms_loci[1] or $real_loci_right < $loops_arms_loci[2]){
					next;
				}
                        }

			my @pair=($real_loci_left,$real_loci_right);
			@pair=sort {$a<=>$b} @pair;
			
			#evaluate the duplex possibility
			if(closest_enriched_pixels($pair[0],$pair[1]) <= $cutoff_to_closest){
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
				foreach ($pair[0]..$pair[1]){
					$resolved_bases{$_}=1;
				}
			}
			else{
				next;
			}
		}
	}

	#update resolved blocks
	#update_resolved_bases();
	#open(RESOLVE,">log.resolved_bases.cycle$run_cycle.bed") || die;
        #my $resolved_bases_blocks=print_resolved_region();
        #print RESOLVE $resolved_bases_blocks;
        #close RESOLVE;
	print USEDLP "$output_sequence_id\t$loop_start\t$loop_end\t$id\n";
	#update resolved regions and print over

	#test block#############################
	my $tmp_prefix="phase1b_cycle$run_cycle";
	visualize_by_IGV($tmp_prefix);
	visualize_by_ct($tmp_prefix);
}

write_constrains("phase1b_final");
visualize_by_dot("phase1b_final");
visualize_by_IGV("phase1b_final");
visualize_by_ct("phase1b_final");

open(FINALRE,">log.resolved_bases.phase1_final.bed") || die;
my $resolved_bases_blocks=print_resolved_region();
print FINALRE $resolved_bases_blocks;
close FINALRE;

sub get_resolved_blocks_in_a_loop{
	my $start=shift;
	my $end=shift;
	$start=$start+1;
	my $resolved_blocks;
	my $begin_to_record;
	my $before_loci;
	foreach my $i ($start..$end){
		if($resolved_bases{$i}){
			next;
		}
		else{
			if($begin_to_record){
				if($i == $before_loci+1){
					$before_loci=$i;
					next;
				}
				else{
					$resolved_blocks.=$before_loci."\n";
					$resolved_blocks.="$output_sequence_id\t$i\t";
					$before_loci=$i;
				}
			}
			else{
				$begin_to_record=1;
				$before_loci=$i;
				$resolved_blocks.="$output_sequence_id\t$i\t";
			}
		}
	}
	$resolved_blocks.=$end."\n";
	return $resolved_blocks;
}

sub shrink_to_farest_unsolved{
	my $start=shift;
	my $end=shift;
	$start=$start+1;
	my $left_most=$start;
	my $right_most=$end;
	print "From $left_most to $right_most before shrinking\n";
	while(1){
		if(!$resolved_bases{$left_most}){
			last;
		}
		$left_most=$left_most+1;
	}
	while(1){
		if(!$resolved_bases{$right_most}){
			last;
		}
		$right_most=$right_most-1;
	}
	print "From $left_most to $right_most after shrinking\n";
	$left_most=$left_most-1;        #to 0-based
	$right_most=$right_most;      #to 0-based
	return ($left_most,$right_most);
}

sub extended_to_farest_bp{
        my $start=shift; #0-based
	my $end=shift;
	$start=$start+1;
        my $left_most=$start;
        my $right_most=$end;

        print "From $left_most to $right_most before extending\n";
	my $all_resolved=1;
	foreach ($start..$end){
		if(!exists $resolved_bases{$_}){
			$all_resolved=0;
			last;
		}
	}
	if($all_resolved){
		print "From $left_most to $right_most is all resolved\n";
		return ("All","All");
	}
	
        foreach my $i ($start..$end){ #to 1-based to fit %dynamic_basepairs
                if(exists $dynamic_basepairs{$i}){
                        my $other_side=$dynamic_basepairs{$i};  #1-based
                        if($other_side < $left_most){
                                $left_most=$other_side;
                        }
                        if($other_side > $right_most){
                                $right_most=$other_side;
                        }
                }
        }

        print "From $left_most to $right_most after exnteding\n";

        $left_most=$left_most-1;        #to 0-based
        $right_most=$right_most;      #to 0-based
        return ($left_most,$right_most);
}



sub extended_loops_region{
	my $start=shift;#0-based, included
	my $end=shift;#0-based, not included
	$end=$end-1;
	my $all_resolved=1;
	foreach ($start..$end){
		if(!exists $resolved_bases{$_}){
			$all_resolved=0;
			last;
		}
	}
	if($all_resolved){
		return ("All","All");
	}

	my $start_need_to_extend;
	my $end_need_to_extend;
	if($resolved_bases{$start}){
		$start_need_to_extend=1;
	}
	if($resolved_bases{$end}){
		$end_need_to_extend=1;
	}
	foreach (1..$maxExtended_solved){
		if($resolved_bases{$start}){
			$start=$start-1;
		}
		else{
			last;
		}
	}
	foreach (1..$maxExtended_solved){
		if($resolved_bases{$end}){
			$end=$end+1;
		}
		else{
			last;
		}
	}
	if($start_need_to_extend){
		$start=$start+1;
	}
	if($end_need_to_extend){
		$end=$end-1;
	}
	$end=$end+1;
	return ($start,$end);
}

sub print_resolved_region{
	my @all_bases=sort {$a<=>$b} keys %resolved_bases;
	my $bedfiles;
	foreach my $i (0..$#all_bases){
		if($i == 0){
			my $zerobased=$all_bases[$i]-1;
			$bedfiles="$output_sequence_id\t$zerobased\t";
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
				$bedfiles.="$output_sequence_id\t$zerobased\t";
			}
		}
	}
	return $bedfiles;
}

sub closest_enriched_pixels{
	my $loci_x=shift;#smaller
	my $loci_y=shift;#bigger
	my $min_distance=30000;
	foreach my $p (keys %enriched_pixels_loci){
		my $dis_x;
		my $dis_y;
		my ($chr_x,$start_x,$end_x,$chr_y,$start_y,$end_y)=split/\s+/,$enriched_pixels_loci{$p};
		if(abs($start_x-$loci_x) > abs($end_x-$loci_x)){
			$dis_x=abs($end_x-$loci_x);
		}
		else{
			$dis_x=abs($start_x-$loci_x);
		}
		if(abs($start_y-$loci_y) > abs($end_y-$loci_y)){
			$dis_y=abs($end_y-$loci_y);
		}
		else{
			$dis_y=abs($start_y-$loci_y);
		}
		my $dis=$dis_x > $dis_y ? $dis_x : $dis_y;
		if($dis < $min_distance){
			$min_distance=$dis;
		}
	}
	return $min_distance;
}

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


sub visualize_by_IGV{
	#visualize by IGV
	my $prefix=shift;
	open(IGV,">virus.$prefix.igv.bed") || die;
	print IGV "track graphType=arc\n";
	foreach (keys %dynamic_constrains){
		my ($loci_left,$loci_right)=split/\s+/,$_;
		$loci_left=$loci_left-1;
		$loci_right=$loci_right-1;
		print IGV "$output_sequence_id\t$loci_left\t$loci_right\n";
	}
	close IGV;
}

sub visualize_by_ct{
	#visualize by ct
	my $prefix=shift;
	open(WCT,">virus.$prefix.whole.ct") || die;
	print WCT "$whole_RNA_length\tCovid19\n";
	my $offset=0;	#for the purpose of start at a given loci
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
	#print DOT substr($genome_seq,$test_region_start,$test_len),"\n";	#test
	#print DOT substr($string,$test_region_start,$test_len),"\n";		#test
}

sub collect_constrain{
	my $start=shift;
	my $end=shift;
	my %tmp_constrain_pairs;
	my %tmp_constrain_singlestrand;
	print $start,"\t",$end,"\tloop region\n";
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

	if(%tmp_constrain_pairs or %tmp_constrain_singlestrand){	
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
			my $relative_loci_left=$loci_left-$start;
			my $relative_right_loci=$loci_right-$start;
			print CON $relative_loci_left," ",$relative_right_loci,"\n";
			substr($testString,$relative_loci_left-1,1)="(";
			substr($testString,$relative_right_loci-1,1)=")";
			
		}
		print $testString,"\n";
		print CON "-1 -1\n";
		print CON "FMN:\n-1\nForbids:\n-1 -1\n";
		close CON;
		return 1;
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
		if($id_of_structure > 5){	#only top 5 structures were used
			last;
		}
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

		if(@paired < 3 or @other < 3){
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
		if($print_pvalue_eachStructure){
			print ">",$id_of_structure,"\t",$#paired+1,"\t",$paired_sum,"\t",$mean_paired,"\t",$median_paired,"\t",$#other+1,"\t",$other_sum,"\t",$mean_other,"\t",$median_other,"\t",$pvalue,"\ttestPvalue\n";
		}
	}
	$/="\n";

	if(%structure){
		#my @sorted_strcuture=sort {$structure{$a}{"pvalue"} <=> $structure{$b}{"pvalue"}} keys %structure;	#min pvalue
		my @sorted_strcuture=sort {$structure{$b}{"sum"} <=> $structure{$a}{"sum"}} keys %structure;		#max sum
		#return ($sorted_strcuture[0],$structure{$sorted_strcuture[0]}{"energy"},$structure{$sorted_strcuture[0]}{"pvalue"},$structure{$sorted_strcuture[0]}{"ct"});

		my $optimal_value=$structure{$sorted_strcuture[0]}{"sum"};

		my @candidate_structure;
		foreach (@sorted_strcuture){
			if($structure{$_}{"sum"} >= $optimal_value){
				push (@candidate_structure,$_);
			}
		}
		
		my @sorted_again_structure=sort {$structure{$a}{"energy"} <=> $structure{$b}{"energy"}} @candidate_structure;	#lowest energy
		#my @sorted_again_structure=(1);	#test
		return ($sorted_again_structure[0],$structure{$sorted_again_structure[0]}{"energy"},$structure{$sorted_again_structure[0]}{"pvalue"},$structure{$sorted_again_structure[0]}{"ct"});
	}
	else{
		return ("NA","NA","NA","NA");
	}

}

	


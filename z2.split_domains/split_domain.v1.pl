#!/usr/bin/perl
die "perl $0 in.matrix\n" if(@ARGV != 2);
my $RICseq_count_matrix=shift;
my $shortest_domain=2000;

#public parameters
my $whole_RNA_length=29903;              #the length of target sequence
my $ruler_used=25;			#resolution
my $score_cutoff=0.01*2974.39;		#747.678 is the factor between score and juicebox vcrt count
#public parameters over

#read RIC-seq count matrix
my %RICseq_matrix;
open(MT,$RICseq_count_matrix) || die;
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
		$sub[$_] = $sub[$_] =~ /nan/i ? 0 : $sub[$_];
                $RICseq_matrix{$sub[0]}{$arr_head[$_]}=$sub[$_];
        }
}

my $last_win=int($whole_RNA_length/$ruler_used);
my $last_win_loci=$last_win*$ruler_used;
my %all_boundary;
$all_boundary{0}=1;
$all_boundary{$last_win_loci}=1;

while(1){
	my $need_to_separate=check();
	#print ">round$_\t$need_to_separate\tneed\n"; 
	my @current_boundary=sort {$a<=>$b} keys %all_boundary;
	#print "@current_boundary\t\n";
	if($need_to_separate){
		my $next_boundary=find_optimal_boundary();
		$all_boundary{$next_boundary}=1;
	}
	else{
		last;
	}
}

print "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\n";
my @final_boundary=sort {$a<=>$b} keys %all_boundary;
foreach (0..$#final_boundary-1){
	my $left=$final_boundary[$_];
	my $right=$final_boundary[$_+1];
	print "hs1\t$left\t$right\ths1\t$left\t$right\tDomain$_\t255\t+\t+\n";
}

sub find_optimal_boundary{
	my $best_boundary;
	my $best_ratio;
	foreach my $c (1..$last_win-1){
		my $candidate_loci=$c*$ruler_used;
		my $near_to_current=away_from_current_boundary($candidate_loci);
		#print $candidate_loci,"\t",$near_to_current,"\taaa\n";
		if($near_to_current){
			next;
		}
		else{
			my @candidate_boundary=keys %all_boundary;
			push (@candidate_boundary,$candidate_loci);
			@candidate_boundary=sort {$a<=>$b} @candidate_boundary;
			my ($count_intra,$sum_intra,$count_intra_significant,$sum_intra_significant,$count_inter,$sum_inter,$count_inter_significant,$sum_inter_significant)=calculate_RICseq_ratio(@candidate_boundary);
			#print $count_intra,"\t",$sum_intra,"\t",$count_intra_significant,"\t",$sum_intra_significant,"\t";
			#print $count_inter,"\t",$sum_inter,"\t",$count_inter_significant,"\t",$sum_inter_significant,"\n";
			my $intra_density=$sum_intra/$count_intra;
			my $inter_density=$sum_inter/$count_inter;
			my $ratio=$intra_density/$inter_density;
			#print $candidate_loci,"\t",$ratio,"\n";
			if($ratio > $best_ratio){
				$best_boundary=$candidate_loci;
				$best_ratio=$ratio;
			}
		}
		
	}
	return $best_boundary;
}

sub calculate_RICseq_ratio{
	my @tmp_boundary=@_;
	my $intra_count;
	my $intra_sum;
	my $intra_significant_count;
	my $intra_significant_sum;
	my $inter_count;
	my $inter_sum;
	my $inter_significant_count;
	my $inter_significant_sum;
	
	foreach my $i (0..$last_win-1){
		my $win_i=$i*$ruler_used;
		foreach my $j ($i..$last_win-1){
			my $win_j=$j*$ruler_used;
			my $is_intra=is_same_domain($win_i,$win_j,@tmp_boundary);
			if($is_intra eq "intra"){
				$intra_count++;
				$intra_sum+=$RICseq_matrix{$win_i}{$win_j};
				if($RICseq_matrix{$win_i}{$win_j} >= $score_cutoff){
					$intra_significant_count++;
					$intra_significant_sum+=$RICseq_matrix{$win_i}{$win_j};
				}
			}
			elsif($is_intra eq "inter"){
				$inter_count++;
				$inter_sum+=$RICseq_matrix{$win_i}{$win_j};
				if($RICseq_matrix{$win_i}{$win_j} >= $score_cutoff){
					$inter_significant_count++;
					$inter_significant_sum+=$RICseq_matrix{$win_i}{$win_j};
				}
			}
			else{
				die;
			}
		}
	}
	return ($intra_count,$intra_sum,$intra_significant_count,$intra_significant_sum,$inter_count,$inter_sum,$inter_significant_count,$inter_significant_sum);

}

sub is_same_domain{
	my $loci_i=shift;
	my $loci_j=shift;
	my @tmp_boundary=@_;
	my $domain_i;
	foreach my $index (0..$#tmp_boundary-1){
		if($loci_i >= $tmp_boundary[$index] and $loci_i < $tmp_boundary[$index+1]){
			$domain_i=$index;
			if($loci_j >= $tmp_boundary[$index] and $loci_j < $tmp_boundary[$index+1]){
				return "intra";
			}
		}
	}
	if(defined($domain_i)){
		return "inter";
	}
	else{
		warn "did not find domain index for $loci_i\n";
		die;
	}
}
	

sub away_from_current_boundary{
	my $candidate=shift;
	my @boundary=sort {$a<=>$b} keys %all_boundary;
	foreach (@boundary){
		if(abs($candidate-$_) < 1*$shortest_domain){
			return 1;
		}
	}
	return 0;
}

sub check{
	my @boundary=sort {$a<=>$b} keys %all_boundary;
	foreach my $i (0..$#boundary-1){
		if($boundary[$i+1] - $boundary[$i] > 2*$shortest_domain){
			return 1;
		}
	}
	return 0;
}




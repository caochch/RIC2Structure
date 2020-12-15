#!/usr/bin/perl
die "perl $0 all_pixels.list loop.list resolution scoreCutoff\n" if(@ARGV != 4);
my $all_pixels_list=shift;
my $all_loop_bedpe=shift;
my $resolution=shift;
my $score_cutoff=shift;#0.01


my %pixels_score;
open(APL,$all_pixels_list) || die;
while(my $line=<APL>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$pixels_score{$sub[0]}{$sub[1]}=$sub[4];
}

my %loop_loci_info;
my %loops_score_info;
open(ALB,$all_loop_bedpe) || die;
my $head=<ALB>;
print $head;
while(my $line=<ALB>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $loop_id=$sub[6];
	$loop_loci_info{$loop_id}=$line;
	my $left_start_index=int($sub[1]/$resolution);
	my $left_end_index=int($sub[2]/$resolution)-1;
	my $right_start_index=int($sub[4]/$resolution);
	my $right_end_index=int($sub[5]/$resolution)-1;
	my $pixel_count;
	my $pixel_significant_count;
	my $pixel_sum;
	my $pixel_significant_sum;
	foreach my $i ($left_start_index..$left_end_index){
		foreach my $j ($right_start_index..$right_end_index){
			my $win_i=$resolution*$i;
			my $win_j=$resolution*$j;
			if($j < $i+4){
				next;
			}
			$pixel_count++;
			$pixel_sum+=$pixels_score{$win_i}{$win_j};
			if($pixels_score{$win_i}{$win_j} >= $score_cutoff){
				$pixel_significant_count++;
				$pixel_significant_sum+=$pixels_score{$win_i}{$win_j};
			}
		}
	}

	#if($loop_id eq "ManualLoops_21"){
	#	print $pixel_count,"\t",$pixel_sum,"\t",$pixel_significant_count,"\t",$pixel_significant_sum,"\n";
	#}

	$loops_score_info{$loop_id}{"Count"}=$pixel_count;
	$loops_score_info{$loop_id}{"Sum"}=$pixel_sum;
	$loops_score_info{$loop_id}{"sigCount"}=$pixel_significant_count;
	$loops_score_info{$loop_id}{"sigSum"}=$pixel_significant_sum;
}

my %already_sorted;
my @running_loops;
my @running_info;

my @test_loops=("ManualLoops_215");

foreach my $id (sort {$loops_score_info{$b}{"sigCount"} <=> $loops_score_info{$a}{"sigCount"}} keys %loops_score_info){
#foreach my $id (@test_loops){
	if($already_sorted{$id}){
		next;
	}
	my @sorted_sons=sort_for_a_loop($id);
	foreach my $s (@sorted_sons){
		if($already_sorted{$s}){
			next;
		}
		else{
			push (@running_loops,$s);
			push (@running_info,$s."\tisSonOf\t".$id);
			$already_sorted{$s}=1;
		}
	}
}

foreach (@running_loops){
	#print $loops_score_info{$_}{"sigCount"},"\tsigCount\n";
	print $loop_loci_info{$_},"\t";
	print $loops_score_info{$_}{"sigCount"},"\t";
	my $info=shift @running_info;
	print $info,"\n";
}

sub sort_for_a_loop{
	my $id=shift;
	my @son_loops=find_son_loops($id);
	my @out;
	if(!@son_loops){
		@out=($id);
		return @out;
	}
	else{
		foreach my $s (sort {$loops_score_info{$b}{"sigCount"} <=> $loops_score_info{$a}{"sigCount"}} @son_loops){
			if($s eq $id){
				next;
			}
			else{
				push (@out,$s);
			}
		}
		push (@out,$id);
		return @out;
	}
}

sub find_son_loops{
	my $id=shift;
	my $grandFather_loci_info=$loop_loci_info{$id};
	my @current_group;

	foreach my $s (sort {$loops_score_info{$b}{"sigCount"} <=> $loops_score_info{$a}{"sigCount"}} keys %loop_loci_info){
		if($already_sorted{$s}){
			next;
		}
		my $son_loci_info=$loop_loci_info{$s};
		if(is_son($grandFather_loci_info,$son_loci_info) ne "son"){
			next;
		}
		my $is_suit=1;
		foreach my $f_id (@current_group){
			my $father_loci_info=$loop_loci_info{$f_id};
			my $index=is_conflict($father_loci_info,$son_loci_info);
			#print $s,"\t",$f_id,"\t",$index,"\n";
			if($index > 0.5){
				$is_suit=0;
				last;
			}
		}
		if($is_suit){
			push (@current_group,$s);
		}
		
	}
	return @current_group;	
}

sub is_conflict{
	my $f_loci=shift;
	my $s_loci=shift;
	my @f_loci_info=split/\s+/,$f_loci;
	my @s_loci_info=split/\s+/,$s_loci;
	my @f_four=($f_loci_info[1],$f_loci_info[2],$f_loci_info[4],$f_loci_info[5]);
	my @s_four=($s_loci_info[1],$s_loci_info[2],$s_loci_info[4],$s_loci_info[5]);

	#simple relationships
	if($s_four[0] >= $f_four[3] or  $s_four[3] <= $f_four[0]){	#un-related
		return 0;
	}
	if($s_four[1] <= $f_four[0] and $s_four[2] >= $f_four[3]){	#as a father
		return 0;
	}
	if($s_four[0] >= $f_four[1] and $s_four[3] <= $f_four[2]){	#as a son
		return 0;
	}
	
	#complexed
	my $all_pairs;
	my $conflict_pairs;
	foreach my $i ($s_four[0]..$s_four[1]-1){
		foreach my $j ($s_four[2]..$s_four[3]-1){
			if(pair_conflict($i,$j,@f_four)){
				$conflict_pairs++;
			}
			$all_pairs++;
		}
	}
	
	my $conflict_index=$conflict_pairs/$all_pairs;
	return $conflict_index;

}

sub pair_conflict{
	my $left_site=shift;
	my $right_site=shift;
	my @four=@_;
	if(($left_site >= $four[0] and $left_site <= $four[1]) or ($left_site >= $four[2] and $left_site <= $four[3]) or ($right_site >= $four[0] and $right_site <= $four[1]) or ($right_site >= $four[2] and $right_site <= $four[3])){
		if(($left_site >= $four[0] and $left_site <= $four[1]) and ($right_site >= $four[2] and $right_site <= $four[3])){
			return 0;	#not conflict
		}
		else{
			return 1;	#conflict
		}
	}
	else{
		return 0;#not conflict
	}
}


sub is_son{
	my $f_loci=shift;
	my $s_loci=shift;
	my @f_loci_info=split/\s+/,$f_loci;
	my @s_loci_info=split/\s+/,$s_loci;
	
	my @f_four=($f_loci_info[1],$f_loci_info[2],$f_loci_info[4],$f_loci_info[5]);
	my @s_four=($s_loci_info[1],$s_loci_info[2],$s_loci_info[4],$s_loci_info[5]);

	if($s_four[0] >= $f_four[0] and $s_four[3] <= $f_four[3]){
		if($s_four[0] >= $f_four[1] and $s_four[3] <= $f_four[2]){
			return "son";
		}
	}
	return "notSon";
}



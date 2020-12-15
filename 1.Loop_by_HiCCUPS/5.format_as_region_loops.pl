#!/usr/bin/perl
die "perl $0 cluster_of_pixels.bedpe\n" if(@ARGV != 1);
my $cluster_of_pixels_bedpe=shift;

open(CPB,$cluster_of_pixels_bedpe) || die;
while(my $line=<CPB>){
	if($line=~/^#/){
		print $line;
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
	$sub[1]=$four[0];
	$sub[2]=$four[3];
	$sub[4]=$four[0];
	$sub[5]=$four[3];

	$new=join"\t",@sub;
	print $new,"\n";
	
}


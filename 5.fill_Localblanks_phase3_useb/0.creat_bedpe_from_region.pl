#!/usr/bin/perl
die "perl $0 in.region.bed\n" if(@ARGV != 1);
my $region_bed=shift;

my $id;
open(RB,$region_bed) || die;
while(my $line=<RB>){
	chomp $line;
	$id++;
	my @sub=split/\s+/,$line;
	print "hs1\t$sub[1]\t$sub[2]\t";
	print "hs1\t$sub[1]\t$sub[2]\t";
	print "Solved_$id\t255\t+\t+\n";
}



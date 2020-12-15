#!/usr/bin/perl
die "perl $0 merged_loops.bedpe\n" if(@ARGV != 2);
my $hiccups_bedpe=shift;
my $resolution=shift;

my $window=5000;

my $factor=$window/$resolution;

my $id;
my $pixel_id;

print "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\n";

open(HB,$hiccups_bedpe) || die;
while(my $line=<HB>){
	if($line=~/^#/){
		next;
	}
	chomp $line;
	my @sub=split/\s+/,$line;
	$sub[0]="hs1";
	$sub[1]=int($sub[1]/$factor);
	$sub[2]=int($sub[2]/$factor);
	$sub[3]="hs1";
	$sub[4]=int($sub[4]/$factor);
	$sub[5]=int($sub[5]/$factor);
	my $new=join"\t",@sub[0..5];
	$pixel_id++;
	print $new,"\tEnrichPixel_",$pixel_id,"\t255\t+\t+\n";
}

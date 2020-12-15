#!/usr/bin/perl
die "perl $0 out1.all_pixels.score.list virus.resolution5.vcrt.matrix resolution\n" if(@ARGV != 3);
my $mine_matrix=shift;
my $juicebox_matrix=shift;
my $resolution=shift;

my $window=5000;

my $factor=$window/$resolution;

my %juicebox;
open(JM,$juicebox_matrix) || die;
while(my $line=<JM>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $win_a=$sub[0]/$factor;
	my $win_b=$sub[1]/$factor;
	$juicebox{$win_a}{$win_b}=$sub[2];
	$juicebox{$win_b}{$win_a}=$sub[2];
}

open(MM,$mine_matrix) || die;
while(my $line=<MM>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($juicebox{$sub[0]}{$sub[1]}){
		print $line,"\t",$juicebox{$sub[0]}{$sub[1]},"\n";
	}
	else{
		print $line,"\t","NA\n";
	}
}

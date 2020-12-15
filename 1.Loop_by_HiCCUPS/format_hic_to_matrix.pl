#!/usr/bin/perl
die "perl inMatrix Chr ruler resolution\n" if(@ARGV != 3);
my $in_matrix=shift;
my $chr_size=shift;       #28s
my $ruler=shift;
my $resolution=$ruler;

my $window=5000;

my $factor=$window/$resolution;


my %matrix;
open(IN,$in_matrix) || die;
while(my $line=<IN>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$sub[0]=$sub[0]/$factor;
	$sub[1]=$sub[1]/$factor;
	$sub[2] = $sub[2]=~/NaN/ ? 0 : $sub[2];
	$matrix{$sub[0]}{$sub[1]}=$sub[2];
	$matrix{$sub[1]}{$sub[0]}=$sub[2];
}

my @all;
foreach (0..$chr_size){
	push (@all,$ruler*$_);
}
	

print "loci\t";

foreach (@all){
	print $_,"\t";
}
print "\n";

foreach my $row (@all){
	print $row,"\t";
	foreach my $col (@all){
		print $matrix{$row}{$col}+0,"\t";
	}
	print "\n";
}

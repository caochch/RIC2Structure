#!/usr/bin/perl
die "perl $0 in.fa in.igv.bed\n" if(@ARGV != 2);
my $in_fa=shift;
my $igv_bed=shift;

my $refseq;
open(IF,$in_fa) || die;
while(my $line=<IF>){
        chomp $line;
        if($line=~/^>/){
                next;
        }
        $refseq.=$line;
}

my %allowed_basepair=("AT" => 1,"TA" => 1,"CG" => 1,"GC" => 1,"GT"=>1,"TG"=>1);

open(IB,$igv_bed) || die;
<IB>;
while(my $line=<IB>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $left=$sub[1]+1;
	my $right=$sub[2]+1;
        if($right < $left){
       		die;
        }
	my $this_pair=substr($refseq,$left-1,1).substr($refseq,$right-1,1);
	my $outside_pair=substr($refseq,$left-2,1).substr($refseq,$right,1);
	my $inside_pair=substr($refseq,$left,1).substr($refseq,$right-2,1);
	print $allowed_basepair{$this_pair}+0,"\t",$allowed_basepair{$inside_pair}+0,"\t",$allowed_basepair{$outside_pair}+0,"\n";
}

	

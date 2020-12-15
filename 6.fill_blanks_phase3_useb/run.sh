ln -s ../5.fill_Localblanks_phase3_useb/virus.phase3_final.CON ./
rm -rf tmp* *igv.bed log* debug.txt *ct

perl 1.fill_blanks.useRNAstructure.v3.pl virus.phase3_final.CON ../z2.split_domains/out1.domains.bedpe > debug.txt

read.table(file="tmp.paired.RICseq.signal.list") -> intact
read.table(file="tmp.other.RICseq.signal.list") -> random
t.test(intact,random,alternative = "greater") -> test
cat(test$p.value,"\t")
cat(mean(intact$V1),"\t")
cat(median(intact$V1),"\t")
cat(mean(random$V1),"\t")
cat(median(random$V1),"\t")

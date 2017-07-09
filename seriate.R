
list.of.packages <- c("seriation")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library("seriation")

f <- file("stdin")
open(f)
a_read = read.csv(f)
close(f)

nums = a_read

d <- dist(nums)
# g_methods = c("QAP_BAR","OLO_ward","OLO_average","OLO","HC_ward","HC_average","GW_ward","SPIN_NH","TSP")
o = seriate(d, "OLO")
r = get_order(o)
cat("this is the start of the data\n")
write.table(r,"",row.names=FALSE,col.names=FALSE,sep=',')

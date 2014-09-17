setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/HutteriteTTP")
results = read.table('TAP2_TTP+eqtlwgeno+imp.txt', as.is=T, header=T)
map = read.table('HuttInfoKey.INFO', as.is=T, header=T)
results_RSID = merge(results, map, by.x = 'stopBP', by.y = 'POS', all.x=T, all.y=F)
write.table(results_RSID, 'TAP2_TTP=eqtlwgeno+imp_rsID.txt', sep ='\t', row.names = F)

results = read.table('HLAF_TTP+eqtlwgeno+imp.txt', as.is=T, header=T)
results_RSID = merge(results, map, by.x = 'stopBP', by.y = 'POS', all.x=T, all.y=F)
write.table(results_RSID, 'HLAF_TTP=eqtlwgeno+imp_rsID.txt', sep ='\t', row.names = F)
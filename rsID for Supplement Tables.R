setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/HutteriteTTP")
results = read.table('TAP2_TTP+eqtlwgeno+imp.txt', as.is=T, header=T)
map = read.table('HuttInfoKey.INFO', as.is=T, header=T)
results_RSID = merge(results, map, by.x = 'stopBP', by.y = 'POS', all.x=T, all.y=F)
write.table(results_RSID, 'TAP2_TTP=eqtlwgeno+imp_rsID.txt', sep ='\t', row.names = F)

results = read.table('HLAF_TTP+eqtlwgeno+imp.txt', as.is=T, header=T)
results_RSID = merge(results, map, by.x = 'stopBP', by.y = 'POS', all.x=T, all.y=F)
write.table(results_RSID, 'HLAF_TTP=eqtlwgeno+imp_rsID.txt', sep ='\t', row.names = F)

#Updated with complete rsID list
setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/HutteriteTTP")
results = read.table('TAP2_TTP+eqtlwgeno+imp.txt', as.is=T, header=T)
map = read.table('HuttInfoKey_Updated.txt', as.is=T, header=T)
key = read.table('HuttKey_onlyupdated.txt', as.is=T, header=T)
results_RSID = merge(results, map, by.x = 'stopBP', by.y = 'POS', all.x=T, all.y=F)
results_updated = merge(results_RSID, key, by.x= 'MarkerID', by.y = 'Marker', all.x=T, all.y=F)
write.table(results_updated, 'TAP2_TTP=eqtlwgeno+imp_rsID_updated.txt', sep ='\t', row.names = F)

results = read.table('HLAF_TTP+eqtlwgeno+imp.txt', as.is=T, header=T)
results_RSID = merge(results, map, by.x = 'stopBP', by.y = 'POS', all.x=T, all.y=F)
results_updated = merge(results_RSID, key, by.x= 'MarkerID', by.y = 'Marker', all.x=T, all.y=F)
write.table(results_updated, 'HLAF_TTP=eqtlwgeno+imp_rsID_updated.txt', sep ='\t', row.names = F)
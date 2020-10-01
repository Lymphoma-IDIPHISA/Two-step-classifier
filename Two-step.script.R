# Two-step genetic classification example

MD = read.table("/path/to/example.txt", sep = "\t", header = T, row.names = T)
gene = colnames(MD)
subtypes = c("MCD.MYD88", "BN2.NOTCH2", "N1", "EZB.BCL2", "ST2.SGK1")
MD$MCD.MYD88 = c(0)
MD$BN2.NOTCH2 = c(0)
MD$N1 = c(0)
MD$EZB.BCL2 = c(0)
MD$ST2.SGK1 = c(0)
MD$Subtype = NA

for (l in 1:nrow(MD)) {
  genesmu = as.vector(MD[l,])
  p = which(genesmu == 1)
  genesmuta = gene[p]
  
  if (MD$MYD88[l] == 1) {MD$MCD.MYD88[l] = MD$MCD.MYD88[l] + 1}
  if (MD$PIM1[l] == 1) {MD$MCD.MYD88[l] = MD$MCD.MYD88[l] + 1}
  if (MD$CD79B[l] == 1) {MD$MCD.MYD88[l] = MD$MCD.MYD88[l] + 1}
  
  if (MD$EZH2[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  if (MD$BCL2[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  if (MD$CREBBP[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  if (is.na(MD$BCL2.fusion[l]) == F && MD$BCL2.fusion[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  
  if (MD$NOTCH2[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}
  if (MD$BCL10[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}
  if (MD$TNFAIP3[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}
  if (is.na(MD$BCL6.fusion[l]) == F && MD$BCL6.fusion[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}
  
  if (MD$TET2[l] == 1) {MD$ST2.SGK1[l] = MD$ST2.SGK1[l] + 1}
  if (MD$SGK1[l] == 1) {MD$ST2.SGK1[l] = MD$ST2.SGK1[l] + 1}
  if (MD$SOCS1[l] == 1) {MD$ST2.SGK1[l] = MD$ST2.SGK1[l] + 1}
  
  p = MD[l, (ncol(MD)-5):(ncol(MD)-1)]
  z = which(p == 4)
  if (length(z) == 1 && is.na(MD$Subtype[l]) == T && length(z) != 2) {MD$Subtype[l] = subtypes[z]}
  
  y = which(p == 3)
  if (length(y) == 1 && is.na(MD$Subtype[l]) == T && length(y) != 2 && length(z) != 2) {MD$Subtype[l] = subtypes[y]}
  
  x = which(p == 2)
  if (length(x) == 1 && is.na(MD$Subtype[l]) == T && length(x) != 2 && length(y) != 2) {MD$Subtype[l] = subtypes[x]}
  
  b = which(p == 1)
  if (length(b) == 1 && is.na(MD$Subtype[l]) == T && length(b) != 2 && length(x) != 2) {MD$Subtype[l] = subtypes[b]}
}


for (l in 1:nrow(MD)) {
  genesmu = as.vector(MD[l,])
  
  if (MD$PRDM1[l] == 1) {MD$MCD.MYD88[l] = MD$MCD.MYD88[l] + 1}
  if (MD$BTG1[l] == 1) {MD$MCD.MYD88[l] = MD$MCD.MYD88[l] + 1}
  if (MD$PIM2[l] == 1) {MD$MCD.MYD88[l] = MD$MCD.MYD88[l] + 1}
  if (MD$CD58[l] == 1) {MD$MCD.MYD88[l] = MD$MCD.MYD88[l] + 1}
  
  if (MD$TNFRSF14[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  if (MD$KMT2D[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  if (MD$IRF8[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  if (MD$EP300[l] == 1) {MD$EZB.BCL2[l] = MD$EZB.BCL2[l] + 1}
  
  if (MD$CD70[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}
  if (MD$DTX1[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}
  if (MD$UBE2A[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}
  if (MD$CCND3[l] == 1) {MD$BN2.NOTCH2[l] = MD$BN2.NOTCH2[l] + 1}

  if (MD$STAT3[l] == 1) {MD$ST2.SGK1[l] = MD$ST2.SGK1[l] + 1}
  
  if (MD$NOTCH1[l] == 1) {MD$Subtype[l] = "N1"}
  
  p = MD[l, (ncol(MD)-5):(ncol(MD)-1)]
  b = which(p == 5)
  if (length(b) == 1 && is.na(MD$Subtype[l]) == T && length(b) != 2) {MD$Subtype[l] = subtypes[b]}
  
  z = which(p == 4)
  if (length(z) == 1 && is.na(MD$Subtype[l]) == T && length(z) != 2 && length(b) != 2) {MD$Subtype[l] = subtypes[z]}
  
  x = which(p == 3)
  if (length(x) == 1 && is.na(MD$Subtype[l]) == T && length(x) != 2 && length(z) != 2) {MD$Subtype[l] = subtypes[x]}
  
  y = which(p == 2)
  if (length(y) == 1 && is.na(MD$Subtype[l]) == T && length(y) != 2 && length(x) != 2) {MD$Subtype[l] = subtypes[y]}
}

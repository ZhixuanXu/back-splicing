rm(list=ls());gc();
setwd('working passway')
datDir <- 'data pasway'
sps <- c('hu','ma','mo')
tis <- c('Brain','Heart','Kidney','Liver','Lung','Prostate','SkeletalMuscle','SpinalCord','Spleen','Stomach','Testis')
parameter <- read.table(file='parameter_fixApt.xls',sep="\t",strip.white = TRUE, stringsAsFactors=FALSE,header=TRUE,na.strings="NA",fill=TRUE)
infType <- c('bsIncDnrNum')
args<-commandArgs(T)
for (i in c(as.numeric(args[1]))){
	dat <- read.table(file=paste(datDir,sps[i],'TpmRdsSjs_byJunc_s.xls',sep=''),sep="\t",strip.white = TRUE, stringsAsFactors=FALSE,header=TRUE,na.strings="NA",fill=TRUE)
	dat <- dat[!is.na(dat[3]) & dat[2] == 'protein_coding',]
	dat$bsLen <- abs(dat[,6]-dat[,7])+1
	#for (j in 1:11){
	for (j in c(as.numeric(args[2]))){
		tisPar <- parameter[parameter[,1]==sps[i] & parameter[,2] == tis[j],]
		tisDat <- dat[dat[,(j-1)*2+10]>0,c(1:8,(j-1)*2+9,(j-1)*2+10,31)]
		bsGens <- unique(tisDat[tisDat[,4]=='bsj',1])
		nr <- 1
		newBsInfs <- as.data.frame(matrix(ncol=20,nrow=0))
		colnames(newBsInfs) <- c(colnames(tisDat)[1:8],'acceptorCrd',
								'aptBsRds','aptRds','aptItnLen','bsLen','bsIncDnrNum','dnrItnLen','dnrRds', 
								'bsLenIncDnrNum_P', 'simulated_reads',
								'adj_bsLenIncDnrNum_P', 'adj_simulated_reads') 
		newBsInfs2 <- newBsInfs
		for (z in 1:length(bsGens)){
			genDat <- tisDat[tisDat[,1] == bsGens[z],]
			genBsDat <- genDat[genDat[,4]=='bsj',]
			genLsDat <- genDat[genDat[,4]=='lsj',]
			
			std <- genDat[1,8]
			if (std == '+'){
				bsApts <- sort(unique(genDat[genDat[,4]=='bsj',6]))
				dnrs <- sort(unique(c(unique(genDat[genDat[,4]=='bsj',7]),unique(genDat[genDat[,4]=='lsj',6]))))
				
				for (k in 1:length(bsApts)){
					aptBsRds <- sum(genBsDat[genBsDat[,6]==bsApts[k],10])
					aptRds <- c()
					aptItnLen <- 0
					bsLen <- c()
					bsIncDnrNum <- c()
					dnrItnLen <- 0
					dnrRds <- c()
					aptRds <- sum(genDat[genDat[,6]==bsApts[k] | genDat[,7]==bsApts[k],10])
					genLsDat_apt <- genLsDat[genLsDat[,7]==bsApts[k],]
					if (length(genLsDat_apt[,1])>0){
						aptItnLen1 <- c()
						for (n in 1:length(genLsDat_apt[,1])){
							aptItnLen1 <- c(aptItnLen1,rep(genLsDat_apt[n,11],genLsDat_apt[n,10]))
						}
						aptItnLen <- mean(aptItnLen1)
					}
					dnrs2 <- dnrs[dnrs>bsApts[k]]
					if (length(dnrs2)>0){
						for (m in 1:length(dnrs2)){
								bsLen <- abs(dnrs2[m] - bsApts[k])+1
								bsIncDnrNum <- length(dnrs2[dnrs2<dnrs2[m] & dnrs2>bsApts[k]])
								genLsDat_dnr <- genLsDat[genLsDat[,6]==dnrs2[m],]
								if (length(genLsDat_dnr[,1])>0){
									dnrItnLen1 <- c()
									for (n in 1:length(genLsDat_dnr[,1])){
										dnrItnLen1 <- c(dnrItnLen1,rep(genLsDat_dnr[n,11],genLsDat_dnr[n,10]))
									}
									dnrItnLen <- mean(dnrItnLen1)
								}
								dnrRds <- sum(genDat[genDat[,7]==dnrs2[m] | genDat[,6]==dnrs2[m],10]) 
							newBsInfs[nr,1:2] <- as.character(genDat[1,1:2,drop=T])
							newBsInfs[nr,4] <- 'bsj'
							newBsInfs[nr,5] <- genDat[1,5]
							newBsInfs[nr,8] <- genDat[1,8]
							newBsInfs[nr,6:7] <- c(bsApts[k],dnrs2[m])
							newBsInfs[nr,3] <- paste(genDat[1,5],bsApts[k],dnrs2[m],genDat[1,8],sep=',')
							newBsInfs[nr,9:16] <- c(bsApts[k],aptBsRds,aptRds,aptItnLen,bsLen,bsIncDnrNum,dnrItnLen,dnrRds)
							gPar <- tisPar[tisPar[,3] %in% 'bsIncDnrNum',]
							if (length(gPar[gPar[,5] %in% newBsInfs[nr,14],7])>0){
							newBsInfs[nr,17] <- gPar[gPar[,5] %in% newBsInfs[nr,14],7]
							}else{
								newBsInfs[nr,17] <- 0
							}
							nr <- nr + 1
						}
					}
				}
			}else if (std == '-'){
				bsApts <- sort(unique(genDat[genDat[,4]=='bsj',7]))
				dnrs <- sort(unique(c(unique(genDat[genDat[,4]=='bsj',6]),unique(genDat[genDat[,4]=='lsj',7]))))
				for (k in 1:length(bsApts)){
					aptBsRds <- sum(genBsDat[genBsDat[,7]==bsApts[k],10])
					aptRds <- c()
					aptItnLen <- 0
					bsLen <- c()
					bsIncDnrNum <- c()
					dnrItnLen <- 0
					dnrRds <- c()
					aptRds <- sum(genDat[genDat[,7]==bsApts[k] | genDat[,6]==bsApts[k],10])
					genLsDat_apt <- genLsDat[genLsDat[,6]==bsApts[k],]
					if (length(genLsDat_apt[,1])>0){
						aptItnLen1 <- c()
						for (n in 1:length(genLsDat_apt[,1])){
							aptItnLen1 <- c(aptItnLen1,rep(genLsDat_apt[n,11],genLsDat_apt[n,10]))
						}
						aptItnLen <- mean(aptItnLen1)
					}
					dnrs2 <- dnrs[dnrs<bsApts[k]]
					if (length(dnrs2)>0){
						for (m in 1:length(dnrs2)){
								bsLen <- abs(dnrs2[m] - bsApts[k])+1
								bsIncDnrNum <- length(dnrs2[dnrs2>dnrs2[m] & dnrs2<bsApts[k]])
								genLsDat_dnr <- genLsDat[genLsDat[,7]==dnrs2[m],]
								if (length(genLsDat_dnr[,1])>0){
									dnrItnLen1 <- c()
									for (n in 1:length(genLsDat_dnr[,1])){
										dnrItnLen1 <- c(dnrItnLen1,rep(genLsDat_dnr[n,11],genLsDat_dnr[n,10]))
									}
									dnrItnLen <- mean(dnrItnLen1)
								}
								dnrRds <- sum(genDat[genDat[,6]==dnrs2[m] | genDat[,7]==dnrs2[m],10])
							newBsInfs[nr,1:2] <- as.character(genDat[1,1:2,drop=T])
							newBsInfs[nr,4] <- 'bsj'
							newBsInfs[nr,5] <- genDat[1,5]
							newBsInfs[nr,8] <- genDat[1,8]
							newBsInfs[nr,6:7] <- c(dnrs2[m],bsApts[k])
							newBsInfs[nr,3] <- paste(genDat[1,5],dnrs2[m],bsApts[k],genDat[1,8],sep=',')
							newBsInfs[nr,9:16] <- c(bsApts[k],aptBsRds,aptRds,aptItnLen,bsLen,bsIncDnrNum,dnrItnLen,dnrRds)
							gPar <- tisPar[tisPar[,3] %in% 'bsIncDnrNum',]
							if (length(gPar[gPar[,5] %in% newBsInfs[nr,14],7])>0){ 
							newBsInfs[nr,17] <- gPar[gPar[,5] %in% newBsInfs[nr,14],7]
							}else{
								newBsInfs[nr,17] <- 0
							}
							nr <- nr + 1
						}
					}
				}
			}
		}
		for (z1 in 1:length(bsGens)){  #length(bsGens)
			simBsGen <- newBsInfs[newBsInfs[,1] == bsGens[z1],]
			simApts <- unique(simBsGen[,9])
			for (z2 in 1:length(simApts)){
				simAptDat <- simBsGen[simBsGen[,9] == simApts[z2],]
				simAptDat[,19] <- simAptDat[,17]/sum(simAptDat[,17])
				aptBsRds <- simAptDat[1,10]
				rdNum <- runif(aptBsRds,0,1)
				
				ori_p <- 0; adj_p <- 0; 
				for (rn in 1:length(simAptDat[,1])){
					ori_p <- ori_p + simAptDat[rn,17]
					adj_p <- adj_p + simAptDat[rn,19]
					if (rn == 1){
						simAptDat[rn,18] <- length(rdNum[rdNum<=ori_p])
						simAptDat[rn,20] <- length(rdNum[rdNum<=adj_p])
					}else{
						simAptDat[rn,18] <- length(rdNum[(rdNum<=ori_p) & (rdNum > (ori_p - simAptDat[rn,17]))])
						simAptDat[rn,20] <- length(rdNum[(rdNum<=adj_p) & (rdNum > (adj_p - simAptDat[rn,19]))])
					}
				}
				newBsInfs2 <- rbind(newBsInfs2,simAptDat)
			}
		}
		write.table(newBsInfs2,file=paste(sps[i],tis[j],'_simBS_fixApt.xls',sep=''),row.names=F,sep="\t",quote=F)
	}
}








library(TitanCNA)

args = commandArgs(trailingOnly=TRUE)
tumor_ref_alt = args[1]
tumWig = args[2]
normWig = args[3]
gcWig = args[4]
mapWig = args[5]
nClusters = args[6]
outFile = args[7]

data <- loadAlleleCounts(tumor_ref_alt, genomeStyle = "NCBI")

cnData <- correctReadDepth(tumWig = tumWig,
                          normWig = normWig,
                          gcWig = gcWig,
                          mapWig = mapWig,
                          genomeStyle = 'NCBI')

logR <- getPositionOverlap(data$chr, data$posn, cnData)

data$logR <- log(2^logR)

data <- filterData(data, c(1:22, "X", "Y"), minDepth = 10, maxDepth = 200,
                   positionList = NULL, centromere = NULL, centromere.flankLength = 10000)

numClusters <- nClusters
params <- loadDefaultParameters(copyNumber = 5,
                                numberClonalClusters = numClusters,
                                symmetric = TRUE, hetBaselineSkew = 0,
                                alleleEmissionModel = "binomial", data = data)

params$ploidyParams$phi_0 <- 2 # for diploid

convergeParams <- runEMclonalCN(data, params,
                                maxiter = 3, maxiterUpdate = 50,
                                useOutlierState = FALSE, txnExpLen = 1e15,
                                txnZstrength = 5e5,
                                normalEstimateMethod = "map",
                                estimateS = TRUE, estimatePloidy = TRUE)

optimalPath <- viterbiClonalCN(data, convergeParams)

results <- outputTitanResults(data, convergeParams, optimalPath,
                              filename = NULL, posteriorProbs = FALSE,
                              subcloneProfiles = TRUE, correctResults = TRUE,
                              proportionThreshold = 0.05,
                              proportionThresholdClonal = 0.05,
                              is.haplotypeData = FALSE)

convergeParams <- results$convergeParam ## use corrected parameters
results <- results$corrResults ## use corrected results

segs <- outputTitanSegments(results, id = "test", convergeParams,
                            filename = NULL, igvfilename = NULL)
# Plotting
par(mfrow=c(2,1))
#logR
ploidy <- tail(convergeParams$phi, 1)
normal <- tail(convergeParams$n, 1)
plotCNlogRByChr(results, segs = segs, chr = 5, ploidy = ploidy, normal = normal,
                ylim = c(-2, 2), cex = 0.25, xlab = "", main = "Chr 5")
#CN Plot
plotSegmentMedians(segs, chr=5, resultType = "LogRatio", plotType = "CopyNumber",
                   plot.new = TRUE, ylim = c(0, 4), main="Chr 5")

segs <- segs[!is.na(segs$Cellular_Frequency),]
sizes <- segs$End_Position.bp.-segs$Start_Position.bp.
segs <- segs[sizes>0,]

cols <- colnames(segs)
cols[c(3,4,15)] <- c('Start_Position(bp)','End_Position(bp)','Clonal_Frequency')
colnames(segs) <- cols

write.table(segs,outFile,sep='\t',row.names = F,col.names = T,quote = F)

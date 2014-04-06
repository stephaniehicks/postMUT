#!/usr/bin/Rscript --max-ppsize=150000
#
# Example: .scripts/postMUT.R <filename>


############################################################################
## Set parameters for estimation
##
## N.cores = Number of cores available to speed up estimation (e.g. 1, 4, 8)
## total.sims = Total number of simulations to run (must be divisible by N.cores)
## n.random.start.xy = Number of random starts in postMUT (simple)
## n.random.start.xyz = Number of random starts in postMUT
##
## N.sims = DO NOT CHANGE (Number of simulations to run per core)
## N.algo = DO NOT CHANGE (Number of in silico methods (SIFT, PPH2, Xvar -> N.algo = 3))
## N.par.xy = DO NOT CHANGE (Number of parameters in postMUT (simple) model)
## N.par.xyz = DO NOT CHANGE (Number of parameters in postMUT model) 
############################################################################

N.cores = 4
total.sims = 100
n.random.start.xy = 10
n.random.start.xyz = 100

# Do not change the following: 
N.sims = total.sims / N.cores 
N.algo = 3
N.par.xy = 2*N.algo + 1 
N.par.xyz = 2*N.algo + 3

############################################################################
## Load R packages, outside functions and create file names
##
## R packages: ggplot2, parallel
############################################################################

library(ggplot2)
library(parallel)
source("scripts/functions_EMAlgo.R")

args <- commandArgs(T)
if(length(args) < 1)
	stop("The script needs 1 input arguments.")
input.file <- args[1]

### Create input and output file names
input.file.name <- paste(input.file, "postMUT-in.txt", sep="-")
output.file.name <- paste(input.file, "postMUT-out.txt", sep="-")
output.file.name.BM <- paste(input.file, "BM-postMUT-out.txt", sep="-")
output.file.name.ABM <- paste(input.file, "ABM-postMUT-out.txt", sep="-")
output.file.name.ABMab <- paste(input.file, "ABMab-postMUT-out.txt", sep="-")
figure.file.name <- paste(input.file, "postMUT-plot.pdf", sep="-")

start.file <- "input_files/postMUT_input"
end.file.parameter <- "output_files/postMUT_output/par_est"
end.file.post <- "output_files/postMUT_output/post_est"
end.file.post.all <- "output_files/postMUT_output/post"
end.file.figure <- "output_files/postMUT_output/plots"


############################################################################
## Import data, clean data, prepare data for estimation
##
## Column names 1-5: CHROM, POS, REF_AA, MUT_AA, GENOTYPE
## Column names 6-7: SIFT_pred, SIFT_score
## Column names 8-9: Xvar_pred, Xvar_score
## Column names 10-11: PPH2_HD_pred, PPH2_HD_score
############################################################################

data <- read.table(paste(start.file, input.file.name, sep="/"), header = TRUE, sep="\t")

# Clean up predictions
data$SIFT_pred <- ifelse(data$SIFT_pred == "DAMAGING",1, ifelse(data$SIFT_pred == "TOLERATED", 0, NA))
data$PPH2_HD_pred <- ifelse(data$PPH2_HD_pred == "benign",0, ifelse(data$PPH2_HD_pred == "possiblydamaging" | data$PPH2_HD_pred == "probablydamaging", 1, NA))
data$Xvar_pred <- ifelse(data$Xvar_pred == "high" | data$Xvar_pred == "medium", 1, ifelse(data$Xvar_pred == "low" | data$Xvar_pred == "neutral", 0, NA))

# Combine predictions and filter for rows (i.e. mutations) with predictions from all three in silico methods
data_pred = data.frame(data[,1:5],"SIFT" = data$SIFT_pred, "Xvar" = data$Xvar_pred, "PPH2_HD" = data$PPH2_HD_pred, "SIFT_score" = data$SIFT_score, "Xvar_score" = data$Xvar_score, "PPH2_HD_score" = data$PPH2_HD_score)
data_pred = na.omit(data_pred)




############################################################################
## Estimate postMUT (simple) and postMUT parameters
##
############################################################################

N.size = nrow(data_pred[,6:8]); 
sum.group <- group.counts(data_pred[,6:8], N.algo)


simCR <- function(N.sims, sum.group, N.par.xy, N.par.xyz, N.algo, n.random.start.xy, n.random.start.xyz){
	source("scripts/functions_EMAlgo.R")
	EM.Est.xy = array(0, dim = c(N.par.xy, N.sims))
	EM.Est.xyz = array(0, dim = c(N.par.xyz, N.sims))
	EM.Est.xyz.ab = array(0, dim = c(N.par.xyz - 2, N.sims))

  for(k in 1:N.sims){
     repeat{
      output.xy <- capture.em.xy(sum.group, random.start = TRUE, n.algos = N.algo, ab_1 = FALSE, N.random.start = n.random.start.xy)
      if(length(unlist(output.xy)) == N.par.xy){
        if(all(output.xy[[1]] <= 0.5) & all(output.xy[[2]] >= 0.5)){ break }
      }
    }

   	if(length(unlist(output.xy)) == N.par.xy){
        EM.Est.xy[,k] <- unlist(output.xy)
    } else { 
        EM.Est.xy[,k] <- rep(NA,N.par.xy) 
    }

    repeat{
      output.xyz <- capture.em.xyz(sum.group, random.start = TRUE, n.algos = N.algo, de_1 = FALSE, N.random.start = n.random.start.xyz)
      if(length(unlist(output.xyz)) == N.par.xyz){
        if(all(a_b(output.xyz)[[1]] <= 0.5) & all(a_b(output.xyz)[[2]] >= 0.5) & unlist(output.xyz)[length(unlist(output.xyz)) - 1] > 0.75 ){ break }
      } 
    }

    if(length(unlist(output.xyz)) == N.par.xyz){
      EM.Est.xyz[,k] <- unlist(output.xyz)
      EM.Est.xyz.ab[,k] <- c(unlist(a_b(output.xyz)), unlist(output.xyz)[length(unlist(output.xyz))])
    } else { 
      EM.Est.xyz[,k] <- rep(NA,N.par.xyz)
      EM.Est.xyz.ab[,k] <- rep(NA,N.par.xyz - 2)
    }
 
  } # closes k for loop
  
  out <- rbind(EM.Est.xy, EM.Est.xyz, EM.Est.xyz.ab)
  return(out)
} # closes function



### Invoke use of multiple cores
cl <- makeCluster(N.cores)
sim.kk <- do.call("cbind", clusterCall(cl, simCR, N.sims, sum.group, N.par.xy, N.par.xyz, N.algo, n.random.start.xy, n.random.start.xyz) )
stopCluster(cl)

summary.data <- data.frame("mean" = apply(sim.kk, 1, mean), "sd" = apply(sim.kk, 1, sd))

BM.mat <- summary.data[1:N.par.xy,]
ABM.mat <- summary.data[(N.par.xy + 1):(N.par.xy + N.par.xyz),]
ABM.ab.mat <- summary.data[(N.par.xy + N.par.xyz + 1): (N.par.xy + N.par.xyz + N.par.xy),]


############################################################################
## WRITE Sensitivity and Specificity estimates to output file
##
############################################################################

write.table(BM.mat, file = paste(end.file.parameter, output.file.name.BM, sep="/"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(ABM.mat, file = paste(end.file.parameter, output.file.name.ABM, sep="/"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(ABM.ab.mat, file = paste(end.file.parameter, output.file.name.ABMab, sep="/"), sep="\t", quote = FALSE, row.names = FALSE)

############################################################################
## Calculate Posterior Probabilities
##
############################################################################

post_short = data.frame(pred.mat(N.algo), "postMUT-simple" = round(posterior.xy(pred.mat(N.algo), BM.mat[,1], n.algos = N.algo),3), "postMUT" = round(posterior.xyz(pred.mat(N.algo), ABM.mat[,1], n.algos = N.algo), 3)) 
write.table(post_short, file = paste(end.file.post, output.file.name, sep="/"), sep="\t", quote = FALSE, row.names = FALSE)

post.out.xy = posterior.xy(data_pred[,6:8], BM.mat[,1], n.algos = N.algo) # Posterior Pr(D)
post.out.xyz = posterior.xyz(data_pred[,6:8], ABM.mat[,1], n.algos = N.algo) # Posterior Pr(D)
data_new = data.frame(data_pred, "postMUT-simple" = post.out.xy, "postMUT" = post.out.xyz)
write.table(data_new, file = paste(end.file.post.all, output.file.name, sep="/"), sep="\t", quote = FALSE, row.names = FALSE)


############################################################################
## Plot Sensitivity and Specificity parameter estimates
##
############################################################################

wc_EM.xy <- data.frame(Sensitivity = BM.mat[4:6,1], Specificity = 1-BM.mat[1:3,1])
wc_EM.xyz <- data.frame(Sensitivity = ABM.ab.mat[4:6,1], Specificity = 1-ABM.ab.mat[1:3,1])
cool <- data.frame(rbind(wc_EM.xy,wc_EM.xyz), Algorithm =c("postMUT (simple)","postMUT (simple)","postMUT (simple)","postMUT","postMUT","postMUT"))

pdf(file=paste(end.file.figure, figure.file.name, sep="/"))
plot(1-cool[1,2],cool[1,1], col = 1, xlim = c(0,1), ylim=c(0,1.01), main = " ", xlab = "False Positive Rate", ylab = "True Positive Rate", pch = 16, cex = 1.5, cex.axis = 1.2, cex.lab = 1.2)
points(1-cool[2,2],cool[2,1], col = 1, pch = 18, cex = 1.5)
points(1-cool[3,2],cool[3,1], col = 1, pch = 17, cex = 1.5)

points(1-cool[4,2],cool[4,1], col = 4, pch = 16, cex = 1.5)
points(1-cool[5,2],cool[5,1], col = 4, pch = 18, cex = 1.5)
points(1-cool[6,2],cool[6,1], col = 4, pch = 17, cex = 1.5)

abline(0,1,lty=2,col=80)
legend(.64,.42,c("SIFT", "MutationAsessor","PolyPhen-2"),col=c(1,1,1),pch=c(1,2,5),cex=1.2)
legend(.56,.20,c("postMUT (simple)", "postMUT"),col=c(1,4),lty=c(1,1),lwd=c(4,4),cex=1.2)

dev.off()


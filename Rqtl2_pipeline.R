library(qtl2)
# read de yaml file
daturaf2 <- read_cross2("controlFile.yaml")
summary(daturaf2)
# calculating genotypes probabilities
pr <- calc_genoprob(daturaf2, map=daturaf2$map, error_prob=1.0e-4)
# Recall that the result of calc_genoprob, pr, is a list of three-dimensional arrays (one per chromosome).
names(pr)
# Each 3d array of probabilities is arranged as individuals × genotypes × positions. Have a look at the names of each of the three dimensions for chromosome 19.
dimnames(pr$`12`)
# Performing genome scan
#Firts set covariate
covariate <- read_pheno("/Users/mijaildelacruz/Desktop/newgenotypes_mareymaps/phenotypes_rqtl2.csv", phenocovarfile = "/Users/mijaildelacruz/Desktop/newgenotypes_mareymaps/Phenotypes_covar.csv")
# run scan
out <- scan1(genoprobs = pr, pheno = daturaf2$pheno, addcovar = covariate$POP)
# Check the LOD scores for phenotypes
head(out, n=10)
# The function plot_scan1() can be used to plot the LOD curves. Use the argument lodcolumn to indicate which column to plot.
plot_scan1(out, map = daturaf2$gmap, lodcolumn = "AverageResistance")
# Perform a permutation test
operm <- scan1perm(genoprobs = pr, pheno = daturaf2$pheno, addcovar = covariate$POP, n_perm = 1000)
head(operm, n=10)
# To get estimated significance thresholds, use the function summary().
summary(operm)
# To find peaks above a given threshold LOD value, use the function find_peaks(). It can also provide Bayesian credible intervals by using the argument prob (the nominal coverage for the Bayes credible intervals). Set the argument expand2markers = FALSE to keep from expanding the interval out to typed markers, or exclude this argument if you’d like to include flanking markers.
# Firts set a variable for threshold
thr = summary(operm)
find_peaks(scan1_output = out, map = daturaf2$gmap, threshold = thr, prob = 0.95, expand2markers = FALSE)
### NOW TESTING ANOTHER MODEL USING A KINDSHIP MATRIX ###
# Calculating A Kinship Matrix, Linear mixed models (LMMs) consider genome-wide similarity between all pairs of individuals to account for population structure, known kinship and unknown relatedness. 
kinship <- calc_kinship(probs = pr)
kinship[1:5, 1:5]
out_kin<-scan1(genoprobs = pr, pheno = daturaf2$pheno, kinship = kinship, addcovar = covariate$POP)
## PLOTTING BOTH MODELS (scanone without kindship and with kindship matrix)
plot_scan1(out, map = daturaf2$gmap, lodcolumn = "AverageResistance", col = "green", add= TRUE)
plot_scan1(out_kin, map = daturaf2$gmap, lodcolumn = "AverageResistance", col = "blue", add = TRUE)
# Testing for QTL effect
c6eff <- scan1coef(pr[,"6"], daturaf2$pheno[,"AverageResistance"])
dim(c6eff)
head(c6eff)
# PLOTTING QTL EFFECTS
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c6eff, map = daturaf2$gmap, columns=1:3, col=col)
last_coef <- unclass(c6eff)[nrow(c6eff),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
# Additive and dominance effects
c6effB <- scan1coef(pr[,"6"], daturaf2$pheno[,"AverageResistance"],
                    contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))
dim(c6effB)
head(c6effB)
# plotting additive and dominance effect
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c6effB, map = daturaf2$gmap, columns=2:3, col=col)
last_coef <- unclass(c6effB)[nrow(c6effB),2:3] # last two coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
# Plotting the raw phenotypes and genotypes for Average REsistance
g <- maxmarg(pr, daturaf2$gmap, chr=6, pos=47.884, return_char=TRUE)
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, daturaf2$pheno[,"AverageResistance"], ylab="Resistance phenotype", SEmult=2)
# Plotting the raw phenotypes and genotypes for Damage
g <- maxmarg(pr, daturaf2$gmap, chr=11, pos=34.578, return_char=TRUE)
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, daturaf2$pheno[,"D4"], ylab="D4 phenotype", SEmult=2)

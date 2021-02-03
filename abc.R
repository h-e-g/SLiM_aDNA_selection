## Copyright 2021 GEH Institut Pasteur, Gaspard Kerner (gakerner@pasteur.fr)

library(abc)

# Read output parameter after SLiM simulations
dir1 # Points to the file with parameter data (one line per simulation)
dir2 # Points to the file with frequency data (one line per simulation, columns are each time bin)
dir3 # Points to the file with ancestry data (one line per simulation, columns are each time bin)


age_data= read.table(dir1,header=T)
pp = read.table(dir2,header=F,sep="\t")
ppAA = read.table(dir3,header=F,sep="\t")

colnames(pp) = c("Generation_51500","Generation_43000","Generation_42500","Generation_38500","Generation_37000","Generation_34500",
                         "Generation_34000","Generation_33500","Generation_33000","Generation_32500",
                         "Generation_30000","Generation_18500","Generation_15000","Generation_13500",
                         "Generation_11500","Generation_11000","Generation_10700","Generation_10000","Generation_9500",
                         "Generation_9000","Generation_8500","Generation_8000","Generation_7500",
                         "Generation_7000","Generation_6500","Generation_6000","Generation_5500",
                         "Generation_5000","Generation_4500","Generation_4000","Generation_3500",
                         "Generation_3000","Generation_2500","Generation_2000","Generation_1500",
                         "Generation_1000","Generation_500","Generation_0","AFR","EUR","ASI","MIDD","CA")

colnames(ppAA)=c("Generation_51500","Generation_43000","Generation_42500","Generation_38500","Generation_37000","Generation_34500",
                 "Generation_34000","Generation_33500","Generation_33000","Generation_32500",
                 "Generation_30000","Generation_18500","Generation_15000","Generation_13500",
                 "Generation_11500","Generation_11000","Generation_10700","Generation_10000","Generation_9500",
                 "Generation_9000","Generation_8500","Generation_8000","Generation_7500",
                 "Generation_7000","Generation_6500","Generation_6000","Generation_5500",
                 "Generation_5000","Generation_4500","Generation_4000","Generation_3500",
                 "Generation_3000","Generation_2500","Generation_2000","Generation_1500",
                 "Generation_1000","Generation_500","Generation_0","EUR")

sampleSizesCounts=c(371,38,109,59,21,36,32,90,177,88,81,42,39,27,60,44,29,10,6,6,6,6,1,1,2,1,1,2,1,1,1,2,1,1,1,1,1,1)
sampleSizesCounts=rev(sampleSizesCounts)

# Gather time bins into culturally accepted epochs
pp_noslid_freq=pp
average_paleo = apply(pp_noslid_freq[,1:which(colnames(pp_noslid_freq)=="Generation_13500")],1,function(i) 
{sum(i[1:which(colnames(pp_noslid_freq)=="Generation_13500")]*sampleSizesCounts[1:which(colnames(pp_noslid_freq)=="Generation_13500")])/sum(sampleSizesCounts[1:which(colnames(pp_noslid_freq)=="Generation_13500")])})
pp_noslid_freq[,which(colnames(pp_noslid_freq)=="Generation_13500")]=average_paleo
average_mesol = apply(pp_noslid_freq,1,function(i) 
{sum(i[which(colnames(pp_noslid_freq)=="Generation_11500"):which(colnames(pp_noslid_freq)=="Generation_8500")]*sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_11500"):which(colnames(pp_noslid_freq)=="Generation_8500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_11500"):which(colnames(pp_noslid_freq)=="Generation_8500")])})
pp_noslid_freq[,which(colnames(pp_noslid_freq)=="Generation_8500")]=average_mesol
average_earlyneol = apply(pp_noslid_freq,1,function(i) 
{sum(i[which(colnames(pp_noslid_freq)=="Generation_8000"):which(colnames(pp_noslid_freq)=="Generation_6500")]*sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_8000"):which(colnames(pp_noslid_freq)=="Generation_6500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_8000"):which(colnames(pp_noslid_freq)=="Generation_6500")])})
pp_noslid_freq[,which(colnames(pp_noslid_freq)=="Generation_6500")]=average_earlyneol
average_lateneol = apply(pp_noslid_freq,1,function(i) 
{sum(i[which(colnames(pp_noslid_freq)=="Generation_6000"):which(colnames(pp_noslid_freq)=="Generation_4500")]*sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_6000"):which(colnames(pp_noslid_freq)=="Generation_4500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_6000"):which(colnames(pp_noslid_freq)=="Generation_4500")])})
pp_noslid_freq[,which(colnames(pp_noslid_freq)=="Generation_4500")]=average_lateneol
average_bronze = apply(pp_noslid_freq,1,function(i) 
{sum(i[which(colnames(pp_noslid_freq)=="Generation_4000"):which(colnames(pp_noslid_freq)=="Generation_3000")]*sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_4000"):which(colnames(pp_noslid_freq)=="Generation_3000")])/sum(sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_4000"):which(colnames(pp_noslid_freq)=="Generation_3000")])})
pp_noslid_freq[,which(colnames(pp_noslid_freq)=="Generation_3000")]=average_bronze
average_iron = apply(pp_noslid_freq,1,function(i) 
{sum(i[which(colnames(pp_noslid_freq)=="Generation_2500"):which(colnames(pp_noslid_freq)=="Generation_2000")]*sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_2500"):which(colnames(pp_noslid_freq)=="Generation_2000")])/sum(sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_2500"):which(colnames(pp_noslid_freq)=="Generation_2000")])})
pp_noslid_freq[,which(colnames(pp_noslid_freq)=="Generation_2000")]=average_iron
average_middle = apply(pp_noslid_freq,1,function(i) 
{sum(i[which(colnames(pp_noslid_freq)=="Generation_1500"):which(colnames(pp_noslid_freq)=="Generation_500")]*sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_1500"):which(colnames(pp_noslid_freq)=="Generation_500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_freq)=="Generation_1500"):which(colnames(pp_noslid_freq)=="Generation_500")])})
pp_noslid_freq[,which(colnames(pp_noslid_freq)=="Generation_500")]=average_middle
pp_noslid2_freq=pp_noslid_freq[,c(which(colnames(pp_noslid_freq)=="Generation_13500"),which(colnames(pp_noslid_freq)=="Generation_8500"),
                                  which(colnames(pp_noslid_freq)=="Generation_6500"),which(colnames(pp_noslid_freq)=="Generation_4500"),which(colnames(pp_noslid_freq)=="Generation_3000"),
                                  which(colnames(pp_noslid_freq)=="Generation_2000"),which(colnames(pp_noslid_freq)=="Generation_500"),
                                  which(colnames(pp_noslid_freq)=="EUR"))]
colnames(pp_noslid2_freq)=c("Paleo","Mesol","EarlyNeol","LateNeol","Bronze","Iron","Middle","PresEUR")

pp_noslid_AA=ppAA
average_paleo = apply(pp_noslid_AA[,1:which(colnames(pp_noslid_AA)=="Generation_13500")],1,function(i) 
{sum(i[1:which(colnames(pp_noslid_AA)=="Generation_13500")]*sampleSizesCounts[1:which(colnames(pp_noslid_AA)=="Generation_13500")])/sum(sampleSizesCounts[1:which(colnames(pp_noslid_AA)=="Generation_13500")])})
pp_noslid_AA[,which(colnames(pp_noslid_AA)=="Generation_13500")]=average_paleo
average_mesol = apply(pp_noslid_AA,1,function(i) 
{sum(i[which(colnames(pp_noslid_AA)=="Generation_11500"):which(colnames(pp_noslid_AA)=="Generation_8500")]*sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_11500"):which(colnames(pp_noslid_AA)=="Generation_8500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_11500"):which(colnames(pp_noslid_AA)=="Generation_8500")])})
pp_noslid_AA[,which(colnames(pp_noslid_AA)=="Generation_8500")]=average_mesol
average_earlyneol = apply(pp_noslid_AA,1,function(i) 
{sum(i[which(colnames(pp_noslid_AA)=="Generation_8000"):which(colnames(pp_noslid_AA)=="Generation_6500")]*sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_8000"):which(colnames(pp_noslid_AA)=="Generation_6500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_8000"):which(colnames(pp_noslid_AA)=="Generation_6500")])})
pp_noslid_AA[,which(colnames(pp_noslid_AA)=="Generation_6500")]=average_earlyneol
average_lateneol = apply(pp_noslid_AA,1,function(i) 
{sum(i[which(colnames(pp_noslid_AA)=="Generation_6000"):which(colnames(pp_noslid_AA)=="Generation_4500")]*sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_6000"):which(colnames(pp_noslid_AA)=="Generation_4500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_6000"):which(colnames(pp_noslid_AA)=="Generation_4500")])})
pp_noslid_AA[,which(colnames(pp_noslid_AA)=="Generation_4500")]=average_lateneol
average_bronze = apply(pp_noslid_AA,1,function(i) 
{sum(i[which(colnames(pp_noslid_AA)=="Generation_4000"):which(colnames(pp_noslid_AA)=="Generation_3000")]*sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_4000"):which(colnames(pp_noslid_AA)=="Generation_3000")])/sum(sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_4000"):which(colnames(pp_noslid_AA)=="Generation_3000")])})
pp_noslid_AA[,which(colnames(pp_noslid_AA)=="Generation_3000")]=average_bronze
average_iron = apply(pp_noslid_AA,1,function(i) 
{sum(i[which(colnames(pp_noslid_AA)=="Generation_2500"):which(colnames(pp_noslid_AA)=="Generation_2000")]*sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_2500"):which(colnames(pp_noslid_AA)=="Generation_2000")])/sum(sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_2500"):which(colnames(pp_noslid_AA)=="Generation_2000")])})
pp_noslid_AA[,which(colnames(pp_noslid_AA)=="Generation_2000")]=average_iron
average_middle = apply(pp_noslid_AA,1,function(i) 
{sum(i[which(colnames(pp_noslid_AA)=="Generation_1500"):which(colnames(pp_noslid_AA)=="Generation_500")]*sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_1500"):which(colnames(pp_noslid_AA)=="Generation_500")])/sum(sampleSizesCounts[which(colnames(pp_noslid_AA)=="Generation_1500"):which(colnames(pp_noslid_AA)=="Generation_500")])})
pp_noslid_AA[,which(colnames(pp_noslid_AA)=="Generation_500")]=average_middle
pp_noslid2_AA=pp_noslid_AA[,c(which(colnames(pp_noslid_AA)=="Generation_13500"),which(colnames(pp_noslid_AA)=="Generation_8500"),
                              which(colnames(pp_noslid_AA)=="Generation_6500"),which(colnames(pp_noslid_AA)=="Generation_4500"),which(colnames(pp_noslid_AA)=="Generation_3000"),
                              which(colnames(pp_noslid_AA)=="Generation_2000"),which(colnames(pp_noslid_AA)=="Generation_500"),
                              which(colnames(pp_noslid_AA)=="EUR"))]
colnames(pp_noslid2_AA)=c("Paleo","Mesol","EarlyNeol","LateNeol","Bronze","Iron","Middle","PresEUR")

# Use these parameters, frequencies and ancestries to perform abc estimations

touse=which((400-age_data[,10])*250>500 & (400-age_data[,10])*250<10000) # use only simulations with an onset of selection between 500 and 10,000 years ago
summstat_sim=cbind(pp_noslid2[touse,],pp_noslid2_AA[touse,]) # simulated summary statistics (frequencies + ancestries)
summstat_emp=cbind(emp_noslid2,emp_noslid2_AA) # empirical summary statistics (input real values here)
age_param=age_data[touse,10] # simulated onset of selection
sel_param=age_data[touse,2] # simulated selection coefficient
lin_nonTol <- abc(target=summstat_emp, param=data.frame(AGE=age_param,SEL=sel_param), sumstat=summstat_sim, tol=1000/nrow(pp_noslid2), hcorr = FALSE, 
                     method = "loclinear", transf=c("none") ) # abc analysis

lin_nonTol$adj.values[,1] = (400-lin_nonTol$adj.values[,1])*250
est_nonTol=estimate_mode( lin_nonTol$adj.values[,1]) # point estimation for the onset of selection
IC_min_nonTol=quantile( (lin_nonTol$adj.values[,1]) , 0.025 ) # Lower bound of the confidence interval for the onset of selection
IC_max_nonTol=quantile( (lin_nonTol$adj.values[,1]) , 0.975 ) # Upper bound of the confidence interval for the onset of selection

est_nonTol=estimate_mode( lin_nonTol$adj.values[,2]) # point estimation for the selection coefficient
IC_min_nonTol=quantile( (lin_nonTol$adj.values[,2]) , 0.025 ) # Lower bound of the confidence interval for the selection coefficient
IC_max_nonTol=quantile( (lin_nonTol$adj.values[,2]) , 0.975 ) # Upper bound of the confidence interval for the selection coefficient




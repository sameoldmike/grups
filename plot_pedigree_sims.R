
# usage:  Rscript scriptname.R 
#             data_dir=path/ 
#             regex=str 
#             [max=float min=float] 
#             [plotval=float | float,float]
#             [plotdist=path/regex]         # ex: plotdist=./analyses/block_jackknife/PWD_from_stdin_2015-09-22.py.only_targets2.only_SNPs.only_transversions.mindepth2_2.ESP2MDR.ESP29MDR.blocksize10M.*.out
#             [label=str] 
#             [range=float]         [1.5]
#             [alpha=float]         [0.01]
#             [noPrint=int]                 # index of the relationship not to print (e.g. 1 for inbred-inbred relationship)
#             [w=float]                     # plot width
#             [h=float]                     # plot height
#             [nperms=int]                  # number of MC permutations to perform 
#             [heatMap]                     # also print a heatmap
#             [BCstringency=int]    [0]
#             [maxSample=int]       [inf]
#             [printAllReps]                # print the mean genetic distance for every replicate to the screen
#


#library(perm)
library(reshape)
library(ggplot2)


# parse command-line arguments
args <- commandArgs(TRUE)
print(args)

plotmax = -1
plotmin = -1
plotval = -1
plotval_std = -1
plotdist = ''
label = '.'
range = 0
violin = 1
boxplot_border_color = "white"
boxplot_filled_color = "white"
alpha = 0.01
noPrint = -1
w = 5 
h = 5
nperms = 10000
heatMap = 0
BCstringency = 0
maxSample=-1
printAllReps = 0
for (i in seq(1, length(args))){
    t = strsplit(args[i], '=')
    if (t[[1]][1] == 'regex'){
        regex = t[[1]][2]
    } else if (t[[1]][1] == 'data_dir'){
        cmdarg = 'data_dir='
        data_dir = t[[1]][2]
        #print(paste(cmdarg, t[[1]][2], sep=''))
    } else if (t[[1]][1] == 'max'){
        plotmax = as.numeric(t[[1]][2])
    } else if (t[[1]][1] == 'min'){
        plotmin = as.numeric(t[[1]][2])
    } else if (t[[1]][1] == 'noPrint'){
        temp = t[[1]][2]
        temp = strsplit(t[[1]][2], ',')
        noPrint = as.numeric(temp[[1]][1])
        for (Y in seq(2, length(temp[[1]]))){ 
            noPrint = c(noPrint, as.numeric(temp[[1]][Y]))
        }
    } else if (t[[1]][1] == 'maxSample'){
        maxSample = as.numeric(t[[1]][2])     
    } else if (t[[1]][1] == 'BCstringency'){
        BCstringency = as.numeric(t[[1]][2])      
    } else if (t[[1]][1] == 'plotval'){
        temp = t[[1]][2]
        temp = strsplit(t[[1]][2], ',')
        plotval = as.numeric(temp[[1]][1])
        window_size = plotval*0.05
        if (length(temp[[1]]) == 2){
            plotval_std = as.numeric(temp[[1]][2])
            window_size = plotval_std
        }
    } else if (t[[1]][1] == 'plotdist'){
        temp = t[[1]][2]
        # this saves a list of absolute file paths to temp2
        temp2 = list.files(path=dirname(temp), pattern=glob2rx(basename(temp)), full.names=TRUE)
        for (Y in seq(1, length(temp2))){
            temp3 = na.omit(read.table(temp2[Y], header=FALSE))
            if (Y == 1){
                plotdist = temp3
            }
            else{
                plotdist = c(plotdist, temp3)
            }
        }
    } else if (t[[1]][1] == 'range'){
        range = as.numeric(t[[1]][2])
    } else if (t[[1]][1] == 'nperms'){
        nperms = as.numeric(t[[1]][2])
    } else if (t[[1]][1] == 'w'){
        w = as.numeric(t[[1]][2])
    } else if (t[[1]][1] == 'h'){
        h = as.numeric(t[[1]][2])
    } else if (t[[1]][1] == 'alpha'){
        alpha = as.numeric(t[[1]][2])
    } else if (t[[1]][1] == 'label'){
        label = paste('.', t[[1]][2], '.', sep="")
    } else if (t[[1]][1] == 'heatMap'){
        heatMap = 1
        library(gplots)
        library(RColorBrewer)
        library(vegan)
    } else if (t[[1]][1] == 'printAllReps'){
        printAllReps = 1
    } else{
        cat('Error parsing command-line args\n')
    }
}


# parse file names
all_files = list.files(path=data_dir, pattern=glob2rx(regex))
print(glob2rx(regex))
print(all_files)
cat(paste('Number of files to process: ', length(all_files), '\n'))


# process and print selected files
for (gg in seq(1, length(all_files))){
    filename_base = paste(data_dir, substr(all_files[gg], 1, nchar(all_files[gg])-4), sep="")
    runID = strsplit(filename_base, "_")
    runID = runID[[1]][length(runID[[1]])-1]
    cat(paste('reading ', filename_base, ".out\n", sep=''))
    data                = read.table(file=paste(filename_base, ".out", sep=''), header=F)
    cat(paste('reading ', filename_base, ".numvars\n", sep=''))
    data_numvars        = read.table(file=paste(filename_base, ".numvars", sep=''), header=F)
    cat(paste('reading ', filename_base, ".numSNPscov\n", sep=''))
    data_SNPscov        = read.table(file=paste(filename_base, ".numSNPscov", sep=''), header=F)
    cat(paste('reading ', filename_base, ".labels\n", sep=''))
    data_labels         = read.table(file=paste(filename_base, ".labels", sep=''), header=F)
    
    numreps = floor(length(data[,1]) / 22)
    cat(paste('Number of replicates detected in file', gg, '=', numreps, '\n'))
    # randomly sample replicates for analysis
    if (maxSample != -1){
        if (numreps >= maxSample){
            to_sample = maxSample
        } else{
            to_sample = numreps
        }
    } else{
        to_sample = numreps
    }
    
    label = paste(label, 'reps', to_sample, sep = '')
    label = paste(label, '.w', w, '.h', h, '.', sep='')
    

    # parse files and perform analysis
    if (numreps < 1){
        cat('Input file does not contain a whole replicate. Skipping.\n')
    } 
    if (numreps > 1){
        cat(paste('Randomly sampling ', to_sample, ' of ', numreps, ' reps for analysis\n', sep=''))
        
        data2 = matrix(ncol=length(data[1,]), nrow=to_sample)
        sums_keeper = matrix(ncol=length(data[1,]), nrow=to_sample)
        sums_SNPscov_keeper = matrix(ncol=length(data[1,]), nrow=to_sample)
        pwdiff_expectations = matrix(ncol=length(data[1,]), nrow=to_sample)
        # extract total number of variants used in simulations (assumes it's the same in each relationship)
        numvar_sum = 0
        for (row in seq(1, 22)){
            numvar_sum = numvar_sum + data_numvars[row, 1]
        }

        # calculate genome-wide, mean pwdiff of all 22 chromsomes for each relationship
        rep = 1
        possible_blockstarts = seq(1, numreps*22, by=22)
        sampled_blockstarts = sample(possible_blockstarts, to_sample, replace = FALSE)
        for (blockstart in sampled_blockstarts){
            for (relationship in 1:length(data[1,])){
                chrom = 1   
                sums = 0
                sums_SNPscov = 0
                sums_numvars = 0
                sum_length = 0
                for (element in seq(blockstart, blockstart+21)){
                    sums          = sums         + (data[element,relationship])
                    sums_SNPscov  = sums_SNPscov + data_SNPscov[element,relationship]
                    sums_numvars  = sums_numvars + data_numvars[element,relationship]
                    sum_length = sum_length + 1
                    chrom = chrom + 1
                }
                data2[rep,relationship]                 = sums / sums_SNPscov
                sums_keeper[rep,relationship]           = sums
                sums_SNPscov_keeper[rep,relationship]   = sums_SNPscov
            }
            rep = rep + 1
        }
    
        
    
        # parse labels, define relationships and printing order
        labels = as.matrix(data_labels[1,])
        print(paste("parsed sims label:", labels, sep=" "))
        potential_labels_4 = t(as.matrix(c("inbred-inbred", "father-father", "father-child1", "child1-child2", "mother-cousin", "child2-gchild", "gchild-halfsib", "cousin-gchild", "father-mother")))
        labels_relationships = -1
        
        labels_printorder = seq(1, length(labels))
        labels_relationships = data_labels
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "white", "#D55E00", "#CC79A7", "black") # need to handle different label lengths
        labels_relationshipsTOPLOT = labels_relationships
        labels_printorderTOPLOT = labels_printorder
        if (noPrint != -1){
            labels_printorderTOPLOT = labels_printorder[-noPrint]
        } 
        print(paste("assigned rels label:", labels_relationships, sep=" "))
        print(paste("TOPLOT rels label:", labels_relationshipsTOPLOT[labels_printorderTOPLOT], sep=" "))
    

    
        # print some stats for comparison with observed data
        colMeans_sums_keeper            = colMeans(sums_keeper)
        colMeans_sums_SNPscov_keeper    = colMeans(sums_SNPscov_keeper)
        colMeans_data2                  = colMeans(data2)
        
        if (plotval == -1){
            cat(paste("\nRelationship", "MeanPWD/site", "MeanPWDs", "MeanOverlaps\n", sep="\t"))
            for (t in seq(1, length(labels))){
                cat(paste(labels_relationships[t], colMeans_data2[t], colMeans_sums_keeper[t], colMeans_sums_SNPscov_keeper[t], "\n", sep=" "))
            }
        }
        else{
            cat(paste("\nRelationship", "MeanPWD/site", "MeanPWDs", "MeanOverlaps", "Error", "numreps\n", sep="\t"))
            for (t in seq(1, length(labels))){
                cat(paste(labels_relationships[t], colMeans_data2[t], colMeans_sums_keeper[t], colMeans_sums_SNPscov_keeper[t], abs(plotval-colMeans_data2[t])/plotval, to_sample, "\n", sep="\t"))
            }
        }
        cat("\n")
        
        
        if (printAllReps == 1){
            cat(paste('Pairwise genetic distance from every replicate being considered:\n'))
            print(labels_relationships)
            print(data2)
        }
        

        # calculate Bhattacharyya coefficient for each pair of relationships
        cat("\n")
        cat("Matrix of Bhattacharyya coefficients for pairwise relationship distributions\n")
        BC_matrix = matrix(0, length(labels_relationships), length(labels_relationships))
        colnames(BC_matrix) <- labels_relationships
        rownames(BC_matrix) <- labels_relationships
        X = 1
        while (X <= length(labels_relationships)){
            Y = 1
            while (Y <= length(labels_relationships)){
                merged_data = c(data2[,X], data2[,Y])
                merged_data_range = max(merged_data)-min(merged_data)
                merged_data_hist = hist(merged_data, breaks=length(merged_data)/10, freq=FALSE)
                histX = hist(data2[,X], breaks=merged_data_hist$breaks, freq=FALSE)
                histY = hist(data2[,Y], breaks=merged_data_hist$breaks, freq=FALSE)
                histXPr = histX$counts / sum(histX$counts)
                histYPr = histY$counts / sum(histY$counts)
                numBCbins = length(histXPr)
                # estimate Bhattacharyya co-efficient
                for (i in seq(1, numBCbins)){
                    BC_matrix[X,Y] = BC_matrix[X,Y] + sqrt(histXPr[i]*histYPr[i])
                }
                Y = Y + 1
            } 
            X = X+1 
        }
        cat("\n")
        print(BC_matrix)



        # test normality of distributions
        kstest = rep(NA, length(labels_relationships))
        names(kstest) <- labels_relationships[labels_printorder]
        op <- options(warn = (-1)) # suppress warnings         
        for (t in 1:length(labels_relationships)){
            kstest[t] = ks.test(data2[,t], rnorm(length(data2[,t]), mean = mean(data2[,t]), sd = sd(data2[,t])))$p.val
        }
        options(op) # reset default value for warnings
        print(paste("Can a normal distribution be rejected at alpha = ", alpha, "?", sep = ""))
        print(kstest)
        print(kstest < alpha)
        cat("\n")
        
        
        # find shortest distance (z-score) from observation to mean of each relationship
        if (plotval != -1){
            obsDistZ = rep(NA, length(labels_relationships))
            names(obsDistZ) <- labels_relationships
            obsProb = rep(NA, length(labels_relationships))
            names(obsProb) <- labels_relationships
            for (t in 1:length(labels_relationships)){
                obsDistZ[t] = abs(plotval - mean(data2[,t]))/sd(data2[,t])
                obsProb[t] = pnorm(-obsDistZ[t])
            }
            print("Z-scores")
            print(obsDistZ)
            cat("\n")
            print("Probability of observation within simulated relationship distribution")
            print(obsProb)
            cat("\n")
            print(paste("Most likely relationship Z-score:", obsDistZ[which(min(obsDistZ) == obsDistZ)], labels_relationships[which(min(obsDistZ) == obsDistZ)], sep=" "))
            print(paste("Odds ratio of obs within most likely relationship (", labels_relationships[which(min(obsDistZ) == obsDistZ)], ") vs. other relationship", sep=""))
            best_odds = obsProb[which(min(obsDistZ) == obsDistZ)]/(1-obsProb[which(min(obsDistZ) == obsDistZ)])
            for (t in 1:length(labels_relationships)){
                if (t != which(min(obsDistZ) == obsDistZ)){
                    these_odds = obsProb[t]/(1-obsProb[t])
                    cat(paste(labels_relationships[t], best_odds/these_odds, "\n", sep="\t"))
                }
            }
            cat("\n")
        
            # calculate matrix of odds ratios for all relationships
            ORs_matrix = matrix(NA, length(labels_relationships), length(labels_relationships))
            colnames(ORs_matrix) <- labels_relationships[labels_printorder]
            rownames(ORs_matrix) <- labels_relationships[labels_printorder]
            X = 1
            while (X <= length(labels_relationships)){
                Y = 1
                while (Y <= length(labels_relationships)){
                    ORs_matrix[X,Y] = (obsProb[X]/(1-obsProb[X])) / (obsProb[Y]/(1-obsProb[Y]))
                    Y = Y + 1
                } 
                X = X + 1 
            }
            cat("Odds ratios matrix\n")
            print(ORs_matrix)
            cat("\n")
        }
        
        # calculate ratios of variances for each distribution
        varratios_matrix = matrix(NA, length(labels_relationships), length(labels_relationships))
        colnames(varratios_matrix) <- labels_relationships[labels_printorder]
        rownames(varratios_matrix) <- labels_relationships[labels_printorder]
        X = 1
        while (X <= length(labels_relationships)){
            Y = X + 1
            while (Y <= length(labels_relationships)){
                varratios_matrix[X,Y] = var(data2[,X]) / var(data2[,Y])
                Y = Y+1
            } 
            X = X+1 
        }
        cat("Variance ratios matrix\n")
        print(varratios_matrix)

        # print mean and variance of each relationship
        cat("\nRelationship\tMean\tVariance\n")
        for (X in seq(1, length(labels_relationships))){
            cat(paste(labels_relationships[X], '\t', mean(data2[,X]), '\t', var(data2[,X]), '\n'))
        }
        cat("\n")
                

        # calculate overlapping for each relationship pair assuming normality
        overlap_areas = matrix(NA, length(labels_relationships), length(labels_relationships))
        colnames(overlap_areas) <- labels_relationships
        rownames(overlap_areas) <- labels_relationships
        OR1v2 = matrix(NA, length(labels_relationships), length(labels_relationships))
        colnames(OR1v2) <- labels_relationships
        rownames(OR1v2) <- labels_relationships
        PrOR1v2gt10 = matrix(NA, length(labels_relationships), length(labels_relationships))
        colnames(PrOR1v2gt10) <- labels_relationships
        rownames(PrOR1v2gt10) <- labels_relationships
        PrOR1v2gt100 = matrix(NA, length(labels_relationships), length(labels_relationships))
        colnames(PrOR1v2gt100) <- labels_relationships
        rownames(PrOR1v2gt100) <- labels_relationships
        X = 1
        while (X <= length(labels_relationships)){
            Y = 1
            while (Y <= length(labels_relationships)){
                u1 = mean(data2[,X])
                sd1 = sd(data2[,X])
                u2 = mean(data2[,Y])
                sd2 = sd(data2[,Y])
                R1 = ((u1*sd2*sd2) - (u2*sd1*sd1) + (sd1*sd2*sqrt(((u1-u2)*(u1-u2)) + (((sd2*sd2) - (sd1*sd1))*log((sd2*sd2)/(sd1*sd1))))))/((sd2*sd2)-(sd1*sd1))
                R2 = ((u1*sd2*sd2) - (u2*sd1*sd1) - (sd1*sd2*sqrt(((u1-u2)*(u1-u2)) + (((sd2*sd2) - (sd1*sd1))*log((sd2*sd2)/(sd1*sd1))))))/((sd2*sd2)-(sd1*sd1))
                X1 = min(R1, R2)
                X2 = max(R1, R2)
                # find overlapping area
                overlap_areas[X,Y]  = pnorm((X1 - u1)/sd1) + pnorm((X2 - u2)/sd2) - pnorm((X1 - u2)/sd2) - pnorm((X2 - u1)/sd1) + 1
                if (is.nan(overlap_areas[X,Y]) == FALSE){
                    if (overlap_areas[X,Y] > 1){
                        overlap_areas[X,Y] = NA
                    }
                }
                # find ORs for 100 randomly drawn values     
                obs = rnorm(100, mean = u1, sd = sd1)
                vals = rep(NA, length(obs))
                for (m in 1:length(obs)){
                    vals[m] = ((1 - pnorm((obs[m] - u1)/sd1)) / pnorm((obs[m] - u1)/sd1)) / (pnorm((obs[m] - u2)/sd2) / (1 - pnorm((obs[m] - u2)/sd2)))
                }
                OR1v2[X,Y] = mean(vals);
                PrOR1v2gt10[X,Y] = length(which(vals > 10)) / length(vals)
                PrOR1v2gt100[X,Y] = length(which(vals > 100)) / length(vals)
                Y = Y+1
            } 
            X = X+1 
        }
        cat("Matrixof overlap areas for relationship pairs\n")
        print(overlap_areas)
        cat("\n")
        cat("Matrix of Pr(OR > 100) for relationship pairs\n")
        print(PrOR1v2gt100)
        cat("\n")
        
        




        # PLOT HEATMAP

        # if making a heatmap of BC
        if (heatMap == 1){
            pdf(paste(filename_base, ".BC_heatmap", label, "pdf", sep=""))
            
            temp_BC_matrix = BC_matrix[labels_printorderTOPLOT,labels_printorderTOPLOT]
            temp_BC_matrix[lower.tri(temp_BC_matrix)] = NA
            diag(temp_BC_matrix) = NA
            if (BCstringency == 0){
                col_breaks = c(seq(0.00, 0.50, length=10),  # for color 1
                               seq(0.50, 0.95, length=9),   # for color 2
                               seq(0.95, 1.00, length=1))   # for color 3
                heatmap.2(temp_BC_matrix, density.info="none", trace="none", dendrogram='none', col=colorRampPalette(c("green", "yellow", "red"), space = "rgb")(n = length(col_breaks)-1), breaks=col_breaks, Rowv=FALSE, Colv=FALSE, sepcolor="white", colsep=1:length(labels_printorderTOPLOT), rowsep=1:length(labels_printorderTOPLOT), c(0.025, 0.025))            
            } else if (BCstringency == 1){
                col_breaks = c(0.000, 0.010, 0.050, 1.000) 
                myCol = c("green", "yellow", "red")
                heatmap.2(temp_BC_matrix, density.info="none", trace="none", dendrogram='none', key=FALSE, col=myCol, breaks=col_breaks, Rowv=FALSE, Colv=FALSE, sepcolor="white", colsep=1:length(labels_printorderTOPLOT), rowsep=1:length(labels_printorderTOPLOT), c(0.025, 0.025))
                legend("left", fill = myCol, legend = c("<1%", "1-5%", ">5%"))
            }
            dev.off()
        }
    
    

        # PLOT SIMULATION RESULTS

        # do a violin plot
        if (violin == 1){
            breaklabels = paste(rep('V',length(labels_printorderTOPLOT)), seq(1, length(labels_printorderTOPLOT)), sep='')
            GV = ggplot(data=melt(as.data.frame(data2[,labels_printorderTOPLOT])), aes(y=value, x=variable, fill=factor(variable))) + geom_violin(alpha = 0.5, scale="width") + scale_fill_manual(values=cbPalette[labels_printorderTOPLOT], breaks=breaklabels, labels=labels_relationshipsTOPLOT[labels_printorderTOPLOT]) + theme(legend.title=element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + labs(y="Pairwise differences per site")
            # incorporate plotmax and plotmin                                                                                                                                                                + scale_fill_discrete(breaks=breaklabels, labels=labels_relationshipsTOPLOT[labels_printorderTOPLOT])    
            
            # add a mini boxplot to show interquartile range
            GV = GV + geom_boxplot(width=0.15, outlier.size = 1, alpha=1.0)
            
            # plot directly oberserved values if provided
            if (plotval != -1){
                GV = GV + geom_hline(yintercept = plotval, size=0.5, linetype = 2)
            }
            if (plotval_std != -1){
                GV = GV + geom_hline(yintercept = plotval-plotval_std, size = 0.5)
                GV = GV + geom_hline(yintercept = plotval+plotval_std, size = 0.5)
            }
            if (plotdist != ''){
                plotdist_mean = mean(plotdist[,2])
                plotdist_std = sd(plotdist[,2])                
                GV = GV + geom_hline(yintercept = plotdist_mean, size=1, linetype = 2)
                GV = GV + geom_hline(yintercept = plotdist_mean-plotdist_std, size = 0.5)
                GV = GV + geom_hline(yintercept = plotdist_mean+plotdist_std, size = 0.5)
            }
            if (plotmin != -1){
                GV = GV + ylim(plotmin, plotmax)
            }
            
            cat(paste('plotting to: ', filename_base, ".results_plot", label, "pdf\n", sep=""))
            ggsave(paste(filename_base, ".results_plot", label, "pdf", sep=""), plot = GV, width = w, height = h)
        }

    }
}
cat(paste('...done!', '\n'))

        
packages_to_install = c("readxl", "writexl", "dplyr", "reshape2", "GGally", 
                        "gridExtra", "grid", "readr", "corrplot", "plot.matrix",
                        "aricode", "infotheo", "heatmap3", "pheatmap", "lattice", 
                        "NPRED", "csvread", "plotrix", "soiltexture", "stringr", 
                        "installr", "resemble", "prospectr", "magrittr", "doParallel", 
                        "parallel", "foreach", "ggplot2", "tidyr", "pls", "ChemoSpec", "MASS", "Johnson")
# Install and load packages
# install.packages(packages_to_install, dependencies = TRUE)
# Load the installed packages
sapply(packages_to_install, require, character.only = TRUE)

##FUNCTIONS
# Function to calculate the minimum indices
get_lowest_indices <- function(data, no_neighbor = 5) {
  min_indices <- apply(data, 2, function(col) sort(col, index.return = TRUE)$ix[1:no_neighbor])
  min_indices <- as.numeric(min_indices)
  unique_min_indices <- unique(min_indices)
  return(unique_min_indices)
}


# Function to calculate R², RMSE, RPD, and RPIQ
calculate_statistics <- function(observed, predicted) {
  r_squared <- cor(observed, predicted)^2
  rmse <- sqrt(mean((observed - predicted)^2))
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  deviation <- observed - mean_observed
  rpd <- sd(observed) / rmse
  iqr_observed <- IQR(observed)
  rpiq <- iqr_observed / rmse
  bias <- mean_observed - mean_predicted
  result <- list(R_squared = r_squared,RMSE = rmse,RPD = rpd,RPIQ = rpiq,Bias = bias)
  return(result)
}

calculate_statistics_result <- function(observed_values1,Y2_hat1, ncomp1) {
  cal_statistics_result <- NULL
  for (i in seq(1, ncomp1, by = 1)) {
    predicted_values1 <- Y2_hat1[, 1, i] 
    cal_statistics_result[[i]] <- calculate_statistics(observed_values1, predicted_values1)
  }
  return(cal_statistics_result)
}

#Stratified Function
stratified <- function(df, group, size, select = NULL, 
                       replace = FALSE, bothSets = FALSE) {
  if (is.null(select)) {
    df <- df
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(df)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) df[[x]] %in% select[[x]])
    df <- df[rowSums(temp) == length(select), ]
  }
  df.interaction <- interaction(df[group], drop = TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
      message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
              paste(n, collapse = ", "), "),\n.Names = c(",
              paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = replace), ])
  set1 <- do.call("rbind", temp)
  
  if (isTRUE(bothSets)) {
    set2 <- df[!rownames(df) %in% rownames(set1), ]
    list(SET1 = set1, SET2 = set2)
  } else {
    set1
  }
}

setwd("D:/Academics/PhD/SEM 11/Paper_3_Global2Local/Working_files/Refined/")
setwd("F:/CRT/Refined/")
setwd("C:/Refined/")
#OSSL library
OSSL = read_csv("OSSL.csv")

#reformatting the dataset to adjust for existing datastructure
foo = OSSL
foo1 = as.data.frame(foo[,2:8])
foo1$spc = -log10(data.matrix(foo[,60:length(foo)]))
rownames(foo1$spc) = as.character(seq(1, 64323, by = 1))
colnames(foo1$spc) = as.character(seq(502, 2500, by = 2))
rm(OSSL,foo)
OSSL = foo1
rm(foo1)

#Karnataka Data
foo = read_excel("Karnataka.xlsx", sheet = "Sheet1")
foo1 = as.data.frame(foo[,1:7])
foo2 = -log10(data.matrix(foo[,8:length(foo)]))
rownames(foo2) = as.character(seq(1, 497, by = 1))
colnames(foo2) = as.character(seq(400, 2500, by = 1))
#pre-process the spectra: resample it to a resolution of 2 nm 
old_wavs = foo2 %>% colnames() %>% as.numeric()
new_wavs = seq(502, 2500, by = 2)
foo1$spc = foo2 %>% 
  resample(wav = old_wavs, new.wav = new_wavs)
kar = foo1
rm(foo,foo1,foo2)

#Berambadi Data
foo = read_excel("Berambadi_IP.xlsx", sheet = "Sheet1")
foo1 = as.data.frame(foo[,1:3])
foo2 = -log10(data.matrix(foo[,4:length(foo)]))
rownames(foo2) = as.character(seq(1, 275, by = 1))
colnames(foo2) = as.character(seq(350, 2500, by = 1))
#pre-process the spectra: resample it to a resolution of 2 nm 
old_wavs = foo2 %>% colnames() %>% as.numeric()
new_wavs = seq(502, 2500, by = 2)
foo1$spc = foo2 %>% 
  resample(wav = old_wavs, new.wav = new_wavs)
Ber = foo1
rm(foo,foo1,foo2)

#Belgium Data
foo = read_csv("Belgium.csv")
foo1 = as.data.frame(foo[,701:703])
foo2 = data.matrix(foo[,1:700])
rownames(foo2) = as.character(seq(1, 825, by = 1))
colnames(foo2) = as.character(seq(1100, 2498, by = 2))
foo1$SOC[foo1$SOC == -9999] <- NA
foo1$CEC[foo1$CEC == -9999] <- NA
foo1$N[foo1$N == -9999] <- NA
foo1$spc = foo2 
Belgium = foo1
rm(foo,foo1,foo2)

#Data provision
IP = OSSL
IP = kar
IP = Ber
IP = Belgium
#Component selection
property = "SOC"
summary(IP[,property])

imp_no_comp = list()
imp_no_samples = list()
IP_X1 = IP$spc[!is.na(IP[,property]),]
IP_Y1 = IP[,property][!is.na(IP[,property])]
summary(IP_Y1)
IP_X_O = IP_X1[order(IP_Y1), ]
IP_Y_O = IP_Y1[order(IP_Y1)]
#For clay, silt, sand the values were made to selected in range of 0-100
#IP_Y1[IP_Y1 < 0 | IP_Y1 > 100] <- NA
#For SOC the values were made to selected >0.1
# IP_Y1[IP_Y1 < 0.1] <- NA
# IP_Y = IP_Y1[!is.na(IP_Y1)]
# IP_X = IP_X1[!is.na(IP_Y1),]
# summary(IP_Y)
# IP_X_O = IP_X[order(IP_Y), ]
# IP_Y_O = IP_Y[order(IP_Y)]
#For Avail_P the values were made to selected in range of 0-200
#IP_Y1[IP_Y1 < 0.5 | IP_Y1 > 200] <- NA
#IP_Y = IP_Y1[!is.na(IP_Y1)]
#IP_X = IP_X1[!is.na(IP_Y1),]
#summary(IP_Y)
#IP_X_O = IP_X[order(IP_Y), ]
#IP_Y_O = IP_Y[order(IP_Y)]
#Split 7:3
X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), ]
Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)]
X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE), ]
Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)]
rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 5:1
# X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 2:1
# X1 = IP_X_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 1:1
# X1 = IP_X_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

# #For Karnataka, Berambadi, Belgium
# #Split 5:1
# X1 = IP_X_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 2:1
# X1 = IP_X_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)
# #Split 1:1
# X1 = IP_X_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), ]
# Y1 = IP_Y_O[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)]
# X2 = IP_X_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), ]
# Y2 = IP_Y_O[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)]
# rm(IP_X,IP_X1,IP_Y,IP_Y1,IP_X_O,IP_Y_O)

####1NEAREST NEIGHBOUR####
##Choosing Dissimilarity measures##
##Training##
foo = X1
foo1 = Y1

#Optional
optimal_sel =  list(method = "opc", value = min(length(foo1)-4,40))
#optimal_sel =  list(method = "manual", value = 27)
#no_components = 27
pls_tr_opc = ortho_projection(Xr = foo, Yr = foo1, method = "pls", pc_selection = optimal_sel)
pls_tr_opc
pls_comp_num = as.matrix(as.numeric(pls_tr_opc$opc_evaluation[,1]))
pls_comp_value = as.matrix(as.numeric(pls_tr_opc$opc_evaluation[,2]))
old_par = par("mfrow")
par(mfrow = c(1,1))
matplot(x = pls_comp_num, y = pls_comp_value, 
        xlab = "No of PLS components",
        ylab = "RMSD of Yr",
        type = "h", lty = 1, col = "#FF1A00CC")
title("method=pls.opc")
#Finding which k has the minimum rmse
no_components = which.min(as.numeric(pls_tr_opc$opc_evaluation[,2]))
print(paste("Number of components:", no_components))
imp_no_comp = append(imp_no_comp, no_components)

#Comparing all the dissimilarity measures#
#PC dissimilarity with default settings (variance-based no. of components)
pcad = dissimilarity(foo, diss_method = "pca", scale = TRUE)
#PLS dissimilarity with default settings (variance-based no of components)
plsd = dissimilarity(foo, diss_method = "pls", Yr = foo1,scale = TRUE)
#PC dissimilarity with optimal selection of components
opc_sel =  list(method = "opc", value = min(length(foo1)-4,40))
o_pcad = dissimilarity(foo,diss_method = "pca",Yr = foo1,pc_selection = opc_sel, scale = TRUE)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(foo,diss_method = "pls",Yr = foo1,pc_selection = opc_sel, scale = TRUE)
#Correlation dissimilarity 
cd = dissimilarity(foo, diss_method = "cor", scale = TRUE)
#Moving window correlation dissimilarity 
mcd = dissimilarity(foo, diss_method = "cor", ws = 51, scale = TRUE)
#Euclidean dissimilarity 
ed = dissimilarity(foo, diss_method = "euclid", scale = TRUE)
#Cosine dissimilarity 
cosd = dissimilarity(foo, diss_method = "cosine", scale = TRUE)
#Spectral information divergence/dissimilarity 
sinfd = dissimilarity(foo, diss_method = "sid", scale = TRUE)

#Evaluations
Y_matrix = as.matrix(foo1)
ev = NULL
ev[["pcad"]] = sim_eval(pcad$dissimilarity, side_info = Y_matrix)
ev[["plsd"]] = sim_eval(plsd$dissimilarity, side_info = Y_matrix)
ev[["o_pcad"]] = sim_eval(o_pcad$dissimilarity, side_info = Y_matrix)
ev[["o_plsd"]] = sim_eval(o_plsd$dissimilarity, side_info = Y_matrix)
ev[["cd"]] = sim_eval(cd$dissimilarity, side_info = Y_matrix)
ev[["mcd"]] = sim_eval(mcd$dissimilarity, side_info = Y_matrix)
ev[["ed"]] = sim_eval(ed$dissimilarity, side_info = Y_matrix)
ev[["cosd"]] = sim_eval(cosd$dissimilarity, side_info = Y_matrix)
ev[["sinfd"]] = sim_eval(sinfd$dissimilarity, side_info = Y_matrix)

#Tabulating and plotting
statistics_result = NULL
for (label in names(ev)) {
  observed_values = ev[[label]]$first_nn[, 1]
  predicted_values = ev[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_calib = bind_rows(statistics_result, .id = "Methods")
print(r_calib)
write.csv(r_calib, "r_calib.csv", row.names = FALSE)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
par(mfrow = c(3, 3))
p = sapply(names(ev), 
           FUN = function(x, label, labs = c("Clay, %", "Clay (1-NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev))])
             title(label)
             grid()
             abline(0, 1)
             text(
               (max(xy[,1])-10), (max(xy[,2])-10),
               labels = paste("RMSD:", round(x[[label]]$eval[1],3), "\nR:", round(x[[label]]$eval[2],3), "\nR2:", round(x[[label]]$eval[2]^2,3)),
               pos = 1,
               col = "black", cex = 1
             )
           },
           x = ev)
par(mfrow = c(1, 1))

##Testing##
#PC dissimilarity with default settings (variance-based no. of components)
pcad = dissimilarity(Xr=foo, Xu=X2, diss_method = "pca", scale = TRUE)
#PLS dissimilarity with default settings (variance-based no of components)
plsd = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls", scale = TRUE)
#PC dissimilarity with optimal selection of components
opc_sel = list("opc", 40)
o_pcad = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pca",pc_selection = opc_sel, scale = TRUE)
#PLS dissimilarity with optimal selection of components
o_plsd = dissimilarity(Xr=foo, Yr = foo1, Xu = X2, Yu = Y2, diss_method = "pls",pc_selection = opc_sel, scale = TRUE)
#Correlation dissimilarity 
cd = dissimilarity(Xr=foo, Xu=X2, diss_method = "cor", scale = TRUE)
#Moving window correlation dissimilarity 
mcd = dissimilarity(Xr=foo, Xu=X2, diss_method = "cor", ws = 51, scale = TRUE)
#Euclidean dissimilarity 
ed = dissimilarity(Xr=foo, Xu=X2, diss_method = "euclid", scale = TRUE)
#Cosine dissimilarity 
cosd = dissimilarity(Xr=foo, Xu=X2, diss_method = "cosine", scale = TRUE)
#Spectral information divergence/dissimilarity 
sinfd = dissimilarity(Xr=foo, Xu=X2, diss_method = "sid", scale = TRUE)

#Evaluations
ev = NULL
observed_values = Y2
ev[["pcad"]]$first_nn = cbind(observed_values, foo1[apply(pcad$dissimilarity, 2, which.min)])
ev[["plsd"]]$first_nn = cbind(observed_values, foo1[apply(plsd$dissimilarity, 2, which.min)])
ev[["o_pcad"]]$first_nn = cbind(observed_values, foo1[apply(o_pcad$dissimilarity, 2, which.min)])
ev[["o_plsd"]]$first_nn = cbind(observed_values, foo1[apply(o_plsd$dissimilarity, 2, which.min)])
ev[["cd"]]$first_nn = cbind(observed_values, foo1[apply(cd$dissimilarity, 2, which.min)])
ev[["mcd"]]$first_nn = cbind(observed_values, foo1[apply(mcd$dissimilarity, 2, which.min)])
ev[["ed"]]$first_nn = cbind(observed_values, foo1[apply(ed$dissimilarity, 2, which.min)])
ev[["cosd"]]$first_nn = cbind(observed_values, foo1[apply(cosd$dissimilarity, 2, which.min)])
ev[["sinfd"]]$first_nn = cbind(observed_values, foo1[apply(sinfd$dissimilarity, 2, which.min)])

statistics_result = NULL
for (label in names(ev)) {
  observed_values = ev[[label]]$first_nn[, 1]
  predicted_values = ev[[label]]$first_nn[, 2]
  statistics_result[[label]] = calculate_statistics(observed_values, predicted_values)
}
r_calib = bind_rows(statistics_result, .id = "Methods")
print(r_calib)
write.csv(r_calib, "r_valid.csv", row.names = FALSE)

colours1 = c("red", "blue","green","yellow","pink","violet","orange","magenta","cyan")
par(mfrow = c(3, 3))
p = sapply(names(ev), 
           FUN = function(x, label, labs = c("Clay, %", "Clay (1-NN), %")) {
             xy = x[[label]]$first_nn[,1:2]
             plot(xy[,1], xy[,2], xlab = labs[1], ylab = labs[2], col = colours1[match(label,names(ev))])
             title(label)
             grid()
             abline(0, 1)
             text(
               (max(xy[,1])-10), (max(xy[,2])-10),
               labels = paste("RMSD:", round(statistics_result[[label]]$RMSE,3), "\nR2:", round(statistics_result[[label]]$R_squared,3)),
               pos = 1,
               col = "black", cex = 1
             )
           },
           x = ev)
par(mfrow = c(1, 1))
rm(foo,foo1)

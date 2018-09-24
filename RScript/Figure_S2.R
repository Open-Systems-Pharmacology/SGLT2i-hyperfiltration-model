#This script simulates the QSP Diabetes model and generates Supplementary Figure S2.

#The below loaded libraries must be installed before running this script.
library(MoBiToolboxForR)
library(parallel)

#Define here the path to the folder with the experimental data extracted from literature.
#For the description of the datasets, refer to the main text.
expDataFolder = "..\\ExpData\\";
#Define here the path to the folder with the simulation files.
simFolder = "..\\Models\\xml\\";
#Define here the path to the folder where the output figure will be stored.
figureFolder = "..\\Figures\\";
decSym = '.';

#Simulate the single ascending dose study as reported by Komoroski et al. 2009
sim_SAD = function(){
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, "simFolder")
clusterEvalQ(cl, library(MoBiToolboxForR))

modelNames = c("Placebo", "Dapagliflozin_2.5mg", "Dapagliflozin_5mg", "Dapagliflozin_10mg", "Dapagliflozin_20mg", 
               "Dapagliflozin_50mg", "Dapagliflozin_100mg", "Dapagliflozin_250mg", "Dapagliflozin_500mg");

executeSim = function(simName){
  initStruct <- list();
  currDCI = initSimulation(XML = paste0(simFolder, simName, ".xml"), ParamList = initStruct, whichInitParam = "none");
  currDCI = processSimulation(DCI_Info = currDCI);
  results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
  return(results_UGE[length(results_UGE)])
}

results = parSapply(cl, modelNames, executeSim);
stopCluster(cl)
return(results);
}

#Plot UGE as a functino of dapagliflozin dose.
plot_UGE_vs_Dapa = function(results){
  #The first value, 2.5/2, is a dummy for placebo - no dose is applied, but the value is used in the axis capture.
  dapaDose = c((2.5/2), 2.5, 5, 10, 20, 50, 100, 250, 500)
  #UGE data from Komoroski et al. 2009
cumulativeGlucoseExcretionData = c(2.25, 33.44, 32.22, 78.60, 65.35, 118.49, 196.5, 226, 303.24);

plot(dapaDose, cumulativeGlucoseExcretionData, ylim = c(0, 350), log="x", xlab = "Dapagliflozin dose [mg], log", ylab = "Cumulative UGE over 120 h [g]",
     lwd=1, pch = 4,
     xaxt = "n");
points(dapaDose, results*180e-6, pch =1, type ="b", lwd=1);
axis(1, at = dapaDose, labels = c("Placebo", "2.5", "5", "10", "20", "50", "100", "250", "500"));
legend("topleft", c("Measured", "Simulated"),
       lty=c(NA, 1),
       pch=c(4,21),
       bty="n", y.intersp = 1.1)
}

ltyArr = c(3, 4, 1, 2)
plot_QSP = function(){
  #Name of the model xml-file.
  modelName = "Dapa_14d";
  #The list "initStruct" stores the paths of the parameters to be initialized (and changed later on).
  initStruct <- list();
  initStruct <- initParameter(path_id = "*|Applications|100mg PO_14d|Dapaglifllozin_Capsule|Dose", initializeIfFormula = "always", 
                              initStruct = initStruct)
  #Initialize the simulation with the parameters defined in "initStruct".
  currDCI = initSimulation(XML = paste0(simFolder, modelName, ".xml"), ParamList = initStruct);

  #Simulate without dapagliflozin
  currDCI = setParameter(0, path_id = "*|Applications|100mg PO_14d|Dapaglifllozin_Capsule|Dose", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the palcebo results.
  results_daytime_plac = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_plac = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_plac = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 10mg dapagliflozin
  currDCI = setParameter(10e-6, path_id = "*|Applications|100mg PO_14d|Dapaglifllozin_Capsule|Dose", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_dapa = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_dapa = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_dapa = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  
  par(mar=c(4,4,0.4,4))
  idx = 1;
  plot(results_GFR_plac[, 1][], results_GFR_plac[, 2][]*1e3,
       log="", type="l",
       axes = FALSE,
       xlab = "",
       ylab = "",
       xlim = c(13*24*60-3*60, 14*24*60-3*60),
       ylim=c(80, 140), lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_GFR_dapa[, 1][], results_GFR_dapa[, 2][]*1e3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  axis(2)
  mtext("GFR [mL/min]", side = 2, line = 2.5)
  
  par(new=TRUE);
  idx = idx + 1;
  plot(results_G_plac[, 1][], results_G_plac[, 2][]*1e-3,
       log="", type="l",
       axes = FALSE,
       xlab = "",
       ylab = "",
       xlim = c(13*24*60-3*60, 14*24*60-3*60),
       ylim=c(3, 10), lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_G_dapa[, 1][], results_G_dapa[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  axis(4)
  mtext("Glucose concentration [mmol/L]", side = 4, line = 2.5)
  axis(1, at = (c(6, 9, 12, 15, 18, 21, 24)-9)*60+(13*24*60), labels = c(6, 9, 12, 15, 18, 21, 24))
  mtext("Time [h]", side = 1, line = 3)
  
  legend("topleft", c("GFR placebo", "GFR 10 mg dapa"),
         lty=c(3,4),
         bty="n", y.intersp = 1.1)
  
  legend("bottomleft", c("Glucose placebo", "Glucose 10 mg dapa"),
         lty=c(1,2),
         bty="n", y.intersp = 1.1)
}

results_SAD = sim_SAD();

# png(file=paste0(figureFolder, "Figure_S2.png"),
#     width=16,
#     height=8,
#     units = "cm",
#     res=300,
#     pointsize=8)
pdf(file=paste0(figureFolder, "Figure_S2.pdf"),
    width=6.3,
    height=3.15,
    pointsize=8)

par(mfrow = c(1, 2), cex=1, oma=c(0,0,0,0))
par(mar=c(4,4,0.4,0.1))

plot_UGE_vs_Dapa(results_SAD)
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "A",cex=1,font=1.5, xpd=T)
plot_QSP();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "B",cex=1,font=1.5, xpd=T)
dev.off();

#This script generates Figure 4.

#The below loaded libraries must be installed before running this script.
library(MoBiToolboxForR)
library(parallel)

#Define here the path to the folder with the experimental data extracted from literature.
#For the description of the datasets, refer to the main text.
expDataFolder = "..\\ExpData\\";
#Define here the path to the folder with the simulation files.
simFolder = "..\\xml\\";
#Define here the path to the folder where the output figure will be stored.
figureFolder = "..\\Figures\\";
decSym = '.';

#These are the applied glucose concentrations.
xVals_target =  c(5.6, 8.3, 11.1, 13.9, 16.7, 19.4, 22.2, 25, 27.8, 30.5)
#These are the time points at which glomerular filtration rate (GFR) and urinary glucose excretion (UGE) are
#read out.
xTime = c(40, 80, 120, 160, 200, 240, 280, 320, 360, 400);
#Time at which SHC simulation without dapagliflozin administration starts.
xTime_offset = 240;

#Name of the model xml-file.
modelName = "Dapagliflozin_2013_SHC_healthy";

ltyArr = c(1, 2, 3, 1, 1)

#Plot UGE results
plotUGE = function(results, idx){
  mw = 180;
  convFac = 1e-6*mw;
  yVals = c();
  for (time in xTime){
    yVals = c(yVals, results[time + xTime_offset, 2]);
  }
  yAxtMarks = c(0, 5, 10, 15, 20, 25, 30);
  
  plot(xVals_target, yVals*convFac,
       log="", type="l", xlab="", ylab="",
       ylim=c(0, 35), yaxt="n", xaxt = "n", lwd=1, lty=ltyArr[idx])
  
  if (idx == 4){
    points(xVals_target, yVals*convFac,
    type="p", lwd=1, pch = 4);
  }
  if (idx == 5){
    points(xVals_target, yVals*convFac,
    type="p", lwd=1, pch = 21);
  }
  if (idx == 1){
    axis(2, at = yAxtMarks, labels = yAxtMarks)
    axis(1, at = yAxtMarks, labels = yAxtMarks)
    title(xlab="Plasma glucose [mmol/L]", ylab="UGE over 40 minutes [g]") 
  }
}

#Plot GFR results
plotGFR = function(results, idx){
  plot(0:440, results[, 2][0:440+(xTime_offset-40)]*1e3,
       log="", type="l", xlab="", ylab="", main = "",
       ylim=c(80, 160), xaxt="n", yaxt = "n", lwd=1, lty=ltyArr[idx]);

  if (idx == 4){
    points(seq(0, 400, 40)+20, results[, 2][seq(0, 400, 40)+20+(xTime_offset-40)]*1e3,
           type="p", lwd=1, pch = 4);
  }
  if (idx == 5){
    points(seq(0, 400, 40)+20, results[, 2][seq(0, 400, 40)+20+(xTime_offset-40)]*1e3,
           type="p", lwd=1, pch = 21);
  }
  if (idx == 1){
    ylabs = c(80, 100, 120, 140, 160)
    title(xlab = "Plasma glucose [mmol/L]", ylab = "GFR [mL/min]")
    axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
    axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
    axis(2, at = ylabs, labels = ylabs);
  }
}

sim_SHC = function(){
  #The list "initStruct" stores the paths of the parameters to be initialized (and changed later on).
  initStruct <- list();
  #Initialize the simulation with the all parameters that are not defined by a formula.
  currDCI = initSimulation(XML = paste0(simFolder, modelName, ".xml"),  whichInitParam = "allNonFormula")#ParamList = initStruct);
  def_SGLT2 = getParameter(path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI)$Value;
  def_SGLT1 = getParameter(path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI)$Value;
  
    #Set simulation time.
  currDCI = setSimulationTime(timepoints = {0 : (xTime_offset + 450)}, DCI_Info = currDCI);
  
  idx = 1;
    #Run the simulation.
    currDCI = processSimulation(DCI_Info = currDCI);
    #Get the GFR results.
    results_GFR = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
    #Get the UGE results.
    results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
    
    par(mfg = c(1, 1));
    plotUGE(results_UGE, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "A",cex=1,font=1.5, xpd=T)
    
    par(mfg = c(1, 2));
    plotGFR(results_GFR, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "B",cex=1,font=1.5, xpd=T)
    
    sglt2Vals = c(0.5, 1);
    sglt1vals = c(1);
    for (parVal in sglt2Vals){
      currDCI = setParameter(1 - parVal, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
      #Run the simulation.
      currDCI = processSimulation(DCI_Info = currDCI);
      #Get the GFR results.
      results_GFR = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
      #Get the UGE results.
      results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
      idx = idx + 1;
      par(mfg = c(1, 1));
      plotUGE(results_UGE, idx);
      par(mfg = c(1, 2));
      plotGFR(results_GFR, idx);
    }
    
    currDCI = setParameter(1, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
    for (parVal in sglt1vals){
      currDCI = setParameter(def_SGLT1 * (1 - parVal), path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI);
      #Run the simulation.
      currDCI = processSimulation(DCI_Info = currDCI);
      #Get the GFR results.
      results_GFR = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
      #Get the UGE results.
      results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
      idx = idx + 1;
      par(mfg = c(1, 1));
      plotUGE(results_UGE, idx);
      par(mfg = c(1, 2));
      plotGFR(results_GFR, idx);
    }
    
    currDCI = setParameter(0, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
      currDCI = setParameter(0, path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI);
      #Run the simulation.
      currDCI = processSimulation(DCI_Info = currDCI);
      #Get the GFR results.
      results_GFR = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
      #Get the UGE results.
      results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
      idx = idx + 1;
      par(mfg = c(1, 1));
      plotUGE(results_UGE, idx);
      par(mfg = c(1, 2));
      plotGFR(results_GFR, idx);
    
    par(mfg = c(1, 1));
    legend("topleft", c("SGLT2 50% inhibition", "SGLT2 100% inhibition", "SGLT1 100% inhibition", "SGLT1 and 2 100% inhibition", "Placebo"),
           lty=c(2, 3, 1, 1, 1),
           pch=c(NA, NA, 4, 21, NA),
           bty="n", y.intersp = 1.1)
}

plot_QSP_healthy = function(){
  #Name of the model xml-file.
  modelName = "Dapa_14d";
  #The list "initStruct" stores the paths of the parameters to be initialized (and changed later on).
  initStruct <- list();
  initStruct <- initParameter(path_id = "*|Applications|100mg PO_14d|Dapaglifllozin_Capsule|Dose", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", initializeIfFormula = "always", 
                              initStruct = initStruct)
  #Initialize the simulation with the parameters defined in "initStruct".
  currDCI = initSimulation(XML = paste0(simFolder, modelName, ".xml"), ParamList = initStruct);
  currDCI = setSimulationTime(timepoints = {0 : (14*24*60)}, DCI_Info = currDCI);
  def_SGLT2 = getParameter(path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI)$Value;
  def_SGLT1 = getParameter(path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI)$Value;
  
  #Simulate without inhibition
  currDCI = setParameter(0, path_id = "*|Applications|100mg PO_14d|Dapaglifllozin_Capsule|Dose", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the palcebo results.
  results_daytime_plac = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_plac = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_plac = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 50% SGLT2 inhibition
  currDCI = setParameter(0.5, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT2_05 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT2_05 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT2_05 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 100% SGLT2 inhibition
  currDCI = setParameter(0, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT2_0 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT2_0 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT2_0 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 100% SGLT1 inhibition
  currDCI = setParameter(1, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  currDCI = setParameter(def_SGLT1 * 0, path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT1_0 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT1_0 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT1_0 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 100% SGLT1 and 2 inhibition
  currDCI = setParameter(0, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  currDCI = setParameter(def_SGLT1 * 0, path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT12_0 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT12_0 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT12_0 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);

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
  points(results_GFR_SGLT2_05[, 1][], results_GFR_SGLT2_05[, 2][]*1e3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_GFR_SGLT2_0[, 1][], results_GFR_SGLT2_0[, 2][]*1e3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_GFR_SGLT1_0[, 1][], results_GFR_SGLT1_0[, 2][]*1e3,
         type="l",
         lwd=1, lty=1);
  points(results_GFR_SGLT1_0[, 1][seq(0, length(results_GFR_SGLT1_0[, 1][]), 80)], 
         results_GFR_SGLT1_0[, 2][seq(0, length(results_GFR_SGLT1_0[, 1][]), 80)]*1e3,
         type="p",
         lwd=1, lty=1, pch = 4);
  
  idx = idx + 1;
  points(results_GFR_SGLT12_0[, 1][], results_GFR_SGLT12_0[, 2][]*1e3,
         type="l",
         lwd=1, lty=1);
  points(results_GFR_SGLT12_0[, 1][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)],
         results_GFR_SGLT12_0[, 2][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)]*1e3,
         type="p",
         lwd=1, lty=1, pch = 21);
  
  axis(2)
  mtext("GFR [mL/min]", side = 2, line = 2.5)
  
  par(new=TRUE);
  idx = 1;
  plot(results_G_plac[, 1][], results_G_plac[, 2][]*1e-3,
       log="", type="l",
       axes = FALSE,
       xlab = "",
       ylab = "",
       xlim = c(13*24*60-3*60, 14*24*60-3*60),
       ylim=c(3, 20), lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_G_SGLT2_05[, 1][], results_G_SGLT2_05[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_G_SGLT2_0[, 1][], results_G_SGLT2_0[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_G_SGLT1_0[, 1][], results_G_SGLT1_0[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  points(results_G_SGLT1_0[, 1][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)], 
         results_G_SGLT1_0[, 2][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)]*1e-3,
         type="p",
         lwd=1, lty=1, pch = 4);
  
  idx = idx + 1;
  points(results_G_SGLT12_0[, 1][], results_G_SGLT12_0[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  points(results_G_SGLT12_0[, 1][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)],
         results_G_SGLT12_0[, 2][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)]*1e-3,
         type="p",
         lwd=1, lty=1, pch = 21);
  
  axis(4, at = c(3, 5, 10, 15, 20))
  mtext("Glucose concentration [mmol/L]", side = 4, line = 2.5)
  axis(1, at = (c(6, 9, 12, 15, 18, 21, 24)-9)*60+(13*24*60), labels = c(6, 9, 12, 15, 18, 21, 24))
  mtext("Time [h]", side = 1, line = 3)
}

plot_QSP_T2DM = function(){
  #Name of the model xml-file.
  modelName = "Dapa_14d";
  #The list "initStruct" stores the paths of the parameters to be initialized (and changed later on).
  initStruct <- list();
  initStruct <- initParameter(path_id = "*|Applications|100mg PO_14d|Dapaglifllozin_Capsule|Dose", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|S_I", initializeIfFormula = "always", 
                              initStruct = initStruct)
  #Initialize the simulation with the parameters defined in "initStruct".
  currDCI = initSimulation(XML = paste0(simFolder, modelName, ".xml"), ParamList = initStruct);
  currDCI = setSimulationTime(timepoints = {0 : (14*24*60)}, DCI_Info = currDCI);
  def_SGLT2 = getParameter(path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI)$Value;
  def_SGLT1 = getParameter(path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI)$Value;
  
  #Simulate without inhibition
  currDCI = setParameter(0, path_id = "*|Applications|100mg PO_14d|Dapaglifllozin_Capsule|Dose", DCI_Info = currDCI);
  currDCI = setParameter(0.12, path_id = "*|Organism|S_I", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the palcebo results.
  results_daytime_plac = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_plac = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_plac = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 50% SGLT2 inhibition
  currDCI = setParameter(0.5, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT2_05 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT2_05 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT2_05 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 100% SGLT2 inhibition
  currDCI = setParameter(0, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT2_0 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT2_0 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT2_0 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 100% SGLT1 inhibition
  currDCI = setParameter(1, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  currDCI = setParameter(def_SGLT1 * 0, path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT1_0 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT1_0 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT1_0 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
  #Simulate 100% SGLT1 and 2 inhibition
  currDCI = setParameter(0, path_id = "*|Organism|Kidney|Cortex|PCT|SGLT2|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  currDCI = setParameter(def_SGLT1 * 0, path_id = "*|Organism|Kidney|Medulla|PST|SGLT1|Relative expression \\(normalized\\)", DCI_Info = currDCI);
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the dapa results.
  results_daytime_SGLT12_0 = getSimulationResult(path_id = "*|Events|Daytime|Daytime", DCI_Info = currDCI);
  results_GFR_SGLT12_0 = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
  results_G_SGLT12_0 = getSimulationResult(path_id = "*|Organism|PeripheralVenousBlood|Glucose|Plasma", DCI_Info = currDCI);
  
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
  points(results_GFR_SGLT2_05[, 1][], results_GFR_SGLT2_05[, 2][]*1e3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_GFR_SGLT2_0[, 1][], results_GFR_SGLT2_0[, 2][]*1e3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_GFR_SGLT1_0[, 1][], results_GFR_SGLT1_0[, 2][]*1e3,
         type="l",
         lwd=1, lty=1);
  points(results_GFR_SGLT1_0[, 1][seq(0, length(results_GFR_SGLT1_0[, 1][]), 80)], 
         results_GFR_SGLT1_0[, 2][seq(0, length(results_GFR_SGLT1_0[, 1][]), 80)]*1e3,
         type="p",
         lwd=1, lty=1, pch = 4);
  
  idx = idx + 1;
  points(results_GFR_SGLT12_0[, 1][], results_GFR_SGLT12_0[, 2][]*1e3,
         type="l",
         lwd=1, lty=1);
  points(results_GFR_SGLT12_0[, 1][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)],
         results_GFR_SGLT12_0[, 2][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)]*1e3,
         type="p",
         lwd=1, lty=1, pch = 21);
  
  axis(2)
  mtext("GFR [mL/min]", side = 2, line = 2.5)
  
  par(new=TRUE);
  idx = 1;
  plot(results_G_plac[, 1][], results_G_plac[, 2][]*1e-3,
       log="", type="l",
       axes = FALSE,
       xlab = "",
       ylab = "",
       xlim = c(13*24*60-3*60, 14*24*60-3*60),
       ylim=c(3, 20), lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_G_SGLT2_05[, 1][], results_G_SGLT2_05[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_G_SGLT2_0[, 1][], results_G_SGLT2_0[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(results_G_SGLT1_0[, 1][], results_G_SGLT1_0[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  points(results_G_SGLT1_0[, 1][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)], 
         results_G_SGLT1_0[, 2][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)]*1e-3,
         type="p",
         lwd=1, lty=1, pch = 4);
  
  idx = idx + 1;
  points(results_G_SGLT12_0[, 1][], results_G_SGLT12_0[, 2][]*1e-3,
         type="l",
         lwd=1, lty=ltyArr[idx]);
  points(results_G_SGLT12_0[, 1][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)],
         results_G_SGLT12_0[, 2][seq(0, length(results_GFR_SGLT12_0[, 1][]), 80)]*1e-3,
         type="p",
         lwd=1, lty=1, pch = 21);
  
  axis(4, at = c(3, 5, 10, 15, 20))
  mtext("Glucose concentration [mmol/L]", side = 4, line = 2.5)
  axis(1, at = (c(6, 9, 12, 15, 18, 21, 24)-9)*60+(13*24*60), labels = c(6, 9, 12, 15, 18, 21, 24))
  mtext("Time [h]", side = 1, line = 3)
}

# png(file=paste0(figureFolder, "Figure_4.png"),
#     width=17.3,
#     height=11.46,
#     units = "cm",
#     res=300,
#     pointsize=8)

pdf(file=paste0(figureFolder, "Figure_4.pdf"),
    width=6.8,
    height=4.5,
    pointsize=8)

par(mfrow = c(2, 2), cex=1, oma=c(0,0,0,0))
par(mar=c(4,4,0.4,0.1))

sim_SHC()

par(mfg = c(2, 1));
plot_QSP_healthy();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "C",cex=1,font=1.5, xpd=T)
par(mfg = c(2, 2));
plot_QSP_T2DM();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "D",cex=1,font=1.5, xpd=T)
dev.off();

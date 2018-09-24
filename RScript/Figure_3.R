#This script simulates the stepped hyperglycemic clamp (SHC) with variation of selected parameters
#and generates Figure 3.

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

#Read experimental data
dataDeFronzo_healthy_placebo = read.table(paste0(expDataFolder, "DeFronzo_2013_SHC_Healthy_placebo.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataDeFronzo_healthy_dapa = read.table(paste0(expDataFolder, "DeFronzo_2013_SHC_Healthy_dapa.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)

#These are the applied glucose concentrations, reported in Lu et al.2014
xVals_placebo = c(5.8, 10.3, 12, 15.1, 17.2, 19, 22.9, 25.2, 28.5, 30);
xVals_dapa =    c(5.5, 9, 12.5, 14.6, 17.1, 19.6, 22.1, 24.1, 28.5, 30.4);
xVals_target =  c(5.6, 8.3, 11.1, 13.9, 16.7, 19.4, 22.2, 25, 27.8, 30.5)
#These are the time points at which glomerular filtration rate (GFR) and urinary glucose excretion (UGE) are
#read out.
xTime = c(40, 80, 120, 160, 200, 240, 280, 320, 360, 400);
#Time at which SHC simulation without dapagliflozin administration starts.
xTime_offset_placebo = 240;
#Time at which SHC simulation with dapagliflozin administration starts.
xTime_offset_dapa = 10140;

#Name of the model xml-file.
modelName = "Dapagliflozin_2013_SHC_healthy";

pchArr = c(1, 4, 5)
ltyArr = c(1, 2, 3)

#Plot GFR results without dapagliflozin adminsitration
plotGFR_placebo = function(results, idx){
  plot(0:440, results[, 2][0:440+(xTime_offset_placebo-40)]*1e3,
       log="", type="l", xlab="", ylab="", main = "",
       ylim=c(80, 160), xaxt="n", yaxt = "n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx]);

  ylabs = c(80, 120, 160)
  if (idx == 1){
    axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
    axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
    axis(2, at = ylabs, labels = ylabs);

    title(ylab="GFR [mL/min]") 
  }
}

#Plot GFR results with dapagliflozin adminsitration
plotGFR_dapa = function(results, idx){
  plot(0:440, results[, 2][0:440+(xTime_offset_dapa-40)]*1e3,
       log="", type="l", xlab="", ylab="", main = "",
       ylim=c(80, 160), xaxt="n", yaxt = "n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx]);
  
  ylabs = c(80, 120, 160)
  if (idx == 1){
    axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
    axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
    axis(2, at = ylabs, labels = ylabs);
  }
}

#Plot UGE results without dapagliflozin adminsitration
plotUGE_placebo = function(results, idx){
  mw = 180;
  convFac = 1e-6*mw;
  yVals = c();
  for (time in xTime){
    yVals = c(yVals, results[time + xTime_offset_placebo, 2]);
  }
  yAxtMarks = c(0, 5, 10, 15, 20, 25, 30);
  
  plot(xVals_placebo, yVals*convFac,
       log="", type="l", xlab="", ylab="",
       ylim=c(0, 35), yaxt="n", xaxt = "n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  
  if (idx == 1){
    points(xVals_placebo, dataDeFronzo_healthy_placebo$Amount*convFac, pch = pchArr[idx], lwd=1, type="p", lty=ltyArr[idx])
    arrows(xVals_placebo, dataDeFronzo_healthy_placebo$Amount*convFac - dataDeFronzo_healthy_placebo$Error*convFac,
           xVals_placebo, dataDeFronzo_healthy_placebo$Amount*convFac + dataDeFronzo_healthy_placebo$Error*convFac,
           length=0.05, angle=90, code=3, lwd=1)
    axis(2, at = yAxtMarks, labels = yAxtMarks)
    axis(1, at = yAxtMarks, labels = yAxtMarks)
    
    title(ylab="UGE [g]") 
  }
}

#Plot UGE results with dapagliflozin adminsitration
plotUGE_dapa = function(results, idx){
  mw = 180;
  convFac = 1e-6*mw;
  yVals = c();
  for (time in xTime){
    yVals = c(yVals, results[time + xTime_offset_dapa, 2]);
  }
  yAxtMarks = c(0, 5, 10, 15, 20, 25, 30);
  
  plot(xVals_dapa, yVals*convFac,
       log="", type="l", xlab="", ylab="",
       ylim=c(0, 35), yaxt="n", xaxt = "n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  if (idx == 1){
    points(xVals_dapa, dataDeFronzo_healthy_dapa$Amount*convFac, pch = pchArr[idx], lwd=1, type="p", lty=ltyArr[idx])
    arrows(xVals_dapa, dataDeFronzo_healthy_dapa$Amount*convFac - dataDeFronzo_healthy_dapa$Error*convFac,
           xVals_dapa, dataDeFronzo_healthy_dapa$Amount*convFac + dataDeFronzo_healthy_dapa$Error*convFac,
           length=0.05, angle=90, code=3, lwd=1)
    axis(2, at = yAxtMarks, labels = yAxtMarks)
    axis(1, at = yAxtMarks, labels = yAxtMarks)
  }
}

#Variation of SLGT2 expression
sim_SGLT2 = function(){
  #The list "initStruct" stores the paths of the parameters to be initialized (and changed later on).
  initStruct <- list();
  #Initialize the simulation with the all parameters that are not defined by a formula.
  currDCI = initSimulation(XML = paste0(simFolder, modelName, ".xml"), ParamList = initStruct, whichInitParam = "allNonFormula");
  #Set simulation time.
  currDCI = setSimulationTime(timepoints = {0:(xTime_offset_dapa+450)}, DCI_Info = currDCI);
  
  #Path to the parameter that will be variated
  parPath = "*|SGLT2|Reference concentration";
  #Values of the parameter for simulations.
  parVals = c(1, 0.8, 1.2);
  
  idx = 1;
  #Run simulations with variated parameter value.
  for (parVal in parVals){
    #Set the value of the parameter.
    currDCI = setParameter(value = parVal, path_id = parPath, DCI_Info = currDCI);
    #Run the simulation.
    currDCI = processSimulation(DCI_Info = currDCI);
    #Get the GFR results.
    results_GFR = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
    #Get the mean GFR during SHC without dapagliflozin.
    results_GFR_mean_placebo = getSimulationResult(path_id = "*|Organism|Kidney|DEBUG_GFR_mean_placebo", DCI_Info = currDCI);
    #Get the mean GFR during SHC with dapagliflozin.
    results_GFR_mean_dapa = getSimulationResult(path_id = "*|Organism|Kidney|DEBUG_GFR_mean_dapa", DCI_Info = currDCI);
    #Get the UGE results.
    results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
    
    par(mfg = c(3, 1));
    plotUGE_placebo(results_UGE, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "E",cex=1,font=1.5, xpd=T)
    
    par(mfg = c(3, 2));
    plotUGE_dapa(results_UGE, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "F",cex=1,font=1.5, xpd=T)
    
    par(mfg = c(4, 1));
    plotGFR_placebo(results_GFR, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "G",cex=1,font=1.5, xpd=T)
    
    par(mfg = c(4, 2));
    plotGFR_dapa(results_GFR, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "H",cex=1,font=1.5, xpd=T)
    
    print(paste("Mean GFR with placebo and SGLT2 expression", parVal, ":", results_GFR_mean_placebo[450+xTime_offset_placebo, 2] * 1000));
    print(paste("Mean GFR with dapa and SGLT2 expression", parVal, ":", results_GFR_mean_dapa[450+xTime_offset_dapa, 2]*1000));
    print(paste("GFR difference without and with dapagliflozin [mL/min]:", (results_GFR_mean_placebo[450+xTime_offset_placebo, 2] - results_GFR_mean_dapa[450+xTime_offset_dapa, 2])*1000));
    
    idx = idx + 1;
  }
  par(mfg = c(3, 1));
  legend("topleft", c("SGLT2 100%", "SGLT2 80%", "SGLT2 120%"),
         lty=ltyArr,
         pch=c(NA,NA,NA),
         bty="n", y.intersp = 1.1)
}

#Variation of the proximal tubule (PT) volume
sim_V_PT = function(){
  #The parameter will be multiplied by these values.
  parFolds = c(1, 0.8, 1.2);
  
  initStruct <- list();
  #By initializing formula parameters, these parameters get "fixed", i.e., their values are constant and will not change.
  #e.g., blood flow rate will not be increased with the increased kidney volume, because the formula is replaced by the constant value.
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Volume", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Blood flow rate", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Fluid recirculation flow rate", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Cortex|Glomerulus|Volume", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Cortex|PCT|Volume", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Medulla|PST|Volume", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Medulla|HenleLoop_asc|Volume", initializeIfFormula = "always", 
                              initStruct = initStruct)
  initStruct <- initParameter(path_id = "*|Organism|Kidney|Medulla|HenleLoop_desc|Volume", initializeIfFormula = "always", 
                              initStruct = initStruct)
  currDCI = initSimulation(XML = paste0(simFolder, modelName, ".xml"), ParamList = initStruct);
  currDCI = setSimulationTime(timepoints = {0:(xTime_offset_dapa+450)}, currDCI);
  
  #Get the default values for PCT, PST, and whole kidney volumes.
  v_PCT_init = getParameter(path_id = "*|Organism|Kidney|Cortex|PCT|Volume", DCI_Info = currDCI)$Value;
  v_PST_init = getParameter(path_id = "*|Organism|Kidney|Medulla|PST|Volume", DCI_Info = currDCI)$Value;
  v_kid_init = getParameter(path_id = "*|Organism|Kidney|Volume", DCI_Info = currDCI)$Value;
  idx = 1;
  
  #Run simulations with variated parameter value.
  for (parFold in parFolds){
    #Change the volumes of PCT, PST, and the whole kidney.
    currDCI = setParameter(v_PCT_init * parFold, path_id = "*|Organism|Kidney|Cortex|PCT|Volume", DCI_Info = currDCI);
    currDCI = setParameter(v_PST_init * parFold, path_id = "*|Organism|Kidney|Medulla|PST|Volume", DCI_Info = currDCI);
    currDCI = setParameter(v_kid_init + (v_PCT_init + v_PST_init) * (parFold-1), path_id = "*|Organism|Kidney|Volume", DCI_Info = currDCI);

    #Run the simulation.
    currDCI = processSimulation(DCI_Info = currDCI);
    results_GFR = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
    results_GFR_mean_placebo = getSimulationResult(path_id = "*|Organism|Kidney|DEBUG_GFR_mean_placebo", DCI_Info = currDCI);
    results_GFR_mean_dapa = getSimulationResult(path_id = "*|Organism|Kidney|DEBUG_GFR_mean_dapa", DCI_Info = currDCI);
    results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
    
    par(mfg = c(1, 1));
    plotUGE_placebo(results_UGE, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "A",cex=1,font=1.5, xpd=T)
    
    par(mfg = c(1, 2));
    plotUGE_dapa(results_UGE, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "B",cex=1,font=1.5, xpd=T)
    
    par(mfg = c(2, 1));
    plotGFR_placebo(results_GFR, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "C",cex=1,font=1.5, xpd=T)
    
    par(mfg = c(2, 2));
    plotGFR_dapa(results_GFR, idx);
    par(new = TRUE)
    plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
    text(par("usr")[1] - par("usr")[1]*0.18, par("usr")[4], "D",cex=1,font=1.5, xpd=T)
    
    print(paste("Mean GFR with placebo and PT volume", parFold, "fold:", results_GFR_mean_placebo[450+xTime_offset_placebo, 2] * 1000));
    print(paste("Mean GFR with dapa and PT volume", parFold, "fold:", results_GFR_mean_dapa[450+xTime_offset_dapa, 2]*1000));
    print(paste("GFR difference without and with dapagliflozin [mL/min]:", (results_GFR_mean_placebo[450+xTime_offset_placebo, 2] - results_GFR_mean_dapa[450+xTime_offset_dapa, 2])*1000));
    
    idx = idx + 1;
  }
  par(mfg = c(1, 1));
  legend("topleft", c("V_PT 100%", "V_PT 80%", "V_PT 120%"),
         lty=ltyArr,
         pch=c(NA,NA,NA),
         bty="n", y.intersp = 1.1)
}

# png(file=paste0(figureFolder, "Figure_3.png"),
#     width=17.3,
#     height=14.38,
#     units = "cm",
#     res=300,
#     pointsize=8)
pdf(file=paste0(figureFolder, "Figure_3.pdf"),
    width=6.8,
    height=5.65,
    pointsize=8)

par(mfrow=c(4,2), cex=1,mar=c(2, 4, 0.5, 0.1), oma = c(1,0,1,0))

sim_V_PT();
sim_SGLT2();
# sim_Na();
mtext(c("Placebo","Dapagliflozin"), 3, 0, outer = TRUE, at = c(0.29, 0.79))
mtext(c("Plasma glucose [mmol/L]","Plasma glucose [mmol/L]"), 1, 0, outer = TRUE, at = c(0.29, 0.79))
dev.off();

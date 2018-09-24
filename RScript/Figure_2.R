#This script simulates the stepped hyperglycemic clamp (SHC) and generates Figure 2.

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
xVals_dapa = c(5.5, 9, 12.5, 14.6, 17.1, 19.6, 22.1, 24.1, 28.5, 30.4);
xVals_target = c(5.6, 8.3, 11.1, 13.9, 16.7, 19.4, 22.2, 25, 27.8, 30.5)
#These are the time points at which glomerular filtration rate (GFR) and urinary glucose excretion (UGE) are
#read out.
xTime = c(40, 80, 120, 160, 200, 240, 280, 320, 360, 400);
#Time at which SHC simulation without dapagliflozin administration starts.
xTime_offset_placebo = 240;
#Time at which SHC simulation with dapagliflozin administration starts.
xTime_offset_dapa = 10140;

#Name of the model xml-file.
modelName = "Dapagliflozin_2013_SHC_healthy";
#The list "initStruct" stores the paths of the parameters to be initialized (and changed later on).
initStruct <- list();
#Initialize the simulation with the parameters defined in "initStruct".
currDCI = initSimulation(XML = paste0(simFolder, modelName, ".xml"), ParamList = initStruct, whichInitParam = "allNonFormula");
#Set simulation time.
currDCI = setSimulationTime(timepoints = {0:(xTime_offset_dapa+450)}, DCI_Info = currDCI);
#Run the simulation.
currDCI = processSimulation(DCI_Info = currDCI);
#Get the GFR
results_GFR = getSimulationResult(path_id = "*|Organism|Kidney|GFR", DCI_Info = currDCI);
#Get the mean GFR during the time without dapagliflozin administration.
results_GFR_mean_placebo = getSimulationResult(path_id = "*|Organism|Kidney|DEBUG_GFR_mean_placebo", DCI_Info = currDCI);
#Get the mean GFR during the time with dapagliflozin administration.
results_GFR_mean_dapa = getSimulationResult(path_id = "*|Organism|Kidney|DEBUG_GFR_mean_dapa", DCI_Info = currDCI);
#Get the UGE
results_UGE = getSimulationResult(path_id = "*|Organism|Kidney|Urine|Glucose", DCI_Info = currDCI);
#Get water reabsortion rate in PCT.
results_rWat_PCT = getSimulationResult(path_id = "*|Organism|Kidney|Cortex|PCT|WaterReabsorption", DCI_Info = currDCI);
#Get water reabsortion rate in PST.
results_rWat_PST = getSimulationResult(path_id = "*|Organism|Kidney|Medulla|PST|WaterReabsorption", DCI_Info = currDCI);
#Get water reabsortion rate in the descending limb of Henle Loop.
results_rWat_HL = getSimulationResult(path_id = "*|Organism|Kidney|Medulla|HenleLoop_desc|WaterReabsorption", DCI_Info = currDCI);
#Get sodium concentration in the descending limb of Henle Loop.
results_Na_HL_desc = getSimulationResult(path_id = "*|Organism|Kidney|Medulla|HenleLoop_desc|C_Na", DCI_Info = currDCI);
#Get sodium concentration in the ascending limb of Henle Loop.
results_Na_HL_asc = getSimulationResult(path_id = "*|Organism|Kidney|Medulla|HenleLoop_asc|C_Na", DCI_Info = currDCI);
#Get filtrate flow rate in the ascending limb of Henle Loop.
results_Q_HL = getSimulationResult(path_id = "*|Organism|Kidney|Medulla|HenleLoop_asc|Filtrate flow rate", DCI_Info = currDCI);
#Concentration of sodium in arterial blood plasma
C_Na_AB = getParameter("*|Organism|ArterialBlood|Plasma|C_Na", DCI_Info = currDCI);


#Get sodium reabsortion rate in PCT.
results_NaReabs_PCT = getSimulationResult(path_id = "*|Organism|Kidney|Cortex|PCT|Na_Reabsorption", DCI_Info = currDCI);
#Get sodium reabsortion rate in PST.
results_NaReabs_PST = getSimulationResult(path_id = "*|Organism|Kidney|Medulla|PST|Na_Reabsorption", DCI_Info = currDCI);
#Get water reabsortion rate in the ascending limb of Henle Loop.
results_NaReabs_HL = getSimulationResult(path_id = "*|Organism|Kidney|Medulla|HenleLoop_asc|Na_Reabsorption", DCI_Info = currDCI);


print(paste("Mean GFR with placebo [mL/min]:", results_GFR_mean_placebo[450+xTime_offset_placebo, 2] * 1000));
print(paste("Mean GFR with dapagliflozin [mL/min]:", results_GFR_mean_dapa[450+xTime_offset_dapa, 2]*1000));
print(paste("GFR difference without and with dapagliflozin [mL/min]:", (results_GFR_mean_placebo[450+xTime_offset_placebo, 2] - results_GFR_mean_dapa[450+xTime_offset_dapa, 2])*1000));

pchArr = c(21, 4, NA, NA)
ltyArr = c(1, 2, 3, 4)

#Create a figure of GFR.
plotGFR = function(){
  idx = 1;
  #Plot GFR without dapagliflozin.
  plot(0:440, results_GFR[, 2][0:440+(xTime_offset_placebo-40)]*1e3,
       log="", type="l", xlab="Plasma glucose [mmol/L]", ylab="GFR [mL/min]",
       ylim=c(80, 160), xaxt="n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx]);
  
  #Plot GFR with dapagliflozin.
  idx = idx + 1;
  points(0:440, results_GFR[, 2][0:440+(xTime_offset_dapa-40)]*1e3,
         type="l",
         lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  
  axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
  axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
  
  legend("topleft", c("Placebo", "10 mg dapa"),
         lty=ltyArr,
         pch=c(NA, NA),
         bty="n", y.intersp = 1.1)
}

#Create a figure of UGE.
plotUGE = function(){
  #Molecular weight of glucose.
  mw = 180;
  #Conversion factor from reported mg/dL to simulated ?mol/L.
  convFac = 1e-6*mw;
  #yVals are the amounts of glucose at defined time points.
  yVals_placebo = c();
  yVals_dapa = c();
  for (time in xTime){
    yVals_placebo = c(yVals_placebo, results_UGE[time + xTime_offset_placebo, 2]);
    yVals_dapa = c(yVals_dapa, results_UGE[time + xTime_offset_dapa, 2]);
  }
  yAxtMarks = c(0, 5, 10, 15, 20, 25, 30);
  
  idx = 1;
  plot(xVals_placebo, yVals_placebo*convFac,
       log="", type="l", xlab="Plasma glucose [mmol/L]", ylab="UGE over 40 minutes [g]",
       ylim=c(0, 35), yaxt="n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  points(xVals_placebo, dataDeFronzo_healthy_placebo$Amount*convFac, pch = pchArr[idx], lwd=1, type="p", lty=ltyArr[idx])
  arrows(xVals_placebo, dataDeFronzo_healthy_placebo$Amount*convFac - dataDeFronzo_healthy_placebo$Error*convFac,
         xVals_placebo, dataDeFronzo_healthy_placebo$Amount*convFac + dataDeFronzo_healthy_placebo$Error*convFac,
         length=0.05, angle=90, code=3, lwd=1)
  axis(2, at = yAxtMarks, labels = yAxtMarks)
  
  idx = idx + 1;
  points(xVals_dapa, yVals_dapa*convFac,
         type="l",
         lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  points(xVals_dapa, dataDeFronzo_healthy_dapa$Amount*convFac, pch = pchArr[idx], lwd=1, type="p", lty=ltyArr[idx])
  arrows(xVals_dapa, dataDeFronzo_healthy_dapa$Amount*convFac - dataDeFronzo_healthy_dapa$Error*convFac,
         xVals_dapa, dataDeFronzo_healthy_dapa$Amount*convFac + dataDeFronzo_healthy_dapa$Error*convFac,
         length=0.05, angle=90, code=3, lwd=1)
  
  legend("topleft", c("Placebo", "10 mg dapa"),
         lty=ltyArr,
         pch=pchArr,
         bty="n", y.intersp = 1.1)
}

#Create a plot of total fractional water reabsorption
plotWaterReabsorption = function(){
  idx = 1;
  waterReabsorption = results_rWat_PCT;
  waterReabsorption[, 1] = results_rWat_PCT[, 1];
  #Sum of water reabsorption rates in all nephorn areas in % of total GFR.
  waterReabsorption[, 2] = (results_rWat_PCT[, 2] + results_rWat_PST[, 2] + results_rWat_HL[, 2]) / results_GFR[, 2] * 100;
  
  plot(0:440, waterReabsorption[, 2][0:440+(xTime_offset_placebo-40)],
       log="", type="l", xlab="Plasma glucose [mmol/L]", ylab="Water reabsorption [%]",
       ylim=c(80, 100), xaxt="n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(0:440, waterReabsorption[, 2][0:440+(xTime_offset_dapa-40)],
         type="l",
         lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  
  axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
  axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
  
  legend("topleft", c("Placebo", "10 mg dapa"),
         lty=ltyArr,
         pch=c(NA, NA),
         bty="n", y.intersp = 1.1)
}

#Plot filtrate flow rate in the ascending limb of Henle Loop.
plotQ = function(){
  idx = 1;
  plot(0:440, results_Q_HL[, 2][0:440+(xTime_offset_placebo-40)]*1e3,
       log="", type="l", xlab="Plasma glucose [mmol/L]", ylab="Q_HL_asc [mL/min]",
       ylim=c(0, 20), xaxt="n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(0:440, results_Q_HL[, 2][0:440+(xTime_offset_dapa-40)]*1e3,
         type="l",
         lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  
  axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
  axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
  
  legend("topleft", c("Placebo", "10 mg dapa"),
         lty=ltyArr,
         pch=c(NA, NA),
         bty="n", y.intersp = 1.1)
}

#Plot sodium flow rate out of the ascending limb of Henle Loop.
plotQ_Na_HL_asc = function(){
  idx = 1;
  plot(0:440, (results_Q_HL[, 2][0:440+(xTime_offset_placebo-40)]) * (results_Na_HL_asc[, 2][0:440+(xTime_offset_placebo-40)]),
       log="", type="l", xlab="Plasma glucose [mmol/L]", ylab="Q_Na_HL_asc [µmol/min]",
       ylim=c(0, 1000),
       xaxt="n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(0:440, (results_Q_HL[, 2][0:440+(xTime_offset_dapa-40)]) * (results_Na_HL_asc[, 2][0:440+(xTime_offset_dapa-40)]),
         type="l",
         lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  
  axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
  axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
  
  legend("topleft", c("Placebo", "10 mg dapa"),
         lty=ltyArr,
         pch=c(NA, NA),
         bty="n", y.intersp = 1.1)
}

#Create a plot of total fractional sodium reabsorption
plotNaReabsorption = function(){
  idx = 1;
  NaReabsorption = results_NaReabs_PCT;
  NaReabsorption[, 1] = results_NaReabs_PCT[, 1];
  #Sodium reabsorption rate in ascending Henle Loop in % of total GFR times AB sodium concentration.
  NaReabsorption[, 2] = (results_NaReabs_HL[, 2]) / (results_GFR[, 2] * C_Na_AB$Value) * 100;
  
  plot(0:440, NaReabsorption[, 2][0:440+(xTime_offset_placebo-40)],
       log="", type="l", xlab="Plasma glucose [mmol/L]", ylab="Na reabsorption [%]",
       #ylim=c(10000, 20000),
       #ylim=c(90, 100),
       xaxt="n", lwd=1, pch = pchArr[idx], lty=ltyArr[idx]);
  
  idx = idx + 1;
  points(0:440, NaReabsorption[, 2][0:440+(xTime_offset_dapa-40)],
         type="l",
         lwd=1, pch = pchArr[idx], lty=ltyArr[idx])
  
  axis(1, at = c(0, xTime)+20, labels = c(4, xVals_target), tick = FALSE);
  axis(1, at = c(0, xTime), labels = FALSE, tick = TRUE);
  
  legend("topleft", c("Placebo", "10 mg dapa"),
         lty=ltyArr,
         pch=c(NA, NA),
         bty="n", y.intersp = 1.1)
}

# Plot into a png file.
# png(file=paste0(figureFolder, "Figure_2.png"),
#     width=17.3,
#     height=11.46,
#     units = "cm",
#     res=400,
#     pointsize=8)

pdf(file=paste0(figureFolder, "Figure_2.pdf"),
    width=6.8,
    height=4.5,
    pointsize=8)

par(mfrow = c(2, 2), cex=1, oma=c(0,0,0,0))
par(mar=c(4,4,0.4,0.1))

plotUGE();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "A",cex=1,font=1.5, xpd=T)
plotGFR();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "B",cex=1,font=1.5, xpd=T)

plotQ_Na_HL_asc()
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "C",cex=1,font=1.5, xpd=T)
plotWaterReabsorption()
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "D",cex=1,font=1.5, xpd=T)

dev.off();

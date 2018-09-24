#This script simulates the PBPK model of dapagliflozin and generates Supplemental Figure 1. The *.xml-model files are generated
#from the respective simulations in the MoBi-project "Renal_Dapagliflozin_PK_Balazki".

#The below loaded libraries must be installed before running this script.
library(MoBiToolboxForR)
library(parallel)
library(RColorBrewer)

#Define here the path to the folder with the experimental data extracted from literature.
#For the description of the datasets, refer to the main text.
expDataFolder = "..\\ExpData\\";
#Define here the path to the folder with the simulation files.
simFolder = "..\\xml\\";
#Define here the path to the folder where the output figure will be stored.
figureFolder = "..\\Figures\\";
decSym = '.';

#Define the list of the colors to use.
cols = c("black", brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"));

#Read experimental data
dataBoulton_iv = read.table(paste0(expDataFolder, "Boulton_2013_iv.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataObermeier_Beagle_iv = read.table(paste0(expDataFolder, "Obermeier_2010_Beagle.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataObermeier_Rat_iv = read.table(paste0(expDataFolder, "Obermeier_2010_Rat.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataObermeier_Monkey_iv = read.table(paste0(expDataFolder, "Obermeier_2010_Monkey.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataObermeier_Human_po = read.table(paste0(expDataFolder, "Obermeier_2010_Human.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_2.5 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_2.5mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_5 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_5mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_10 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_10mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_20 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_20mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_50 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_50mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_100 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_100mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_250 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_250mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiSAD_500 = read.table(paste0(expDataFolder, "Komoroski_2009_SAD_500mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiMAD_2.5 = read.table(paste0(expDataFolder, "Komoroski_2009_MAD_2.5mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiMAD_10 = read.table(paste0(expDataFolder, "Komoroski_2009_MAD_10mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiMAD_20 = read.table(paste0(expDataFolder, "Komoroski_2009_MAD_20mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiMAD_50 = read.table(paste0(expDataFolder, "Komoroski_2009_MAD_50mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKomoroskiMAD_100 = read.table(paste0(expDataFolder, "Komoroski_2009_MAD_100mg.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)
dataKasichayanula_2011 = read.table(paste0(expDataFolder, "Kasichayanula_2011.csv"), header = TRUE, sep=";", as.is = TRUE, stringsAsFactors = FALSE)

#Runs the simulation with the file name "simName", and returns the result of the simulation with the path "resultsPath"
executeSim = function(simName, resultsPath){
  #The list "initStruct" stores the paths of the parameters to be initialized (and changed later on).
  initStruct <- list();
  #Initialize the simulation with the parameters defined in "initStruct".
  currDCI = initSimulation(XML = paste0(simFolder, simName, ".xml"), ParamList = initStruct, whichInitParam = "none");
  #Run the simulation.
  currDCI = processSimulation(DCI_Info = currDCI);
  #Get the specified results and return.
  results = getSimulationResult(resultsPath, DCI_Info = currDCI);
  return(results);
}

#Run all the simulations specified in the list "modelNames" in parallel, and return the results.
sim_Models = function(modelNames, resultsPath){
  # Calculate the number of cores
  no_cores = detectCores() - 1;
  
  # Initiate cluster
  cl = makeCluster(no_cores);
  #Make the variable "simFolder" visible to the cluster.
  clusterExport(cl, "simFolder");
  #Make the MoBi-Toolbox visible to the cluster.
  clusterEvalQ(cl, library(MoBiToolboxForR));
  #Run the method "executeSim" in parallel.
  results = parSapply(cl, modelNames, executeSim, resultsPath);
  stopCluster(cl)
  return(results);
}

#Define the path of the result - concentration of dapagliflozin in peripheral venous blood plasma.
resultsPath = "*|Organism|PeripheralVenousBlood|Dapagliflozin|Plasma \\(Peripheral Venous Blood\\)";

#Names of the xml-model files to be simulated.
modelNames = c("Komoroski_2009_SAD_100mg",
               "Komoroski_2009_SAD_10mg",
               "Komoroski_2009_SAD_2.5mg",
               "Komoroski_2009_SAD_20mg",
               "Komoroski_2009_SAD_250mg",
               "Komoroski_2009_SAD_500mg",
               "Komoroski_2009_SAD_50mg",
               "Komoroski_2009_SAD_5mg",
               "Komoroski_2009_MAD_100mg",
               "Komoroski_2009_MAD_10mg",
               "Komoroski_2009_MAD_2.5mg",
               "Komoroski_2009_MAD_20mg",
               "Komoroski_2009_MAD_50mg",
               "Kasichayanula_2011_healthy",
               "Obermeier_2010_Human_po",
               "Boulton_iv",
               "Obermeier_2010_Beagle_iv",
               "Obermeier_2010_Monkey_iv",
               "Obermeier_2010_rat_iv"
               );

#Simulate the models.
results = sim_Models(modelNames, resultsPath);

#Define converstion factor from simulated ?mol/l to observed ng/ml. mw: molecular weight of dapagliflozin.
mw = 408.87;
convFac = 1e-6*mw*1e9*1e-3;

pchArr = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
ltyArr = c(1, 2, 4, 5, 6, 7, 8, 10)

#Create a plot of the single ascending dose (SAD) simulations. Plot simulted and observed concentrat-time profiles.
plotSAD = function(){
  modelNames = c("Komoroski_2009_SAD_2.5mg",
                 "Komoroski_2009_SAD_5mg",
                 "Komoroski_2009_SAD_10mg",
                 "Komoroski_2009_SAD_20mg",
                 "Komoroski_2009_SAD_50mg",
                 "Komoroski_2009_SAD_100mg",
                 "Komoroski_2009_SAD_250mg",
                 "Komoroski_2009_SAD_500mg"
                 );
  
  #define y-axis marks
  yAxtMarks = c(1, 10, 100, 1000, 10000)
  #turn warnings off
  oldw <- getOption("warn")
  options(warn = -1)
  
  idx = 1;
  #Plot simulated results first. The time values are divided by 60 to convert from minutes to hours.
  #The concentration values are multiplied by the conversion factor to convert to ng/ml.
  plot(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
       log="y", type="l", xlab="Time [h]", ylab="Plasma dapagliflozin [ng/ml]",
       xlim=c(0,24),
       ylim=c(min(yAxtMarks), max(yAxtMarks)), yaxt="n", lwd=1,
       col = cols[idx]);
  #Plot observed data points incl. SD
  points(dataKomoroskiSAD_2.5$Time/60, dataKomoroskiSAD_2.5$Concentration,
         pch = pchArr[idx], lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_2.5$Time/60, dataKomoroskiSAD_2.5$Concentration - dataKomoroskiSAD_2.5$Error,
         dataKomoroskiSAD_2.5$Time/60, dataKomoroskiSAD_2.5$Concentration + dataKomoroskiSAD_2.5$Error,
         length=0.05, angle=90, code=3, lwd=1,
         col = cols[idx])
  axis(2, at = yAxtMarks, labels = yAxtMarks)
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiSAD_5$Time/60, dataKomoroskiSAD_5$Concentration,
         pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_5$Time/60, dataKomoroskiSAD_5$Concentration - dataKomoroskiSAD_5$Error,
         dataKomoroskiSAD_5$Time/60, dataKomoroskiSAD_5$Concentration + dataKomoroskiSAD_5$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiSAD_10$Time/60, dataKomoroskiSAD_10$Concentration,
         pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_10$Time/60, dataKomoroskiSAD_10$Concentration - dataKomoroskiSAD_10$Error,
         dataKomoroskiSAD_10$Time/60, dataKomoroskiSAD_10$Concentration + dataKomoroskiSAD_10$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiSAD_20$Time/60, dataKomoroskiSAD_20$Concentration, pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_20$Time/60, dataKomoroskiSAD_20$Concentration - dataKomoroskiSAD_20$Error,
         dataKomoroskiSAD_20$Time/60, dataKomoroskiSAD_20$Concentration + dataKomoroskiSAD_20$Error,
         length=0.05, angle=90, code=3, lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiSAD_50$Time/60, dataKomoroskiSAD_50$Concentration, pch=pchArr[idx], lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_50$Time/60, dataKomoroskiSAD_50$Concentration - dataKomoroskiSAD_50$Error,
         dataKomoroskiSAD_50$Time/60, dataKomoroskiSAD_50$Concentration + dataKomoroskiSAD_50$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiSAD_100$Time/60, dataKomoroskiSAD_100$Concentration,
         pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_100$Time/60, dataKomoroskiSAD_100$Concentration - dataKomoroskiSAD_100$Error,
         dataKomoroskiSAD_100$Time/60, dataKomoroskiSAD_100$Concentration + dataKomoroskiSAD_100$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiSAD_250$Time/60, dataKomoroskiSAD_250$Concentration,
         pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_250$Time/60, dataKomoroskiSAD_250$Concentration - dataKomoroskiSAD_250$Error,
         dataKomoroskiSAD_250$Time/60, dataKomoroskiSAD_250$Concentration + dataKomoroskiSAD_250$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiSAD_500$Time/60, dataKomoroskiSAD_500$Concentration,
         pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiSAD_500$Time/60, dataKomoroskiSAD_500$Concentration - dataKomoroskiSAD_500$Error,
         dataKomoroskiSAD_500$Time/60, dataKomoroskiSAD_500$Concentration + dataKomoroskiSAD_500$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  #Paint a grey rectangular to highlight concentrations outside the precise assay quantification range.
  rect(0,0.1,25,10,col = rgb(0.5,0.5,0.5,1/4), border = NA)
  rect(0,1000,25,20000,col = rgb(0.5,0.5,0.5,1/4), border = NA)
  
  legend("topright", c("2.5 mg SAD", "5   mg SAD", "10 mg SAD", "20 mg SAD", "50   mg SAD", "100 mg SAD", "250 mg SAD", "500 mg SAD"),
         lty=ltyArr,
         pch=pchArr,
         col = cols,
         ncol = 2,
         bty="n", y.intersp = 1.1)
  
  #turn warnings back on
  options(warn = oldw)
}

#Create a plot of the 14th day of the multiple ascending dose (MAD) simulations. Plot simulted and observed concentrat-time profiles.
plotMAD_day14 = function(){
  modelNames = c("Komoroski_2009_MAD_2.5mg",
                 "dummy",
                 "Komoroski_2009_MAD_10mg",
                 "Komoroski_2009_MAD_20mg",
                 "Komoroski_2009_MAD_50mg",
                 "Komoroski_2009_MAD_100mg");
  
  #define y-axis marks
  yAxtMarks = c(1, 10, 100, 1000, 10000)
  oldw <- getOption("warn")
  #turn warnings off
  options(warn = -1)
  idx = 1;
  plot(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
       log="y", type="l", xlab="Time [h]", ylab="Plasma dapagliflozin [ng/ml]",
       xlim=c(24*13, 24*14),
       ylim=c(min(yAxtMarks), max(yAxtMarks)), yaxt="n", lwd=1,
       col = cols[idx])
  points(dataKomoroskiMAD_2.5$Time, dataKomoroskiMAD_2.5$Concentration,
         pch = pchArr[idx], lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiMAD_2.5$Time, dataKomoroskiMAD_2.5$Concentration - dataKomoroskiMAD_2.5$Error,
         dataKomoroskiMAD_2.5$Time, dataKomoroskiMAD_2.5$Concentration + dataKomoroskiMAD_2.5$Error,
         length=0.05, angle=90, code=3, lwd=1,
         col = cols[idx])
  axis(2, at = yAxtMarks, labels = yAxtMarks)
  
  idx = idx + 1;
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l", 
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiMAD_10$Time, dataKomoroskiMAD_10$Concentration,
         pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiMAD_10$Time, dataKomoroskiMAD_10$Concentration - dataKomoroskiMAD_10$Error,
         dataKomoroskiMAD_10$Time, dataKomoroskiMAD_10$Concentration + dataKomoroskiMAD_10$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiMAD_20$Time, dataKomoroskiMAD_20$Concentration, pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiMAD_20$Time, dataKomoroskiMAD_20$Concentration - dataKomoroskiMAD_20$Error,
         dataKomoroskiMAD_20$Time, dataKomoroskiMAD_20$Concentration + dataKomoroskiMAD_20$Error,
         length=0.05, angle=90, code=3, lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l", lty=ltyArr[idx], lwd=1,
         col = cols[idx])
  points(dataKomoroskiMAD_50$Time, dataKomoroskiMAD_50$Concentration, pch=pchArr[idx], lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiMAD_50$Time, dataKomoroskiMAD_50$Concentration - dataKomoroskiMAD_50$Error,
         dataKomoroskiMAD_50$Time, dataKomoroskiMAD_50$Concentration + dataKomoroskiMAD_50$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l",
         lty=ltyArr[idx],
         lwd=1,
         col = cols[idx])
  points(dataKomoroskiMAD_100$Time, dataKomoroskiMAD_100$Concentration,
         pch=pchArr[idx],
         lwd=1,
         col = cols[idx])
  arrows(dataKomoroskiMAD_100$Time, dataKomoroskiMAD_100$Concentration - dataKomoroskiMAD_100$Error,
         dataKomoroskiMAD_100$Time, dataKomoroskiMAD_100$Concentration + dataKomoroskiMAD_100$Error,
         length=0.05, angle=90, code=3,
         lwd=1,
         col = cols[idx])
  
  rect(300,0.1,340,1,col = rgb(0.5,0.5,0.5,1/4), border = NA)
  rect(300,1000,340,20000,col = rgb(0.5,0.5,0.5,1/4), border = NA)
  
  
  legend("topright", c("2.5 mg MAD", "10 mg MAD", "20 mg MAD", "50   mg MAD", "100 mg MAD"),
         ncol = 2,
         lty = ltyArr[c(1, 3:length(ltyArr))],
         pch = pchArr[c(1, 3:length(pchArr))],
         col = cols[c(1, 3:length(cols))],
         bty="n", y.intersp = 1.1)
  
  #turn warnings back on
  options(warn = oldw)
}

#Create a plot of the intravenous (iv) simulations. Plot simulted and observed concentrat-time profiles.
plot_IV = function(){
  modelNames = c("Boulton_iv",
                 "Obermeier_2010_Beagle_iv",
                 "Obermeier_2010_Monkey_iv",
                 "Obermeier_2010_rat_iv");
  
  #define y-axis marks
  yAxtMarks = c(0.001, 0.1, 10, 1000, 100000)
  oldw <- getOption("warn")
  #turn warnings off
  options(warn = -1)
  idx = 1;
  plot(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
       log="y", type="l", xlab="Time [h]", ylab="Plasma dapagliflozin [ng/ml]",
       xlim = c(0, 50),
       ylim=c(0.001, 100000),
       yaxt="n", lwd=1,
       col = cols[idx])
  points(dataBoulton_iv$Time, dataBoulton_iv$Concentration, pch = pchArr[idx], lwd=1,
         col = cols[idx])
  arrows(dataBoulton_iv$Time, dataBoulton_iv$Concentration - dataBoulton_iv$Error,
         dataBoulton_iv$Time, dataBoulton_iv$Concentration + dataBoulton_iv$Error,
         length=0.05, angle=90, code=3, lwd=1,
         col = cols[idx])
  axis(2, at = yAxtMarks, labels = c("0.001","0.1","10","1000","100000"))
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l", lty=ltyArr[idx], lwd=1,
         col = cols[idx])
  points(dataObermeier_Beagle_iv$Time, dataObermeier_Beagle_iv$Concentration, pch=pchArr[idx], lwd=1,
         col = cols[idx])
  arrows(dataObermeier_Beagle_iv$Time, dataObermeier_Beagle_iv$Concentration - dataObermeier_Beagle_iv$Error,
         dataObermeier_Beagle_iv$Time, dataObermeier_Beagle_iv$Concentration + dataObermeier_Beagle_iv$Error,
         length=0.05, angle=90, code=3, col = cols[idx], lwd=1)
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l", lty=ltyArr[idx], lwd=1,
         col = cols[idx])
  points(dataObermeier_Monkey_iv$Time, dataObermeier_Monkey_iv$Concentration, pch=pchArr[idx], lwd=1,
         col = cols[idx])
  arrows(dataObermeier_Monkey_iv$Time, dataObermeier_Monkey_iv$Concentration - dataObermeier_Monkey_iv$Error, 
         dataObermeier_Monkey_iv$Time, dataObermeier_Monkey_iv$Concentration + dataObermeier_Monkey_iv$Error, 
         length=0.05, angle=90, code=3, lwd=1,
         col = cols[idx])
  
  idx = idx + 1;
  points(results[[modelNames[idx]]][,1]/60, results[[modelNames[idx]]][,2]*convFac,
         type="l", lty=ltyArr[idx], lwd=1,
         col = cols[idx])
  points(dataObermeier_Rat_iv$Time, dataObermeier_Rat_iv$Concentration, pch=pchArr[idx], lwd=1,
         col = cols[idx])
  #arrows(dataObermeier_Rat_iv$Time, dataObermeier_Rat_iv$Concentration - dataObermeier_Rat_iv$Error, dataObermeier_Rat_iv$Time, dataObermeier_Rat_iv$Concentration + dataObermeier_Rat_iv$Error, length=0.05, angle=90, code=3, col="green", lwd=1)
  
  legend("topright", c("Dog 6.6 mg/kg", "Monkey 6 mg/kg", "Rat 1 mg/kg", "Human 80 µg"),
         lty=c(2, 4, 5, 1),
         pch=c(1, 2, 3, 0),
         col = c(cols[2:4], "black"),
         bty="n", y.intersp = 1.1)
  
  #turn warnings back on
  options(warn = oldw)
}

#Create a simulated-vs-observed plot
simulatedVsObserved = function(){
  #Determine the maximal simulated or reported concentration value to set the limits of the axes.
  maxVal = max(c(
    dataObermeier_Beagle_iv$Concentration..ng.ml.,
    dataObermeier_Monkey_iv$Concentration..ng.ml.,
    dataObermeier_Rat_iv$Concentration..ng.ml.,
    dataObermeier_Human_po$Concentration..ng.ml., dataKomoroskiSAD_2.5$Concentration..ng.ml.,
    dataKomoroskiSAD_5$Concentration..ng.ml., dataKomoroskiSAD_10$Concentration..ng.ml.,
    dataKomoroskiSAD_20$Concentration..ng.ml., dataKomoroskiSAD_50$Concentration..ng.ml.,
    dataKomoroskiSAD_100$Concentration..ng.ml., dataKomoroskiSAD_250$Concentration..ng.ml.,
    dataKomoroskiSAD_500$Concentration..ng.ml., dataKomoroskiMAD_2.5$Concentration..ng.ml.,
    dataKomoroskiMAD_10$Concentration..ng.ml., dataKomoroskiMAD_20$Concentration..ng.ml.,
    dataKomoroskiMAD_50$Concentration..ng.ml., dataKomoroskiMAD_100$Concentration..ng.ml.,
    dataKasichayanula_2011$Concentration..ng.ml.))

    for (result in results){
      maxVal = max(c(maxVal, max(result[,2])*convFac))
    }
  #This list stores all the deviations.
  devs = c();
  
  identity = function(x){
    x;
  }
  dev = 1
  plusDev = function(x){
    x * (1/0.5);
  }
  minusDev = function(x){
    x * 0.5;
  }
  
  #Plot identity line, representing perferct accordance.
  curve(identity(x), from = 0.01, to = maxVal, log="xy", xlab="Observed concentration [ng/ml]", ylab="Simulated concentration [ng/ml]")
  #Plot the 0.5 and 2 fold deviation lines.
  curve(plusDev(x), from = 0.001, to = maxVal, add = TRUE);
  curve(minusDev(x), from = 0.01, to = maxVal*2, add = TRUE);
  
  #The list devCurr stores the deviations of simulated values from the reported ones for the current data set.
  devCurr = c();
    idxPlot = 1;
  #Boulton_iv
  currData = dataBoulton_iv;
  currSim = results[["Boulton_iv"]];
  #For each dataset, iterate through the reported data points.
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    #Find the corresponding data point in the simulation results.
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac;
    #Draw a point with the coordinates equals to y = simulated value, x = observed value.
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    #Calculate the fold deviation of simulated from observed value.
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  #Print the number of points within a certain deviation range.
  print(paste('Bouton_iv: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('Bouton_iv: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #Obermeier_Beagle
  idxPlot = idxPlot+1;
  currData = dataObermeier_Beagle_iv;
  currSim = results[["Obermeier_2010_Beagle_iv"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('Obermeier_beagle: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('Obermeier_beagle: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #Obermeier_Monkey
  idxPlot = idxPlot+1;
  currData = dataObermeier_Monkey_iv;
  currSim = results[["Obermeier_2010_Monkey_iv"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('Obermeier_monkey: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('Obermeier_monkey: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #Obermeier_Rat
  idxPlot = idxPlot+1;
  currData = dataObermeier_Rat_iv;
  currSim = results[["Obermeier_2010_rat_iv"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('Obermeier_rat: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('Obermeier_rat: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #Obermeier_po
  idxPlot = idxPlot+1;
  currData = dataObermeier_Human_po;
  currSim = results[["Obermeier_2010_Human_po"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('Obermeier_hum: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('Obermeier_hum: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD2.5
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_2.5;
  currSim = results[["Komoroski_2009_SAD_2.5mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_2.5: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_2.5: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD5
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_5;
  currSim = results[["Komoroski_2009_SAD_5mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_5: nr of points:', length(currData$Time..min.)));
  print(paste('SAD_5: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_5: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD10
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_10;
  currSim = results[["Komoroski_2009_SAD_10mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_10: nr of points:', length(currData$Time..min.)));
  print(paste('SAD_10: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_10: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD20
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_20;
  currSim = results[["Komoroski_2009_SAD_20mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_20: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_20: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD50
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_50;
  currSim = results[["Komoroski_2009_SAD_50mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_50: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_50: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD100
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_100;
  currSim = results[["Komoroski_2009_SAD_100mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_100: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_100: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD250
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_250;
  currSim = results[["Komoroski_2009_SAD_250mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_250: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_250: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #SAD500
  idxPlot = idxPlot+1;
  currData = dataKomoroskiSAD_500;
  currSim = results[["Komoroski_2009_SAD_500mg"]];
  for (i in 1:length(currData$Time..min.)){
    dataTime = currData$Time..min.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime)==min(abs(currSim[,1] - dataTime)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('SAD_500: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('SAD_500: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #MAD2.5
  idxPlot = idxPlot+1;
  currData = dataKomoroskiMAD_2.5;
  currSim = results[["Komoroski_2009_MAD_2.5mg"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('MAD_2.5: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('MAD_2.5: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #MAD10
  idxPlot = idxPlot+1;
  currData = dataKomoroskiMAD_10;
  currSim = results[["Komoroski_2009_MAD_10mg"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('MAD_10: nr of points:', length(currData$Time..h.)));
  print(paste('MAD_10: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('MAD_10: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #MAD20
  idxPlot = idxPlot+1;
  currData = dataKomoroskiMAD_20;
  currSim = results[["Komoroski_2009_MAD_20mg"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('MAD_20: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('MAD_20: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #MAD50
  idxPlot = idxPlot+1;
  currData = dataKomoroskiMAD_50;
  currSim = results[["Komoroski_2009_MAD_50mg"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('MAD_50: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('MAD_50: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #MAD100
  idxPlot = idxPlot+1;
  currData = dataKomoroskiMAD_100;
  currSim = results[["Komoroski_2009_MAD_100mg"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('MAD_100: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('MAD_100: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  devCurr = c();
  
  #Kasichayanula_2011
  idxPlot = idxPlot+1;
  currData = dataKasichayanula_2011;
  currSim = results[["Kasichayanula_2011_healthy"]];
  for (i in 1:length(currData$Time..h.)){
    dataTime = currData$Time..h.[i];
    dataVal = currData$Concentration..ng.ml.[i];
    idx = which(abs(currSim[,1] - dataTime*60)==min(abs(currSim[,1] - dataTime*60)));
    simVal = currSim[idx,2]*convFac
    points(dataVal, simVal, pch = pchArr[idxPlot],
           col = cols[idxPlot])
    devCurr = c(devCurr, simVal / dataVal);
    devs = c(devs, simVal / dataVal);
  }
  print(paste('Kasichayanula: nr of points:', length(currData$Time..h.)));
  print(paste('Kasi_a: Deviation higher then 2x resp. 0.5x:', (sum(devCurr > (1/0.5)) + sum(devCurr < 0.5))));
  print(paste('Kasi_a: Deviation higher then 1/0.7x resp. 0.7x:', (sum(devCurr > (1/0.7)) + sum(devCurr < 0.7))));
  
  return(devs);
}

pdf(file=paste0(figureFolder, "Supp_Figure_1.pdf"),
           width=6.3,
           height=6.3,
           pointsize=8)

par(mfrow = c(2, 2), cex=1, oma=c(0,0,0,0))
par(mar=c(4,4,0.4,0.1))

plotSAD();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "A",cex=1,font=1.5, xpd=T)
plotMAD_day14();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "B",cex=1,font=1.5, xpd=T)
plot_IV();
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "C",cex=1,font=1.5, xpd=T)
devs = simulatedVsObserved();
print(paste('Total number of points: ', length(devs)))
par(new = TRUE)
plot(10,0, axes = FALSE, pch=NA, xlab = "", ylab = "")
text(par("usr")[1] - par("usr")[1]*0.25, par("usr")[4], "D",cex=1,font=1.5, xpd=T)
dev.off();

print(paste('Deviation higher then 2x resp. 0.5x:', (sum(devs > (1/0.5)) + sum(devs < 0.5))));
print(paste('Deviation higher then 1/0.7x resp. 0.7x:', (sum(devs > (1/0.7)) + sum(devs < 0.7))));
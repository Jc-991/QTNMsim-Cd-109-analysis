# QTNMsim-Cd-109-analysis
The codes within this resository are used to analyse the data generated from the QTNMsim code in this link https://github.com/QTNM/Electron-Tracking/tree/main/examples/Cd109source
The analyseraw.py is used to classify the particles from the raw data and have a general check at the distribution of selected value. eval and gval are values for the expected values of electrons and gammma rays. Here to make it more general, I have made them as inputs of function.
The analysegaussian.py can be used to find the Gaussian distribution for the raw and data analysis for the targeted peaks hence analyse the precision of the calibration. The detailed instructions for each function are in the script.

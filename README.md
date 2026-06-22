Team and Contributors

D-FENSE
Americo Cunha Jr (team leader)
Paulo Esquef (LNCC)
Emanuelle Paixão (LNCC)

2. Repository Structure

Main
|_>  Aggredated Data
     (Data used to train the model)
|_>  LNCC_ARp26_dengue  (model source and results)
     |_> validation1    (tranining for validation 1) 
         |_> matlab    (matlab scritps to generate the results)
         |_> spreadsheets  (CSV files with the results; one for each state)
         |_> plots   (PDF files with plots of the results; one for each state)
     |_> validation2  (see validation 1)
         |_> matlab
         |_> spreadsheets
         |_> plots
     |_> validation3 (see validation 1)
         |_> matlab
         |_> spreadsheets
         |_> plots
     |_> validation4  (see validation 1)
         |_> matlab
         |_> spreadsheets
         |_> plots
         
3. Libraries and Dependencies

Proprietary functions of Matlab toolboxes:
armcov.m   % matlab function for AR(p) model estimation
filter.m  % matlab function for linear convolution
filter2.m % matlab function for centered linear convolution
buffer.m  % matlab function for data organization into columns (block segmentation)
randn.m  % matlab function for random number generation from a Normal Distribution
prctile.m % matlab function for calculating the percentiles of a set of numbers
table.m  % matlab function for managing table creation
writetable.m % matlab function for creating a CSV file from a table

Extra function:
ssa_modPE.m  % Singular Spectrum Analysis method used for lowpass filtering the data (included in folder matlab)

4. Data and Variables
For each state, the input time series aggregates the sum of all dengue cases over a epiweek, in a given state. 
Data is organized in 27 CSV files, one per state, stored in folder 'Aggregated Data'.
Only the time series of dengue cases has been used. No climate data has been used.
  
5. Model Training

Model training only invoved regular estimation of an AR(p) parameters, based on the zero-mean log2 transformation of the original time series of dengue cases.  
For validation 1: all available data has been used, from EW 1 of 2010 to EW 25 of 2022 
For validation 2: all available data has been used, from EW 1 of 2010 to EW 25 of 2023 
For validation 3: all available data has been used, from EW 1 of 2010 to EW 25 of 2024 
For validation 4: all available data has been used, from EW 1 of 2010 to EW 25 of 2025 

CODE (example for validation 1) : 
cc=kcc(1:ind_v1);  % observed cases from EW 1 of 2010 to EW 25 of 2022
cclog=log2(cc); % apply a log2 function to the observed raw data 
mcc=mean(cclog); % calculate and save the mean of cclog
sig=cclog-mcc;  % zero-mean log2 observed data
p=92; % AR model order 
a=armcov(sig,p);  % AR model estimation via the modified 


5. Data Usage Restriction

After estimating the model parameters, we run the model, i.e., excite the model with appropriate randonly generated signals with length which is long enough to provide forecasts
from EW 26 of year 202X up to EW 52 of year 202(X+1). We then lowpass filter the attained reaults and, finally, crop the forecasts to comply with the specified range from EW 41 
of year 202X up to EW 40 of year 202(X+1). 


6. Predictive Uncertainty How are your prediction intervals computed?
We have carriedout a Monte Carlo simulation with 10000 runs, for each state. Out of the 10000 forecast realizations, for a given time instant, we calculate the specified percentiles (via prctile.m),
 produce the specified prediction intervals. 


8. References
  None. 

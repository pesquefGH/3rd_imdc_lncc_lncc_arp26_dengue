% Validation 1: LNCC ARP26 forecast model (high-order AR(p) model)
% Tasks of this script (for a given BR state):

% 1) Read the time series of dengue cases for a given BR state and estimate an
% high-order AR(p) model to the log2 mapping of the time-series.

% 2) Analyze the modeling error, obtained via inverse filtering,
% considering a seasonality of 52 epidemic weeks (EW). Error analysis will
% inform the tuning of a statistical distribution for excitation generation, 
% for time-series forecast purposes.

% 3) Obtain a set of initial conditions (IC) for the AR(p) model.

% 4) Run a Monte Carlo (MC) simulation for time-series forecast using the
% set of IC and randomly generated model excitation realizations.

% 5) Plot the obtained time-series forecast (mean, upper- and lower-bounds
% of the prediction intervals)

SUNDAYS = readtable('sunday_dates.csv');  % reads CSV with Sunday dates
S_dates=SUNDAYS{:,1};   % Sunday dates format YYYY-MM-DD

M = readtable(['IMDC2026_AggregatedData_',UF,'.csv']); 
% reads data related to the selected state

kcc=M{:,3};  % known dengue cases from EW 1 of 2010 to EW 52 of 2025


ind_v1=649;  % index of the end of EW 25 of 2022 (upper bound 
% observed data for validation 1)

cc=kcc(1:ind_v1);  % observed cases from EW 1 of 2010 to EW 25 of 2022

cc(cc==0)=0.1;  % replaces null value in cc with a small value. 
% This is because we will apply 
% a log2 function to the data values. 

lwea=2*52-25;   % number of EW to be forecast, i.e., from EW 26 of 2022 
% up to EW 52 2023

cclog=log2(cc); % apply a log2 function to the observed raw data 
mcc=mean(cclog); % calculate and save the mean of cclog
sig=cclog-mcc;  % zero-mean log2 observed data

% sig is the signal to be analyzed via an AR(p) model

% Represents oscillations in sig via p-order AR(p) model 
p=92; % model order (p is chosen arbitrarily and should be
% larger than 52 (season period))
a=armcov(sig,p);  % AR model estimation via the modified 
% covariance method
a=a(:); % makes it a column vector

e1=filter(a,1,sig);  % obtains the modeling error via inverse filtering
% Obs.: e1 will be used to obtain the initial conditions (IC) of the direct
% filter, during forecast
e=filter2(a,sig,'same'); % centered inverse filtering, i.e., vector 'e' is
% time-synchronous with sig
% 'e' is supposed to be white Gaussian noise
% Sequence 'e' will be used to estimate the typical variance of the model
% excitation used to generate the forecast realizations

eo=e;  % original centered excitation (model error) 

PP=2*52;  % twice the period of 52 EW
% Idea here: we need a forecast of PP-25 weeks ahead. So we will organize
% sequence 'e' into a matrix (column wise) with PP rows, so that the 
% excitation signals are all
% synchronized from EW 1 to EW 52 of the subsequent year. 
% We later trim-off
% the first 25 rows, so that the excitation signals cover 
% EW 26 of a given year to EW 52 of the
% subsequent year. We then calculate a mean excitation from the matrix
% columns.


E=buffer(e,PP); % organizes vector 'e' in a matrix with PP rows.
% Each column of E has an excitation vector related to 2 years (104 EWs)

% E will be used to obtain a set of artificial excitation signals 
% (white Gaussian noise) to the AR model to
% extrapolate the data (time-series forecasting)

E=E'; % transpose matrix E (now with PP columns).
% Each row of E has a excitation vector related to 2 years (104 EWs). 
% For a given year: EW 1 of a given year to EW 52 of the subsequent year


mE=mean(E);  % row vector of the mean excitation, i.e., mean 
% of the ensemble (7) of excitation signals of PP samples
stdE=std(E); % row vector of the standard deviation 
% of the ensemble (7) of excitation signals of PP samples

med_stdE=median(stdE); % median of vector stdE (informs 
% the standard deviation of 
% a typical white gaussian noise to be generated to excite 
% the direct AP(p) model

% Note: we could use a different standard deviation for each EW,
% as they are not about equal  
% across EW. Plot stdE to see the pattern. Nevertheless, we are using
% zero-mean excitation signals
% which are white Gaussian noise realizations with the
% same standard-deviation (med_stdE) across the EWs

stdE2=med_stdE.*ones(size(stdE)); % row vector (PP values) with
% the same value med_stdE

eek=mE(26:end);  % known mean excitation related to EW 26 of 2010 
% to EW 25 of 2022 
eek=eek(:);  % make it a column vector

[aux,zf]=filter(1,a,e1);  % obtains the initial conditions 'zf' of 
% the direct source-model filter

% Note: zf is the set of initial conditions for the direct AR(p) forecast
% filter, to which an artificial excitation will be fed.

MC=10000; % number of runs of a Monte Carlo simulation

% Now, we generate the set of realizations of artificial model
% excitation signals.
RE=randn(MC,PP);  % matrix with random realizations of white 
% Gaussian noise (zero-mean), unit variance 
RE=RE.*(1.0*stdE2); % adjust variance of the noise according to 
% that of the mean excitation 
% We keep a zero-mean realizations and only adjust the variance.

RE(:,1:25)=[];  % removes first 25 columns (since forecast starts 
% from EW 26 onward)

% Each column of matrix RE contains a realization of the excitation to be
% fed to the AR(p) model

CP_v=zeros(MC,lwea);  % blank matrix to store the forecast values of dengue cases

% The forecast below (with a mean measured excitation) is only for 
% comparison purposes
cases_prediction_d=filter(1,a,eek,(zf));  % forecast using the mean 
% excitation (measured)
cases_prediction_d=2.^(cases_prediction_d+mcc); % maps back to original
% scale

% Runs the Monte Carlo Simulation
for kk=1:MC 
    ee=RE(kk,:);  % one realization of the excitation 
    % (artificially generated)
    cases_prediction=filter(1,a,ee,(zf)); % direct filtering to obtain
    % extrapolated signal with case predicted values
    cases_prediction=2.^(cases_prediction+mcc); % maps back to the 
    % original scale
    CP_v(kk,:)=cases_prediction; % store forecast values
end

% Calculate statistics based on the set of forecast sequences

% Obtain approximations of the [50 80 90 95]% prediction intervals 
% (data driven)


set_prctile=[2.5 5 10 25 50 75 90 95 97.5]; % 2.5 to 97.5% percentiles
PP=prctile(CP_v,set_prctile); % calculates the percentiles and stores in PP
lower_95 = PP(1,:);  % 2.5% percentile 
lower_90 = PP(2,:);  % 5% percentile
lower_80 = PP(3,:);  % 10% percentile
lower_50 = PP(4,:);  % 25% percentile
pred = PP(5,:);  % 50% percentile - median prediction
upper_50 = PP(6,:); % 75% percentile
upper_80 = PP(7,:); % 90% percentile 
upper_90 = PP(8,:); % 95% percentile 
upper_95 = PP(9,:); % 97.5% percentile 

% Note: the forecast results (mean, upper and lower bounds) are then 
% filtered by an SSA (Singular Spectral Analysis) reconstruction filter.

L=20; % window length for the SSA filter
nsv=5; % number of selected eigenvalues (ordered)
[pred]=round(ssa_modPE(pred,L,nsv)); % filtered mean forecast
[lower_95]=round(ssa_modPE(lower_95,L,nsv));  % filtered 2.5% quartile
[lower_90]=round(ssa_modPE(lower_90,L,nsv));  % filtered 5% quartile
[lower_80]=round(ssa_modPE(lower_80,L,nsv));  % filtered 10% quartile
[lower_50]=round(ssa_modPE(lower_50,L,nsv));  % filtered 25% quartile
[upper_50]=round(ssa_modPE(upper_50,L,nsv));  % filtered 75% quartile
[upper_80]=round(ssa_modPE(upper_80,L,nsv));  % filtered 90% quartile
[upper_90]=round(ssa_modPE(upper_90,L,nsv));  % filtered 95% quartile
[upper_95]=round(ssa_modPE(upper_95,L,nsv));  % filtered 97.5% quartile


indf_ini=665; % time index of the EW 41 2022
indf_end=716; % time index of the EW 40 2023

% Writes a CSV file with the known and forecast data

date=S_dates(667:718);  % Sunday dates for validation 1


pred_range=indf_end+1-indf_ini; % forecast range
gapf=15;  % gap in samples from EW 26 2022 to EW  40 2022  
pred(1:gapf)=[]; pred=pred(1:pred_range);
lower_95(1:gapf)=[]; lower_95=lower_95(1:pred_range); % 2.5% quartile 
lower_90(1:gapf)=[]; lower_90=lower_90(1:pred_range); % 5% quartile 
lower_80(1:gapf)=[]; lower_80=lower_80(1:pred_range); % 10% quartile 
lower_50(1:gapf)=[]; lower_50=lower_50(1:pred_range); % 25% quartile 
upper_50(1:gapf)=[]; upper_50=upper_50(1:pred_range); % 75% quartile
upper_80(1:gapf)=[]; upper_80=upper_80(1:pred_range); % 90% quartile
upper_90(1:gapf)=[]; upper_90=upper_90(1:pred_range); % 95% quartile
upper_95(1:gapf)=[]; upper_95=upper_95(1:pred_range); % 97.5% quartile

state_code=M{1,2}*ones(size(pred));

T = table(date,pred,lower_50,upper_50,lower_80,upper_80,lower_90,...
    upper_90,lower_95,upper_95);

writetable(T,['..\spreadsheets\v1_ARp26_',num2str(state_code(1)),'.csv'],'Delimiter',',')

% Generate plots for each state and saves in PDF files separately

max75=max(upper_50);
max90=max(upper_80);
max95=max(upper_90);
max97p5=max(upper_95);
EW_index=indf_ini:indf_end;

cases=kcc(EW_index);  % known observed dengue cases

figure
subplot(221)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases
% from EW 41 of 2022 to EW 40 of 2023
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_50,'k','linewidth',2)
plot(EW_index,upper_50,'g','linewidth',2)
legend('observed','forecast','lower\_50','upper\_50')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 50% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max75])

subplot(222)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases
% from EW 41 of 2022 to EW 40 of 2023
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_80,'k','linewidth',2)
plot(EW_index,upper_80,'g','linewidth',2)
legend('observed','forecast','lower\_80','upper\_80')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 80% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max90])

subplot(223)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases
% from EW 41 of 2022 to EW 40 of 2023
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_90,'k','linewidth',2)
plot(EW_index,upper_90,'g','linewidth',2)
legend('observed','forecast','lower\_90','upper\_90')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 90% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max95])

subplot(224)
plot(EW_index,cases,'linewidth',2); % plot known sequence of cases
% from EW 41 of 2022 to EW 40 of 2023
hold on;
plot(EW_index,pred,'r','linewidth',2)
plot(EW_index,lower_95,'k','linewidth',2)
plot(EW_index,upper_95,'g','linewidth',2)
legend('observed','forecast','lower\_95','upper\_95')  
xlabel('Time (EW index)')
ylabel('Number of Cases')
title(['Forecast 95% -  ',UF])
axis([EW_index(1) EW_index(end) 0 2*max97p5])

print(['..\plots\v1_ARp26_',UF],'-dpdf')

close all

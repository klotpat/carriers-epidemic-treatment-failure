function   centraltendency
beta = 0.000775;
f = 0.2;
Time = 140;
P_ext_ana = 0.57280963;
[T_allend]= simdata(beta,f,Time,P_ext_ana);
%xlswrite('file name.xlsx',T_allend',sheet,'column'); in case for
%recording the simulated data
sortdata = sort(T_allend);
sortdata(find(sortdata >= Time))=[]; % remove maximum time 
mean_data = mean(sortdata)
%std_data = std(sortdata)
%Q2dat = median(sortdata)
%TF = isoutlier(sortdata,'quartiles'); % detect outliers by indexing
%outlier = sortdata(TF); % define outliers data
%oir = sum(outlier)/sum(sortdata) % comput influence ratio
%ts = timeseries(sortdata');
%tsiqr = iqr(ts) % interquartile
end







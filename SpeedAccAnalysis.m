function [] = SpeedAccAnalysis()

% Read the data
VarNames = {...
    'Var1'
    'ppantId'
    'taskDuration'
    'practiceK'
    'practiceMeanRt'
    'apoeHaplotype'
    'gender'
    'ethnicity'
    'ageBand'
    'bmi'
    'alcohol'
    'smoker'
    'pastSmoker'
    'ageBand2'
    'e4'
    'dose'
    'loN'
    'loK'
    'loP'
    'hiN'
    'hiK'
    'hiP'
    'loMeanRt'
    'loSdRt'
    'hiMeanRt'
    'hiSdRt'
    'loOutlier'
    'loPercAcc'
    'zLoAcc'
    'hiPercAcc'
    'zHiAcc'
    'diffScore'
    'zDiffScore'};
dataTable = readtable('Greebles_SummaryData.csv',...
    'PreserveVariableNames',true);
dataTable.Properties.VariableNames = VarNames;

%% Remove unused vars
dataTable.Var1 = [];
dataTable.ageBand2 = [];
dataTable.loPercAcc = [];
dataTable.hiPercAcc = [];

%% Extract the data
dataTable.apoeHaplotype = categorical(dataTable.apoeHaplotype);
dose = dummyvar(dataTable.apoeHaplotype)*[0;1;2];
acc = dataTable.hiP;
rt = dataTable.hiMeanRt ./ 1000;

%% Plot
figure;
hold on;
for iD = 0:2
    scatter(acc(dose==iD),rt(dose==iD),50,'LineWidth',2);
end
xlabel('pCorrect');
ylabel('RT (seconds)');
colormap copper;
legend(unique(dataTable.apoeHaplotype),'Location','northwest');

%%
n = nan(3,1);
r = nan(3,1);
ci = nan(3,2);
t = nan(3,1);
p = nan(3,1);
for ii = 1:3
    cD = ii-1;
    n(ii) = sum(dose==cD);
    
    % Functions to calcuate confidence bounds on r values
    df = n(ii) - 2;
    t2r = @(t) sqrt((t.^2)./((t.^2) + df)) .* sign(t);
    r2t = @(r) r ./ sqrt((1-(r.^2))/df);
    rCI = @(r) t2r(r2t(r) + tinv([0.025,0.975],df));
    
    r(ii) = corr(acc(dose==cD),rt(dose==cD));
    ci(ii,:) = rCI(r(ii));
    t(ii) = r2t(r(ii));
    p(ii) = fcdf(t(ii)^2,1,df,'upper');
end
resultsTable = table(n,r,ci,t,p,'RowNames',cellstr(unique(dataTable.apoeHaplotype)));
disp(resultsTable);
return
clear;
close all;
clc;
rng(0);

%% Read in the model estimates
Estimates = readtable('Pc-05c_Conditions.csv');
Estimates.Ambiguity = categorical(Estimates.Ambiguity);
Estimates.Genotype = categorical(Estimates.Genotype);
Estimates.Age = categorical(Estimates.Age);

%% Fiter for average genotype
Estimates = Estimates(Estimates.Genotype==categorical({'ave'}),:);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
%% Frequentist estimates
E = reshape(Estimates.Frequ_Est,4,2);
L = reshape(Estimates.Frequ_Low,4,2) - E;
H = E - reshape(Estimates.Frequ_Hig,4,2);
CI = cat(3,L,H);
[hBar,hErrorbar] = barwitherr(CI,E);
set(hErrorbar(1),'LineWidth',1.5);
set(hErrorbar(2),'LineWidth',1.5);
set(hBar(1),'LineWidth',1.5);
set(hBar(1),'FaceColor',[0.9,0.9,0.9]);
set(hBar(2),'LineWidth',1.5);
set(hBar(2),'FaceColor',[0.5,0.5,0.5]);
hold on;

%% Appearnece
xticklabels(unique(Estimates.Age));
xlabel('Age group','FontSize',18);
ax = gca;
ax.FontSize = 18;
ylim([0.8,1]); %ylim([0,12]);
yticks((8:0.5:10)./10);
set(gca,'YMinorTick','Off');
ylabel('Probability of a correct response','FontSize',18);
legend({'Low ambiguity','High ambiguity'},'FontSize',15,'Location','bestoutside');
title('Age by ambiguity interaction (all genotypes)','FontSize',18);
box off;
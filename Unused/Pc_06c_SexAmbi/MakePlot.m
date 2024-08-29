clear;
close all;
clc;
rng(0);

%% Read in the model estimates
Estimates = readtable('Pc-06c_Conditions.csv');
Estimates.Ambiguity = categorical(Estimates.Ambiguity);
Estimates.Genotype = categorical(Estimates.Genotype);
Estimates.Sex = categorical(Estimates.Sex);
Estimates.Age = categorical(Estimates.Age);

%% Fiter for average genotype and age
Estimates = Estimates(Estimates.Genotype==categorical({'ave'}),:);
Estimates = Estimates(Estimates.Age==categorical({'mean'}),:);

%% Switch female and male around
Estimates = Estimates([2;1;4;3],:);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
%% Frequentist estimates
E = reshape(Estimates.Frequ_Est,2,2);
L = reshape(Estimates.Frequ_Low,2,2) - E;
H = E - reshape(Estimates.Frequ_Hig,2,2);
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
xticklabels(unique(Estimates.Sex));
%xlabel('Age group','FontSize',18);
ax = gca;
ax.FontSize = 18;
ylim([0.8,1]); %ylim([0,12]);
yticks((8:0.5:10)./10);
set(gca,'YMinorTick','Off');
ylabel('Probability of a correct response','FontSize',18);
legend({'Low ambiguity','High ambiguity'},'FontSize',15,'Location','bestoutside');
title('Sex by ambiguity interaction (all genotypes)','FontSize',18);
box off;
clear;
close all;
clc;
rng(0);

%% Read in the model estimates
Estimates = readtable('Pc-05f_Conditions.csv');
Estimates.Ambiguity = categorical(Estimates.Ambiguity);
Estimates.Genotype = categorical(Estimates.Genotype);

%%
FH = figure('units','normalized','outerposition',[0 0 1 1]);
uGenotype = unique(Estimates.Genotype);

%% Frequentist estimates
E = reshape(Estimates.Frequ_pEst,3,2);
L = reshape(Estimates.Frequ_pLow,3,2) - E;
H = E - reshape(Estimates.Frequ_pHig,3,2);
CI = cat(3,L,H);
[hBar,hErrorbar] = barwitherr(CI,E);
set(hErrorbar(1),'LineWidth',1.5);
set(hErrorbar(2),'LineWidth',1.5);
set(hBar(1),'LineWidth',1.5);
set(hBar(2),'LineWidth',1.5);
set(hBar(1),'FaceColor',[0.2,0.4,0.8]);
set(hBar(2),'FaceColor',[0.8,0.4,0.2]);
hold on;

%% Add the subject means
RawData = GetRawData();
xOffset = 0.285;
xJitter = 0.09;

x = dummyvar(RawData.Genotype)*([1;2;3]+0.143);
rng(101);
r = rand(numel(x),1);
xx = [];
yy = [];
for iAmbi = 1:2
    if iAmbi == 1
        y = RawData.LowAmbiguity_accuracy./RawData.LowAmbiguity_n;
        
    else
        y = RawData.HighAmbiguity_accuracy./RawData.HighAmbiguity_n;
    end
    deltaX = (iAmbi-1).*(xOffset) - xOffset;
    deltaX = deltaX + r.*(xJitter*2) - xJitter;
    scatter(x + deltaX, y,15,'+',...
        'MarkerEdgeColor',[0.1,0.1,0.1],'MarkerEdgeAlpha',0.5,'LineWidth',1);
    xx = [xx,x+deltaX];
    yy = [yy,y];
end
for iPoint = 1:size(xx,1)
    plot(xx(iPoint,:),yy(iPoint,:),':','Color','k');
end

%% Appearnece
xticklabels(uGenotype);
xlabel('Genotype','FontSize',18);
ax = gca;
ax.FontSize = 18;
ylim([0.5,1]); %ylim([0,12]);
yticks((5:1:10)./10);
set(gca,'YMinorTick','Off');
ylabel('Probability of a correct response','FontSize',18);
legend({'Low ambiguity','High ambiguity'},'FontSize',15,'Location','bestoutside');
box off;
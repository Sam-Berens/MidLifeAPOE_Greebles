clear;
close all;
clc;
rng(0);

%% Read in the model estimates
Estimates = readtable('RTplot_Estimates.csv');
Estimates.Ambiguity = categorical(Estimates.Ambiguity);
Estimates.Genotype = categorical(Estimates.Genotype);
Estimates.Age = categorical(Estimates.Age);

%% Read in the raw data to plot
SelectData = readtable('RTplot_SelectData.csv');
SelectData.PID = categorical(SelectData.PID);
SelectData.Ambiguity = categorical(SelectData.Ambiguity);
SelectData.Genotype = categorical(SelectData.Genotype);
SelectData.Age = categorical(SelectData.Age);
SelectData.AgeN = double(SelectData.Age);
SelectData = [SelectData(:,1:4),SelectData(:,8),SelectData(:,5:7)];

%%
figure('units','normalized','outerposition',[0 0 1 1]);
uGenotype = unique(SelectData.Genotype);
xOffset = 0.14;
xJitter = 0.08;
FH = [NaN,NaN];
for iPlt = 1:2 % 1==LowAmbi, 2==HigAmbi
    FH(iPlt) = subplot(1,2,iPlt);
    
    %% Frequentist estimates
    sIdx = ((iPlt-1)*8) + 1;
    eIdx = sIdx + 7;
    E = reshape(Estimates.Frequ_Est(sIdx:eIdx),4,2);
    L = reshape(Estimates.Frequ_Low(sIdx:eIdx),4,2) - E;
    H = E - reshape(Estimates.Frequ_Hig(sIdx:eIdx),4,2);
    CI = cat(3,L,H);
    [hBar,hErrorbar] = barwitherr(CI,E);
    set(hErrorbar(1),'LineWidth',1.5);
    set(hErrorbar(2),'LineWidth',1.5);
    set(hBar(1),'LineWidth',1.5);
    set(hBar(2),'LineWidth',1.5);
    set(hBar(1),'FaceColor',[0.2,0.4,0.8]);
    set(hBar(2),'FaceColor',[0.8,0.4,0.2]);
    hold on;
    
    %% Raw data
    if iPlt == 1
        Sambi = SelectData.Ambiguity==categorical({'Low'});
    else
        Sambi = SelectData.Ambiguity==categorical({'High'});
    end
    for iGeno = 1:numel(uGenotype)
        Sgeno = SelectData.Genotype == uGenotype(iGeno);
        S = Sambi & Sgeno;
        deltaX = (iGeno-1).*(xOffset*2) - xOffset;
        deltaX = deltaX + rand(sum(S),1).*(xJitter*2) - xJitter;
        scatter(SelectData.AgeN(S) + deltaX, SelectData.RT(S),5,'+',...
            'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2,'LineWidth',1);
    end
    
    %% Appearnece
    xticklabels(unique(SelectData.Age));
    xlabel('Age group','FontSize',15);
    ax = gca; 
    ax.FontSize = 12; 
    ylim([0,40]); %ylim([0,12]);
    set(gca,'YScale','log');
    yticks([1,2,4,8,16,32]);
    set(gca,'YMinorTick','Off')
    if iPlt == 1
        title('Low ambiguity','FontSize',18);
        ylabel('Response time (seconds)','FontSize',15);
    else
        title('High ambiguity','FontSize',18);
        legend({'e33','e4+'},'FontSize',15);
    end
    box off;
end
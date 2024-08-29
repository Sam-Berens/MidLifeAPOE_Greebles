clear;
close all;
clc;
rng(0);

%% Read in the model estimates
Estimates = readtable('Rt-05l_Conditions.csv');
Estimates.Ambiguity = categorical(Estimates.Condition);
Estimates.Genotype = categorical(Estimates.Genotype);
Estimates.Age = categorical(Estimates.Age);

%% Read in the raw data to plot
SelectData = readtable('Data.csv');
SelectData = SelectData(SelectData.Correct==1,:);
SelectData.PID = categorical(SelectData.PID);
SelectData.Ambiguity = categorical(SelectData.Ambiguity);
SelectData.Genotype = categorical(SelectData.Genotype);
SelectData.Age = categorical(SelectData.Age);
SelectData.AgeN = double(SelectData.Age);
SelectData = [SelectData(:,1:2),SelectData(:,5),SelectData(:,4),SelectData(:,7:8),SelectData(:,3)];
SelectData.RT = SelectData.RT./1000;

%%
figure('units','normalized','outerposition',[0 0 1 1]);
uGenotype = unique(SelectData.Genotype);
xOffset = 0.22;
xJitter = 0.06;
FH = [NaN,NaN];
for iPlt = 1:2 % 1==LowAmbi, 2==HigAmbi
    FH(iPlt) = subplot(1,2,iPlt);
    
    %% Frequentist estimates
    sIdx = ((iPlt-1)*12) + 1;
    eIdx = sIdx + 11;
    E = reshape(Estimates.Frequ_Est(sIdx:eIdx),4,3);
    L = reshape(Estimates.Frequ_Low(sIdx:eIdx),4,3) - E;
    H = E - reshape(Estimates.Frequ_Hig(sIdx:eIdx),4,3);
    CI = cat(3,L,H);
    [hBar,hErrorbar] = barwitherr(CI,E);
    set(hErrorbar(1),'LineWidth',1.5);
    set(hErrorbar(2),'LineWidth',1.5);
    set(hErrorbar(3),'LineWidth',1.5);
    set(hBar(1),'LineWidth',1.5);
    set(hBar(2),'LineWidth',1.5);
    set(hBar(3),'LineWidth',1.5);
    set(hBar(1),'FaceColor',[0.2,0.4,0.8]);
    set(hBar(2),'FaceColor',[0.3,0.8,0.3]);
    set(hBar(3),'FaceColor',[0.8,0.4,0.2]);
    hold on;
    
    %% Raw data
    if iPlt == 1
        Sambi = SelectData.Ambiguity==categorical({'L'});
    else
        Sambi = SelectData.Ambiguity==categorical({'H'});
    end
    for iGeno = 1:numel(uGenotype)
        Sgeno = SelectData.Genotype == uGenotype(iGeno);
        S = Sambi & Sgeno;
        deltaX = (iGeno-1).*(xOffset) - xOffset;
        deltaX = deltaX + rand(sum(S),1).*(xJitter*2) - xJitter;
        scatter(SelectData.AgeN(S) + deltaX, SelectData.RT(S),5,'+',...
            'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2,'LineWidth',1);
    end
    
    %% Appearnece
    xticklabels(unique(SelectData.Age));
    xlabel('Age group','FontSize',18);
    ax = gca; 
    ax.FontSize = 18; 
    ylim([0,72]); %ylim([0,12]);
    set(gca,'YScale','log');
    yticks([1,2,4,8,16,32,64]);
    set(gca,'YMinorTick','Off')
    if iPlt == 1
        title('Low ambiguity','FontSize',18);
        ylabel('Response time (seconds)','FontSize',18);
    else
        title('High ambiguity','FontSize',18);
        legend({'e33','e34','e44'},'FontSize',15,'Location','bestoutside');
    end
    box off;
end
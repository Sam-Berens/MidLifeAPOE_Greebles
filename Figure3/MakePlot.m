clear;
close all;
clc;
rng(0);

%% Read in the model estimates
Estimates = readtable('Pc-05f_Conditions.csv');
Estimates.Ambiguity = categorical(Estimates.Condition);
Estimates.Genotype = categorical(Estimates.Genotype);
Estimates.Age = categorical(Estimates.Age);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
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
    
    %% Appearnece
    xticklabels(unique(Estimates.Age));
    xlabel('Age group','FontSize',18);
    ax = gca; 
    ax.FontSize = 18; 
    ylim([0.8,1]); %ylim([0,12]);
    yticks((8:0.5:10)./10);
    set(gca,'YMinorTick','Off')
    if iPlt == 1
        title('Low ambiguity','FontSize',18);
        ylabel('Probability of a correct response','FontSize',18);
    else
        title('High ambiguity','FontSize',18);
        legend({'e33','e34','e44'},'FontSize',15,'Location','bestoutside');
    end
    box off;
end
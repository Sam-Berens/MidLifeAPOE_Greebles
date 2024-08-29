clear;
close all;
clc;
rng(0);

%% Read in the model estimates
Estimates = readtable('Rt-06c_Conditions.csv');
Estimates.Ambiguity = categorical(Estimates.Ambiguity);
Estimates.Genotype = categorical(Estimates.Genotype);
Estimates.Age = categorical(Estimates.Age);
Estimates.Sex = categorical(Estimates.Sex);

%% Fiter condtions
Estimates = Estimates(Estimates.Genotype==categorical({'ave'}),:);
Estimates = Estimates(Estimates.Age~=categorical({'mean'}),:);

%% Place female first
Estimates = Estimates([5;6;7;8;1;2;3;4;13;14;15;16;9;10;11;12],:);

%% Read in the raw data to plot
Data = readtable('Data.csv');
Data.PID = categorical(Data.PID);
Data.Age = arrayfun(@(n)find(n==unique(Data.zAge)),Data.zAge);

%%
figure('units','normalized','outerposition',[0 0 1 1]);
uSex = unique(Estimates.Sex);
xOffset = 0.145;
xJitter = 0.07;
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
    set(hBar(1),'FaceColor',[1,110/255,199/255]);
    set(hBar(2),'FaceColor',[0,1,1]);
    hold on;
    
    %% Raw data
    if iPlt == 1
        Sambi = Data.hAmbi==0;
    else
        Sambi = Data.hAmbi==1;
    end
    for iSex = 1:numel(uSex)
        Ssex = Data.Female == 2-iSex;
        S = Sambi & Ssex;
        deltaX = (iSex-1).*(xOffset).*2 - xOffset;
        deltaX = deltaX + rand(sum(S),1).*(xJitter*2) - xJitter;
        scatter(Data.Age(S) + deltaX, Data.Rt(S),5,'+',...
            'MarkerEdgeColor','k','MarkerEdgeAlpha',0.2,'LineWidth',1);
    end
    
    %% Appearnece
    xticklabels(unique(Estimates.Age));
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
        legend({'Female','Male'},'FontSize',15,'Location','bestoutside');
    end
    box off;
end
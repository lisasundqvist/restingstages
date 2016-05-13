clear all
close all
clc

% Parameters and start values
seasons = 300; % number of seasons
stats = 1000; % number of runs to calculate statistics
startsize = 80; % number of allels starting blooms in the two different strategies
alleles1 = [1:20]; % possible alleles in population 1 and 3
alleles2 = [21:40]; % possible alleles in population 2 and 4
percentsediment= [0.05 0.25 0.75]; % percent allels picked from the sediment
k = [0.6 0.15 0.05]; % rate of decay in the exponetial curve from which the sediment alleles are drawn
history = [10 40 120]; % the number of seasons back in time from which alleles can be regenerated
migrationrate = [0.01 0.05 0.2]; % Migration rate
seasonssplot=1:seasons;
endbloomplot1 = zeros(startsize,seasons);
endbloomplot2 = zeros(startsize,seasons);
endbloomplot3 = zeros(startsize,seasons);
endbloomplot4 = zeros(startsize,seasons);
DifferentiationD1 = zeros(stats,seasons);
DifferentiationD2 = zeros(stats,seasons);
position = 1; % for positioning in figure

for mig = 1:3 %  migration rate is changed within this loop
    
    migrants = ceil(migrationrate(mig)*startsize); % migrationrate times startsize equals number of migrants
    
    for sed = 1:3 % sediment parameters percent regenerated, seasons back and k is changed within this loop
        
        startfrombloom = floor(startsize*(1-percentsediment(sed))); % number of allels starting from last years bloom
        startfromsediment = ceil(startsize*percentsediment(sed)); % number of alleles starting from the sediment
        sedimenttobloom1 = zeros(1,startfromsediment);
        sedimenttobloom2 = zeros(1,startfromsediment);
        
        % Calculations
        
        for statistic = 1:stats

            sediment1 = []; % empty sediment
            sediment2 = []; % empty sediment
            
            for ii = 1:history(sed)
                sediment1(ii,:) = randsample(alleles1,startsize,true); % building up a starting sediment
                sediment2(ii,:) = randsample(alleles2,startsize,true);
            end
            
            pop1 = randsample(alleles1,startsize,true); % starting population 1 with resting cells
            pop2 = randsample(alleles2,startsize,true); % starting population 2 with resting cells
            pop3 = randsample(alleles1,startsize,true); % starting population 3 without resting cells
            pop4 = randsample(alleles2,startsize,true); % starting population 4 without resting cells
            
            for strategi = 1:2; % test two strategiers, with and without resting cells
                
                for t = 1:seasons % number of seasons with one bloom in each
                    if strategi==1 % Strategi 1, forming resting cells
                        
                        % migration
                        migrants1 = randsample(pop1,migrants,true); % picks migrants from population 1
                        migrants2 = randsample(pop2,migrants,true); % picks migrants from population 2
                        pop1(randperm(numel(pop1),numel(migrants2))) = migrants2; % replaces alleles with alleles migrated from population 2
                        pop2(randperm(numel(pop2),numel(migrants1))) = migrants1; % replaces alleles with alleles migrated from population 1
                        
                        % save bloom
                        endbloom1 = pop1;
                        endbloom2 = pop2;
                        
                        % make sediment
                        sediment1 = [pop1;sediment1];
                        sediment1 = sediment1(1:history(sed),:); % keep only the needed depth
                        sediment2 = [pop2;sediment2];
                        sediment2 = sediment2(1:history(sed),:); % keep only the needed depth
                        
                        % set pop to zero
                        pop1 = zeros(1,startsize);
                        pop2 = zeros(1,startsize);
                        
                        % add allels from the privious bloom to the next bloom
                        pop1 = randsample(endbloom1,startfrombloom,true);
                        pop2 = randsample(endbloom2,startfrombloom,true);
                        
                        % add allels from the sediment to the next bloom
                        column1 = randi([1 startsize],1,startfromsediment); % picks columns
                        row1 = floor(-1/k(sed)*log(exp(-k(sed)*1)+rand(1,startfromsediment)*(exp(-k(sed)*(history(sed)+1))-exp(-k(sed)*1))));% picks rows with an exponentiallly decreasing function, from this page http://www.mathworks.com/matlabcentral/newsreader/view_thread/292852
                        column2 = randi([1 startsize],1,startfromsediment); % picks columns
                        row2 = floor(-1/k(sed)*log(exp(-k(sed)*1)+rand(1,startfromsediment)*(exp(-k(sed)*(history(sed)+1))-exp(-k(sed)*1))));% from this page http://www.mathworks.com/matlabcentral/newsreader/view_thread/292852
                        
                        for i = 1:startfromsediment % picks alleles from sediment according to positions generated above
                            sedimenttobloom1(i) = sediment1(row1(i),column1(i));
                            sedimenttobloom2(i) = sediment2(row2(i),column2(i));
                        end
                        
                        pop1(startfrombloom+1:startsize) = sedimenttobloom1; % adds allels to pop
                        pop2(startfrombloom+1:startsize) = sedimenttobloom2;
                        
                        % calculate differentation
                        totpop1 = [endbloom1;endbloom2]; % merges populations to calculate comparable allele frequencies for the two populations
                        countElements1=histc(endbloom1,unique(totpop1)); % counts alleles in pop 1
                        freqA=countElements1/numel(endbloom1); % calculates allele frequencies of pop 1
                        countElements2=histc(endbloom2,unique(totpop1)); % counts alleles in pop 2
                        freqB=countElements2/numel(endbloom2); % calculates allele frequencies of pop 2
                        H1cyst = 1-sum(freqA.^2); % calculates H1
                        H2cyst = 1-sum(freqB.^2); % calculates H2
                        Hscyst = (H1cyst+H2cyst)/2; % calculates Hs
                        Meanfreq1 = (freqA+freqB)/2;
                        Htcyst= 1-(sum(Meanfreq1.^2)); % calculates Ht
                        DifferentiationD1(statistic,t) = ((Htcyst-Hscyst)/(1-Hscyst))*(2/(2-1)); % calculates Jost's D, number of subpopulations = 2
                        
                    else % Strategi 2, not forming resting cells
                        
                        % migration
                        migrants3 = randsample(pop3,floor(migrants),true); % picks migrants from population 3
                        migrants4 = randsample(pop4,floor(migrants),true); % picks migrants from population 4
                        pop3(randperm(numel(pop3),numel(migrants4))) = migrants4; % replaces alleles with alleles migrated from population 4
                        pop4(randperm(numel(pop4),numel(migrants3))) = migrants3; % replaces alleles with alleles migrated from population 3
                        
                        % save bloom
                        endbloom3 = pop3;
                        endbloom4 = pop4;
                        
                        % set pop to zero
                        pop3 = zeros(1,startsize);
                        pop4 = zeros(1,startsize);
                        
                        % add allels from the privious bloom to the next bloom
                        pop3 = randsample(endbloom3,startsize,true);
                        pop4 = randsample(endbloom4,startsize,true);
                        
                        % calculate differentation
                        totpop2 = [endbloom3;endbloom4]; % merges populations to calculate comparable allele frequencies for the two populations
                        countElements3=histc(endbloom3,unique(totpop2)); % counts alleles in pop 3
                        freqC=countElements3/numel(endbloom3); % calculates allele frequencies of pop 3
                        countElements4=histc(endbloom4,unique(totpop2)); % counts alleles in pop 4
                        freqD=countElements4/numel(endbloom4); % calculates allele frequencies of pop 4
                        H1nocyst = 1-sum(freqC.^2); % calculates H1
                        H2nocyst = 1-sum(freqD.^2); % calculates H2
                        Hsnocyst = (H1nocyst+H2nocyst)/2; % calculates Hs
                        Meanfreq2 = (freqC+freqD)/2;
                        Htnocyst= 1-(sum(Meanfreq2.^2)); % calculates Ht
                        DifferentiationD2(statistic,t) = ((Htnocyst-Hsnocyst)/(1-Hsnocyst))*(2/(2-1)); % calculates Jost's D, number of subpopulations = 2
                    end
                    
                end
            end
            countdown = stats-statistic
        end
        
        %CI regular
        SEcyst = std(DifferentiationD1)/sqrt(length(DifferentiationD1(:,1))); % Standard Error
        SEnocyst = std(DifferentiationD2)/sqrt(length(DifferentiationD2(:,1))); % Standard Error
        CIerrorcyst = 1.9623*SEcyst;
        CIerrornocyst = 1.9623*SEnocyst;
        
%         %CI bootstrap
%         capable = @mean; % Bootstrap parameter
%         CIbootcyst = bootci(2000,capable,DifferentiationD1); % BCa confidence interval
%         CIbooterrorcyst = mean(DifferentiationD1)-CIbootcyst(1,:);
%         CIbootnocyst = bootci(2000,capable,DifferentiationD2); % BCa confidence interval
        
        figure(1)
        set(gcf,'Color','w')
        subplot(3,3,position)
        set(gca,'linewidth',2.0,'fontsize',14,'fontname','arial','fontweight','bold','color','w')
        hold on
        plot(seasonssplot,mean(DifferentiationD1),'b','LineWidth',2.0)
        plot(mean(DifferentiationD1)+CIerrorcyst,'b.','MarkerSize',6.0)
        plot(mean(DifferentiationD1)-CIerrorcyst,'b.','MarkerSize',6.0)
        plot(seasonssplot,mean(DifferentiationD2),'r','LineWidth',2.0)
        plot(mean(DifferentiationD2)+CIerrornocyst,'r.','MarkerSize',6.0)
        plot(mean(DifferentiationD2)-CIerrornocyst,'r.','MarkerSize',6.0)
        axis([0 seasons ,-0.1 1])
        %ylabel('Josts D')
        %xlabel('Years')
        
        position = position + 1;
    end
end
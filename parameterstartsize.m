clc 
clear all
close all

% Parameters and start values
stats = 1000; % number of runs to calculate statistics
percentsediment = 0.25; % percent allels picked from the sediment
k = 0.3; % rate of decay in the exponetial curve from which the sediment alleles are drawn
history = 20; % the number of seasons back in time from which alleles can be regenerated
numberofalleles = 20;
alleles = [1:numberofalleles]; % possible alleles
teststartsize = [10 40 80 120 160];% number of allels starting the two different strategies

for ii = 1:5 % startsize is changed within this loop
    
   % Evaluated parameters 
    startsize = teststartsize(ii);
    startfrombloom = floor(startsize*(1-percentsediment)); % number of allels starting from last years bloom
    startfromsediment = ceil(startsize*percentsediment); % number of alleles starting from the sediment
    sedimenttobloom = zeros(1,startfromsediment); % empty vector
    
    % Calculations
    
    for statistic = 1:stats
        
        % startvalues
        pop1 = randsample(alleles,startsize,true); % starting population with resting stages
        pop2 = randsample(alleles,startsize,true); % starting population without resting stages
        endbloom1 = pop1; % saving the population
        endbloom2 = pop2; % saving the population
        endbloomstart1=endbloom1;
        endbloomstart2=endbloom2;
        sediment = []; % predefining matrix
        t=1; % t calculats the number of times within the first while loop
        tt=1; % tt calculats the number of times within the second while loop
        for i = 1:history
        sediment(i,:) = randsample(alleles,startsize,true); % building up a starting sediment
        end
        
        while numel(unique(endbloom1))>0.1*numberofalleles % Strategi 1, forming resting cells
        
            % save bloom
            endbloom1 = pop1;
            
            % make sediment
            sediment = [pop1;sediment];
            sediment = sediment(1:history,:); % keep only the needed depth
            
            % set pop to zero
            pop1 = zeros(1,startsize);
            
            % add allels from the privious bloom to the next bloom
            pop1 = randsample(endbloom1,startfrombloom,true);
           
            % add allels from the sediment to the next bloom
            column = randi([1 startsize],1,startfromsediment); % picks columns
            row = floor(-1/k*log(exp(-k*1)+rand(1,startfromsediment)*(exp(-k*(history+1))-exp(-k*1)))); % picks rows with an exponentiallly decreasing function, from this page http://www.mathworks.com/matlabcentral/newsreader/view_thread/292852
            
            for i = 1:startfromsediment % picks alleles from sediment according to positions generated above
                sedimenttobloom1(i) = sediment(row(i),column(i));
            end
            
            pop1(startfrombloom+1:startsize) = sedimenttobloom1; % adds allels to pop
            t = t+1;% calculates the number of times whithin this while loop
            timeto90lostcyst(statistic,ii) = t; % collects the result number of seasons until 90% of original allels are lost
            
        end
        
        while numel(unique(endbloom2))>0.1*numberofalleles % Strategi 2, not forming resting cells
            
            % save bloom
            endbloom2 = pop2;

            % set pop to zero
            pop2 = zeros(1,startsize);
            
            % add allels from the privious bloom to the next bloom
            pop2 = randsample(endbloom2,startsize,true);
            tt = tt+1; % calculates the number of times whithin this while loop
            timeto90lostnocyst(statistic,ii) = tt; % collects the result number of seasons until 90% of original allels are lost
            
        end
    end
    ii % countdown
end

%CI regular
SEcyst = std(timeto90lostcyst)/sqrt(length(timeto90lostcyst(:,1)));               % Standard Error
SEnocyst = std(timeto90lostnocyst)/sqrt(length(timeto90lostnocyst(:,1)));         % Standard Error
CIerrorcyst = 1.9623*SEcyst;
CIerrornocyst = 1.9623*SEnocyst;

% %CI bootstrap
% capable = @mean;                                        % Bootstrap parameter
% CIbootcyst = bootci(2000,capable,timeto90lostcyst);            % BCa confidence interval
% CIbooterrorcyst = mean(timeto90lostcyst)-CIbootcyst(1,:);
% CIbootnocyst = bootci(2000,capable,timeto90lostnocyst);            % BCa confidence interval

figure(1)
set(gcf,'Color','w')
set(gca,'linewidth',2.0,'fontsize',14,'fontname','arial','fontweight','bold','color','w')
set(gca,'xtick',teststartsize)
hold on
errorbar(teststartsize,mean(timeto90lostcyst),CIerrorcyst,'bx','LineWidth',1.5)
errorbar(teststartsize,mean(timeto90lostnocyst),CIerrornocyst,'rx','LineWidth',1.5)
legend('With resting cells','Without resting cells')
xlabel('Start size')
ylabel('Seasons until 90% of alleles are lost')


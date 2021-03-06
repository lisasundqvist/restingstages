clear all
close all
clc

% Parameters
startsize = 80; % number of allels starting the two different strategies
stats = 1000; % number of runs to calculate statistics
numberofalleles = 20;
alleles = [1:numberofalleles]; % possible alleles
percentsediment = 0.25; % percent allels picked from the sediment
startfrombloom = floor(startsize*(1-percentsediment)); % number of allels starting from last years bloom
startfromsediment = ceil(startsize*percentsediment); % number of alleles starting from the sediment
sedimenttobloom = zeros(1,startfromsediment); % empty vector
testhistory = [0 10 20 40 80 120]; % seasons back that alleles can be picked from the sediment
testk = [0 0.6 0.3 0.15 0.075 0.05]; % adjusting the slope of the exponential function picking alleles from the sediment so that the function looks the same when history is changed

for ii = 1:6 % history and k is changed within this loop
    
    % Evaluated parameters
    history = testhistory(ii); % seasons back that alleles can be picked from the sedimen
    k = testk(ii); % adjusting the slope of the exponential function picking alleles from the sediment so that the function looks the same when history is changed
    
    % Calculations
    
    for statistic = 1:stats % number of simulation for statistic calculations
        
        % startvalues
        pop = randsample(alleles,startsize,true); % starting a population
        endbloom = pop; % saving the population
        sediment = []; % predefining matrix
        t=1; % t calculats the number of times within the while loop
        for i = 1:history % building up a starting sediment
            sediment(i,:) = randsample(alleles,startsize,true); % building up a starting sediment
        end
        
        while numel(unique(endbloom))>1
            
            % save bloom
            endbloom = pop;
            
            % make sediment
            sediment = [pop;sediment];
            sediment = sediment(1:history,:); % keep only the needed depth
            
            % set pop to zero
            pop = zeros(1,startsize);
            
            if history >=1 % when alleles are added from the sediment
                
                % add allels from the privious bloom to the next bloom
                pop = randsample(endbloom,startfrombloom,true);
                
                % add allels from the sediment to the next bloom
                column = randi([1 startsize],1,startfromsediment); % picks columns
                row = floor(-1/k*log(exp(-k*1)+rand(1,startfromsediment)*(exp(-k*(history+1))-exp(-k*1)))); % picks rows with an exponentiallly decreasing function, from this page http://www.mathworks.com/matlabcentral/newsreader/view_thread/292852
                
                for i = 1:startfromsediment % picks alleles from sediment according to positions generated above
                    sedimenttobloom1(i) = sediment(row(i),column(i));
                end
                pop(startfrombloom+1:startsize) = sedimenttobloom1; % adds allels to pop
                
            else
                pop = randsample(endbloom,startsize,true); % new pop from last seasons bloom when alleles are not added from the sediment
                
            end
            
            t = t+1;% calculates the number of times whithin the while loop
            timetofixation(statistic,ii) = t; % collects the result number of seasons until 90% of original allels are lost
            
        end
    end
    ii % countdown
end

% CI regular
SE = std(timetofixation)/sqrt(length(timetofixation(:,1)));  % Standard Error
CIerror = 1.9623*SE;
CI = mean(timetofixation)+CIerror;

% % CI bootstrap
% capable = @mean;                                       % Bootstrap parameter
% CIboot = bootci(2000,capable,timetofixation);            % BCa confidence interval
% CIbooterror = mean(timetofixation)-CIboot(1,:);

figure(1)
set(gcf,'Color','w')
set(gca,'linewidth',2.0,'fontsize',14,'fontname','arial','fontweight','bold','color','w')
hold on
set(gca,'xtick',testhistory)
axis([-5 130,0 4500])
errorbar(testhistory,mean(timetofixation),CIerror,'kx','LineWidth',2.0)
xlabel('Seasons back in sediment')
ylabel('Seasons until fixation')

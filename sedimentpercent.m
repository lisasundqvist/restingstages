clear all
close all
clc

% Parameters
startsize = 90; % number of allels starting the two different strategies
stats = 1000; % number of runs to calculate statistics
numberofalleles = 30;
alleles = [1:numberofalleles]; % possible alleles
k = 0.3; % rate of decay in the exponetial curve from which the sediment alleles are drawn
history = 20; % the number of seasons back in time from which alleles can be regenerated
testpercentsediment = [0 0.05 0.1 0.25 0.5 0.75]; % percent starting from the sediment

for ii = 1:6 % percent sediment is changed within this loop
    
    % Parameters and start values
    percentsediment = testpercentsediment(ii); % the result is evaluated for different percentages of allels picked from the sediment
    startfrombloom = floor(startsize*(1-percentsediment)); % number of allels starting from last seasons bloom
    startfromsediment = ceil(startsize*percentsediment); % number of alleles starting from the sediment
    sedimenttobloom = zeros(1,startfromsediment); % empty vector
    
    % Calculations
    for statistic = 1:stats % number of simulation for statistic calculations
        
        % startvalues
        pop = randsample(alleles,startsize,true); % starting a population
        endbloom = pop; % saving the population
        sediment = []; % predefining matrix
        t=1; % t calculats the number of times within the while loop
        for i = 1:history % building up a starting sediment
            sediment(i,:) = randsample(alleles,startsize,true); 
        end
        
        while numel(unique(endbloom))>=0.1*numberofalleles
            
            % save bloom
            endbloom = pop;
            
            % make sediment
            sediment = [pop;sediment];
            sediment = sediment(1:history,:); % keep only the needed depth
            
            % set pop to zero
            pop = zeros(1,startsize);
            
            % add allels from the privious bloom to the next bloom
            pop = randsample(endbloom,startfrombloom,true);
            
            if startfromsediment >=1 % when alleles are added from the sediment 
                
                % add allels from the sediment to the next bloom
                column = randi([1 startsize],1,startfromsediment); % picks columns
                row = floor(-1/k*log(exp(-k*1)+rand(1,startfromsediment)*(exp(-k*(history+1))-exp(-k*1)))); % picks rows with an exponentiallly decreasing function, from this page http://www.mathworks.com/matlabcentral/newsreader/view_thread/292852
                
                for i = 1:startfromsediment % picks alleles from sediment according to positions generated above
                    sedimenttobloom(i) = sediment(row(i),column(i));
                end
                
                pop(startfrombloom+1:startsize) = sedimenttobloom; % adds allels to pop
                
            else
                pop(startfrombloom+1:startsize) = randsample(endbloom,startfromsediment,true); % adds additional alleles from bloom when allels are not picked from the sediment 
            end
            
            t = t+1; % calculates the number of times whithin the while loop
            timeto90lost(statistic,ii) = t; % collects the result number of seasons until 90% of original allels are lost
            
        end
        
    end
    ii % countdown
end

%CI regular
SE = std(timeto90lost)/sqrt(length(timeto90lost(:,1)));  % Standard Error
CIerror = 1.9623*SE;
CI = mean(timeto90lost)+CIerror;

%CI bootstrap
capable = @mean;                                       % Bootstrap parameter
CIboot = bootci(2000,capable,timeto90lost);            % BCa confidence interval
CIbooterror = mean(timeto90lost)-CIboot(1,:);

testpercent=100*testpercentsediment;

figure(1)
set(gcf,'Color','w')
set(gca,'linewidth',2.0,'fontsize',14,'fontname','arial','fontweight','bold','color','w')
hold on
set(gca,'xtick',testpercent)
axis([-5 80,0 1500])
errorbar(testpercent,mean(timeto90lost),CIerror,'kx','LineWidth',2.0)
xlabel('% regenerated from sediment')
ylabel('Seasons until 90% of alleles are lost')



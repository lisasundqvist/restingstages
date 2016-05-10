clear all
close all
clc

% Parameters and start values
years = 800; % number of seasons
startsize = 90;% number of allels starting the two different strategies
percentsediment = 0.25; % percent allels picked from the sediment
startfrombloom = floor(startsize*(1-percentsediment)); % number of allels starting from last years bloom
startfromsediment = ceil(startsize*percentsediment); % number of alleles starting from the sediment
numberofalleles = 30;
alleles = [1:numberofalleles]; % possible alleles
k = 0.3; % rate of decay in the exponetial curve from which the sediment alleles are drawn
history = 20; % the number of seasons back in time from which alleles can be regenerated
for i = 1:history
   sediment(i,:) = randsample(alleles,startsize,true); % building up a starting sediment
end
pop1 = randsample(alleles,startsize,true); % starting population with resting stages
pop2 = randsample(alleles,startsize,true); % starting population without resting stages

% Predefinitions
yearsplot=1:years;
endbloomplot1 = zeros(startsize,years);
endbloomplot2 = zeros(startsize,years);
sedimenttobloom = zeros(1,startfromsediment);

% Calculations
    for strategi = 1:2; % tests two life history strategies (forming and not forming resting cells)
        for t = 1:years % number of seasons with one bloom in each
            if strategi==1 % forming resting cells
                
                % save bloom
                endbloom1 = pop1;
                endbloomplot1(:,t) = sort(endbloom1);
                
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
                
            else  % startegi 2, not forming resting cells
                
                % save bloom
                endbloom2 = pop2;
                endbloomplot2(:,t) = sort(endbloom2);
                                
                % set pop to zero
                pop2 = zeros(1,startsize);
                
                % add allels from the privious bloom to the next bloom
                pop2 = randsample(endbloom2,startsize,true);
                
            end
        end
    end

figure(1)
subplot(2,1,1)
set(gcf,'Color','w')
set(gca,'linewidth',2.0,'fontsize',14,'fontname','arial','fontweight','bold','color','w')
hold on
axis([0 800,0 90])
imagesc(endbloomplot1)
title('Resting cells')
ylabel('Number of alleles')
subplot(2,1,2)
set(gcf,'Color','w')
set(gca,'linewidth',2.0,'fontsize',14,'fontname','arial','fontweight','bold','color','w')
hold on
axis([0 800,0 90])
imagesc(endbloomplot2)
title('No resting cells')
xlabel('Seasons')
ylabel('Number of alleles')

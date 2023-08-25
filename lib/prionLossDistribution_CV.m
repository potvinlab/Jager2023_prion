function [lostTime,keptTime,lostCell,lostIndex,keptCell,keptIndex,rejectedCell, rejectedIndex,neverCell,neverIndex, frame_index, loss_curve] = prionLossDistribution_CV(dataCell,plots,oldPole,startTime)

% Classify cells that never had the prion (neverIndex), had it but lost it 
% (lostIndex), kept it for the whole trace (keptIndex), and puts the rest
% in rejected. Analysis is done only for the small or old pole aggregates,
% the other category going in the rejected. Also finds when the cells
% lost the prion (lostTime).
% 
% This version uses the lack of movement from spots to classify them as old
% pole or not.
%
% INPUT ARGUMENTS
% dataCell: data from catLineagesCell, with the images of the cells or
%   the spots already found included
% plots: level of plot to show (default 1)
%         1: plot histograms, CDF, residual lifetime
%         2:also show the classification for each trace
% oldPole: analyze old pole aggregate or small aggregates
%       0: small, 1: old pole
% startTime: start time for the analysis, in frames

if nargin ==1
   plots =1; 
   startTime = 1;
elseif nargin ==2
    startTime = 1;
        
end

if ~isfield(dataCell,'spotNo')
    %find spots if not already done
    disp('Finding spots...')
    dataCell = dataSpotFind(dataCell);
    disp('Done finding spots')
end
    


keptCell={};
keptIndex=[];

lostCell = {};
lostIndex= [];

rejectedCell = {};
rejectedIndex = [];

neverCell = {};
neverIndex= [];

lostTime = [];
keptTime = [];

stdallnonzeros = @(x) std(x(x>0),0,'all');
meanallnonzeros = @(x) mean(x(x>0),'all');

%cycle through traces, calling prion status

for i=1:length(dataCell.spotNo)
    ind = find(dataCell.spotNo{i});
    %only keep traces that start at beginning for now
    if dataCell.frameNo{i}(1)<20   && median(dataCell.YFPAvg{i}(ind)) > 10
        
        
        % 3 out of the first 15 frames have spot, start in prion
        
        if sum(dataCell.spotNo{i}(5:20)>=1)>=3 

            
            median_cv = median(dataCell.YFPCV{i}(ind));
            is_old_pole = median_cv > 1.3;
            % only track the prion type (small or old pole) that is passed
            % in argument
            if oldPole == is_old_pole
                
                %did they lose it
                loseAt= strfind(dataCell.spotNo{i}(12:end),[0 0 0 0 0 0 0 0])+11;
                
                if ~isempty(loseAt)
                    %do they have the prion at some time after loss? If so,
                    %reject trace
                    
                    prionAgain = sum(dataCell.spotNo{i}(loseAt(1):end) >=1);
                    
                    if prionAgain < 5
                        lostTime(length(lostTime)+1) = loseAt(1)+dataCell.frameNo{i}(1)-1;
                        lostCell{length(lostCell)+1}=dataCell.spotNo{i};
                        lostIndex(length(lostIndex)+1) = i;
                        
                        
                    else
                        
                        rejectedCell{length(rejectedCell)+1}=dataCell.spotNo{i};
                        rejectedIndex(length(rejectedIndex)+1) = i;
                        
                    end
                    
                    
                else
                    %kept the whole time
                    keptCell{length(keptCell)+1}=dataCell.spotNo{i};
                    keptIndex(length(keptIndex)+1) = i;
                    
                    keptTime(length(keptTime)+1) = length(dataCell.spotNo{i});
                    
                end
                
            end
            
            
            
        else %never had prion
            neverCell{length(neverCell)+1}=dataCell.spotNo{i};
            neverIndex(length(neverIndex)+1) = i;
            
        end
        
        
    else
        rejectedCell{length(rejectedCell)+1}=dataCell.spotNo{i};
        rejectedIndex(length(rejectedIndex)+1) = i;
        
        
    end
    
    
end

lostTime = lostTime-startTime +1;
keptTime = keptTime-startTime +1; 

%Calculate cumulative distribution
fractionAtLeast=zeros(1,max([keptTime lostTime]));
residualLifetime = fractionAtLeast;
for t=1:length(fractionAtLeast)
    %number of cells that kept the prion more than t, divided by the
    %number of cells we followed for more than t
    fractionAtLeast(t)=(sum(keptTime>t)+sum(lostTime>t))/(sum(keptTime>t)+sum(lostTime>1));

end

frame_index = (1:length(fractionAtLeast))-1 ;
loss_curve = fractionAtLeast;


if plots>=1
    
   figure;histogram(lostTime) 
    title('Time lost')
    
    figure;histogram(keptTime)
    title('Kept for at least');
    


    timeIndex = ((1:length(fractionAtLeast))-1)*8/60;

    figure;plot(timeIndex,fractionAtLeast);
    xlabel('Time (h)')
    ylabel('Fraction of cells with prions')

    
    figure;
    timeIndex = ((1:length(fractionAtLeast))-1);
    for t=1:length(residualLifetime)
        residualLifetime(t)=mean(lostTime(lostTime>t)-t);
        
    end
    plot(timeIndex, residualLifetime);
    title(['Residual Lifetime. ' inputname(1)])
    xlabel('Time (frames)')
    ylabel('Residual lifetime (frames)')
    
 


    if plots >=2

    plotPrionProperty(dataCell,'spotNo',length(fractionAtLeast))
    title('spotNo')
    
        plotPrionProperty(dataCell,'YFPAvg',length(fractionAtLeast))
    title('YFPAvg')
    
        plotPrionProperty(dataCell,'Area',length(fractionAtLeast))
    title('Area')
    
    
        plotPrionProperty(dataCell,'YFPPeak',length(fractionAtLeast))
    title('YFPPeak')
    
           plotPrionProperty(dataCell,'peakOverAvg',length(fractionAtLeast))
    title('PeakOverAvg')
    
    end
    
end




if plots >=3
    figure;hold on;
    for i=1:length(keptCell)
         scatter(dataCell.frameNo{keptIndex(i)},repmat(keptIndex(i),1,length(keptCell{i})),[],keptCell{i}>=1,'filled')
        
    end
    title('Cells that kept prion')
    
    figure;hold on;
    for i=1:length(lostCell)
         scatter(dataCell.frameNo{lostIndex(i)},repmat(lostIndex(i),1,length(lostCell{i})),[],lostCell{i}>=1,'filled')
        
    end
    title('Cells that lost prion')
    
        figure;hold on;
    for i=1:length(neverCell)
         scatter(dataCell.frameNo{neverIndex(i)},repmat(neverIndex(i),1,length(neverCell{i})),[],neverCell{i}>=1,'filled')
        
    end
    title('Cells that did not start with prion')
    
    
            figure;hold on;
    for i=1:length(rejectedCell)
         scatter(dataCell.frameNo{rejectedIndex(i)},repmat(rejectedIndex(i),1,length(rejectedCell{i})),[],rejectedCell{i}>=1,'filled')
        
    end
    title('Rejected cells')
    
    
end



    function plotPrionProperty(dataCell,field,maxTime)
        
        meanProp = zeros(1,maxTime);
        stdProp =meanProp;
        for time=1:maxTime
            cellIndices = [keptIndex( keptTime>=time) lostIndex(lostTime>=time)];
            
            if ~strcmp(field,'peakOverAvg')
                
                tempCell = dataCell.(field)(cellIndices);
                tempArray = zeros(1,length(tempCell));
                for j=1:length(tempArray)
                    tempArray(j)=tempCell{j}(time);
                    
                end
            else
                tempCellPeak = dataCell.YFPPeak(cellIndices);
                tempCellAvg = dataCell.YFPAvg(cellIndices);

                tempArray= zeros(1,length(tempCellPeak));
                
                
                for j=1:length(tempArray)
                    tempArray(j)=tempCellPeak{j}(time)./tempCellAvg{j}(time);
                    
                end
                
            end
            meanProp(time) = mean(tempArray);
            stdProp(time) = std(tempArray);
        end
        figure;errorbar(meanProp,stdProp);
    end


end





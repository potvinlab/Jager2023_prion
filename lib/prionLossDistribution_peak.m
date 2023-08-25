function [lostTime,keptTime,lostCell,lostIndex,keptCell,keptIndex,rejectedCell, rejectedIndex,neverCell,neverIndex, frame_index, loss_curve] = prionLossDistribution_peak(dataCell,plots,oldPole,startTime)

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


%cycle through traces, calling prion status

for i=1:length(dataCell.spotNo)
    ind = find(dataCell.spotNo{i});
    %only keep traces that start at beginning for now
    if dataCell.frameNo{i}(1)<20   && median(dataCell.YFPAvg{i}(ind)) > 10
        
        
        % 3 out of the first 15 frames have spot, start in prion
        
        if sum(dataCell.spotNo{i}(5:20)>=1)>=3 
            
            % if tracking old pole aggregates, select peak/avg >5 and vice
            % versa. Don't "reject" the traces that don't follow the
            % criteria
            if oldPole == (median(dataCell.YFPPeak{i}(ind)./dataCell.YFPAvg{i}(ind))> 5)
                
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
    
   
   
    timeIndex = ((1:length(fractionAtLeast))-1)*8;

    figure;plot(timeIndex,fractionAtLeast,  'color', [0.49,0.74,0.80]);%[0.49,0.74,0.80][0.7490   0.6863  0.8000]
    xlabel('Time (min)')
   ylabel('Fraction of cells in the prion state (1-CDF)')

    figure;
    timeIndex = ((1:length(fractionAtLeast))-1)*8;
    for t=1:length(residualLifetime)
        residualLifetime(t)=mean(lostTime(lostTime>t)-t)*8;
        
    end
    plot(timeIndex, residualLifetime);
    title(['Residual Lifetime. ' inputname(1)])
    xlabel('Time (min)')
    ylabel('Residual lifetime (min)')
    



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





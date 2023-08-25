function [alignedTraces] = alignSpotCenAtLossBeforeAfter(dataCell,cenNo,lossTime,lossIndex)
% [alignedTraces] =
% alignSpotCenAtLossBeforeAfter(dataCell,fieldToAlign,lossTime,lossIndex)
MAX_LENGTH=2001;

alignedTraces=zeros(length(lossIndex),MAX_LENGTH);
alignedTraces(:)=NaN;
zeroPosition = (MAX_LENGTH-1)/2 +1;
for i=1:length(lossIndex)
    traceIndices = (1:length(dataCell.spotCen{lossIndex(i)}))+zeroPosition-(lossTime(i)- dataCell.frameNo{lossIndex(i)}(1) +1) ;
    
    averageCenTrace = zeros(1,length(traceIndices));
    for j=1:length(traceIndices)
        centroids = dataCell.spotCen{lossIndex(i)}{j};
        if isempty(centroids)
            averageCenTrace(j) = NaN;
        else
            averageCenTrace(j) = mean(centroids(cenNo:2:end));
        end
            
    end
    alignedTraces(i,traceIndices) = averageCenTrace;
      
  
end


end


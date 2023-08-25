function [alignedTraces] = alignFieldAtLossBeforeAfter(dataCell,fieldToAlign,lossTime,lossIndex)
% [alignedTraces] =
% alignFieldAtLoss(dataCell,fieldToAlign,lossTime,lossIndex)
MAX_LENGTH=2001;

alignedTraces=zeros(length(lossIndex),MAX_LENGTH);
alignedTraces(:)=NaN;
zeroPosition = (MAX_LENGTH-1)/2 +1;
for i=1:length(lossIndex)
    traceIndices = (1:length(dataCell.(fieldToAlign){lossIndex(i)}))+zeroPosition-(lossTime(i)- dataCell.frameNo{lossIndex(i)}(1) +1) ;
    
    alignedTraces(i,traceIndices) = dataCell.(fieldToAlign){lossIndex(i)};
      
  
end


end


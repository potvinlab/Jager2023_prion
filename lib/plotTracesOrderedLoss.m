function [] = plotTracesOrderedLoss(data,lostTime,lostIndex)
%   plotTracesOrderedLoss(data,lostTime,lostIndex)
% 
%   Plot the traces in ordered by their time of loss colored by whether or
%   not they have the prion. 
%    data: contains the traces
%    lostTime: the lostTime of the prion from prionLossDistribution
%    lostIndex: the indices for the traces (from prionLossDistribution)

[lostTime, sortIndex] = sort(lostTime);
lostIndex = lostIndex(sortIndex);

tracesFrame = data.frameNo(lostIndex);
traces = data.spotNo(lostIndex);
maxLength=1;
for i=1:length(tracesFrame)
    if tracesFrame{i}(end)>maxLength
       maxLength= tracesFrame{i}(end);
        
    end
   
end

traceMatrices = zeros(length(traces),maxLength);
traceMatrices(:) = NaN;

for i=1:length(traces)
   traceMatrices(i,tracesFrame{i}) = traces{i};
    
end
yValues = 1:length(traces);
colorToPlot = traceMatrices>0;
colorToPlot = double(colorToPlot);
colorToPlot(isnan(traceMatrices))=NaN;

figure; 
h=imagesc(colorToPlot);
set(h,'alphadata',~isnan(colorToPlot))
xinFrames = get(h,'XData');
set(h,'Xdata',xinFrames*8);
colormap('copper')
xlabel('Time (min)')
ylabel('Cell index')
xlim('Auto')
end


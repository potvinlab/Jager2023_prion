function [cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    fraction_small] = align_trace_part_div_cen(data, minLossTime)

if nargin==1
    minLossTime = 30;
end
startTime = 1;
% Calculate partitioning error and cell cycle phase 
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig] =prionLossDistribution(data,0,1, startTime);

fraction_small = (length(lostTimeNormal)+length(keptIndexNormal))/(length(lostTimeNormal)+length(lostTimeBig)+length(keptIndexNormal)+length(keptIndexBig));
disp(['Number of cells (fraction): ' num2str(length(lostTimeNormal)+length(lostTimeBig)+length(keptIndexNormal)+length(keptIndexBig))])
neverIndex = neverIndexNormal;
lossIndex =  [lostIndexNormal lostIndexBig ];
lossTime =  [lostTimeNormal lostTimeBig ];
autofluo = 0; 
dataQ = prionPartErrorSignTot_keepOutliers(data,lossIndex,lossTime, neverIndex, autofluo);

fieldToAlign = 'cellPhase';
lossTime = lostTimeNormal;
lossIndex = lostIndexNormal;

if ~isfield(dataQ,'cellPhase')

    [ cellCyclePhase,fieldCell,cellCyclePhaseCell ] = cellcycleVsFieldv3( dataQ,'YFPAvg' );
    dataQ.cellPhase = cellCyclePhaseCell;

end

% histogram of cell cycle loss events 
[alignedTraces] = alignFieldAtLossBeforeAfter(dataQ,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

indexCellsNotAtDivision = alignedTraces(:,1001) >0;
alignedCellCycle = alignedTraces;

cell_cycle_at_loss = alignedTraces(:,1001);

disp(['Number of cells (loss): ' num2str(sum(~isnan(cell_cycle_at_loss)))])


% Part error in divisions

fieldToAlign = 'Q';
lossTime = lostTimeNormal;
lossIndex = lostIndexNormal;
%minLossTime = 30;


[alignedTraces] = alignFieldAtLossBeforeAfter(dataQ,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));
tracesInDivisions=nan(size(alignedTraces));
for i=1:size(alignedTraces,1)
    % Only plot divisions
    divisions = alignedCellCycle(i,:);
    beforeLoss = alignedTraces(i,divisions(1:1000)==0);
    atAndAfterLoss = alignedTraces(i,find(divisions(1001:2001)==0)+1000);
    tracesInDivisions(i,1000-length(beforeLoss)+1:1000) = beforeLoss;
    tracesInDivisions(i,1001:1001+length(atAndAfterLoss)-1) = atAndAfterLoss;
end


%at least 10 points
validTimePoints = sum(~isnan(tracesInDivisions),1)>10;

x_part = (-1000:1000) ; 
x_part = x_part(validTimePoints);
part_error_in_div = nanmedian(tracesInDivisions(:,validTimePoints),1);
part_error_in_div_ste = nanstd(tracesInDivisions(:,validTimePoints),[],1)./sqrt(sum(~isnan(tracesInDivisions(:,validTimePoints)),1));

abs_part_error_in_div = nanmedian(abs(tracesInDivisions(:,validTimePoints)),1);
abs_part_error_in_div_ste = nanstd(abs(tracesInDivisions(:,validTimePoints)),[],1)./sqrt(sum(~isnan(tracesInDivisions(:,validTimePoints)),1));
%At cell div
at_cell_div_ind = cell_cycle_at_loss==0;
validTimePoints = sum(~isnan(tracesInDivisions(at_cell_div_ind,:)),1)>10;

x_part_at_div = (-1000:1000) ; 
x_part_at_div = x_part_at_div(validTimePoints);
part_error_in_div_at_div = nanmedian(tracesInDivisions(at_cell_div_ind,validTimePoints),1);
part_error_in_div_at_div_ste = nanstd(tracesInDivisions(at_cell_div_ind,validTimePoints),[],1)./sqrt(sum(~isnan(tracesInDivisions(:,validTimePoints)),1));

%Not at cell div
at_cell_div_ind = cell_cycle_at_loss==0;
validTimePoints = sum(~isnan(tracesInDivisions(at_cell_div_ind,:)),1)>10;

disp(['Number of cells (part): ' num2str(max(sum(~isnan(tracesInDivisions(:,validTimePoints)),1)))])

x_part_not_at_div = (-1000:1000) ; 
x_part_not_at_div = x_part_not_at_div(validTimePoints);
part_error_in_div_not_at_div = nanmedian(tracesInDivisions(~at_cell_div_ind,validTimePoints),1);
part_error_in_div_not_at_div_ste = nanstd(tracesInDivisions(~at_cell_div_ind,validTimePoints),[],1)./sqrt(sum(~isnan(tracesInDivisions(:,validTimePoints)),1));

%Not at cell div

% Centroid positions


data=dataQ;
spotCen = 2;
lossTime = lostTimeNormal;
lossIndex = lostIndexNormal;
%minLossTime = 30;


[alignedTraces] = alignSpotCenAtLossBeforeAfter(data,spotCen,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

%at least 10 points
validTimePoints = sum(~isnan(alignedTraces),1)>10;
averageTrace = nanmean(alignedTraces,1);
time_min = (-1000:1000)*8 ;

time_min_part = time_min;
spot_pos = averageTrace;
spot_pos_ste = nanstd(alignedTraces,[],1)./sqrt(sum(~isnan(alignedTraces),1));

disp(['Number of cells :' num2str(max(sum(~isnan(alignedTraces),1)))])

fraction_at_div = sum(cell_cycle_at_loss==0)./length(cell_cycle_at_loss);
end

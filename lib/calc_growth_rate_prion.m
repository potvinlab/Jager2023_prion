function [doubling_time_no_prion, doubling_time_with_prion]=  calc_growth_rate_prion(dataCell_array)

no_prion_cat = [];
with_prion_cat = [];
doubling_time_no_prion = [];
doubling_time_with_prion = [];
n_time_points_no_prion = 0;
n_time_points_with_prion = 0;
for i = 1:length(dataCell_array)
    if isa(dataCell_array,"struct")
        dataCell = dataCell_array(i);
    else
        dataCell = dataCell_array{i};
    end
    startTime = 1;
    
    [no_prion, with_prion] = calc_single_growth_rate(dataCell);

    doubling_time_no_prion = [doubling_time_no_prion; median(no_prion)];
    n_time_points_no_prion = n_time_points_no_prion+  length(no_prion);

    doubling_time_with_prion = [doubling_time_with_prion; median(with_prion)];
    n_time_points_with_prion = n_time_points_with_prion+  length(with_prion);

end

    disp(['number of time points with prion: ' num2str(n_time_points_with_prion)])
    disp(['number of time points without prion: ' num2str(n_time_points_no_prion)])
end

function [no_prion, with_prion] = calc_single_growth_rate(data)
     startTime = 1;
    [lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
    [lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig] =prionLossDistribution(data,0,1, startTime);
    
    fieldToAlign = 'frameNo';
    lossIndex = neverIndexNormal;
    lossTime = ones(1,length(lossIndex))*50;
    minLossTime = 1;
    [alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));
    frameNoNever = alignedTraces;
    
    fieldToAlign = 'normGrowthRate';
    [alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));
    
    growthRateNever = alignedTraces(frameNoNever>30);
    
    
    fieldToAlign = 'frameNo';
    lossIndex = [lostIndexNormal lostIndexBig keptIndexNormal keptIndexBig];
    lossTime = [lostTimeNormal lostTimeBig keptTimeNormal keptTimeBig];
    [alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));
    
    frameNoLost = alignedTraces;
    
    fieldToAlign = 'normGrowthRate';
    [alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));
    
    delayAfterLoss = 20;
    
    afterLost = alignedTraces(:,1001+delayAfterLoss:2000);
    
    
    afterLost = afterLost(frameNoLost(:,1001+delayAfterLoss:2000)>30);
    
    beforeLost = alignedTraces(:,1:1000);
    beforeLost = beforeLost(frameNoLost(:,1:1000)>30);
    
    

    
    afterLost = 1./afterLost  *data.avgDivDuration*8;
    growthRateNever = 1./growthRateNever*data.avgDivDuration*8;
    beforeLost = 1./beforeLost *data.avgDivDuration*8;

    
    
    no_prion = [growthRateNever; afterLost];
    with_prion = beforeLost;


end

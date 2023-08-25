function [dist] = spot_change_dist(data,startTime)

[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
dist = [];
for i=1:length(lostTimeNormal)
    dist = cat(2,dist,diff(data.spotNo{lostIndexNormal(i)}(startTime:lostTimeNormal(i)+startTime-1)));

end

end

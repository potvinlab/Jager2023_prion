function [loss_rate_per_frame, loss_observed, frames_observed] = calc_loss_rate_per_frame(dataCell,startTime, want_old_pole)


[lostTime,keptTime] =prionLossDistribution(dataCell,0,want_old_pole,startTime);

frames_observed = sum(lostTime(lostTime>0)) + sum(keptTime(keptTime>0));
loss_observed = sum(lostTime>0)
loss_rate_per_frame = loss_observed / frames_observed
end

function [death_rate_never, death_rate_combined, death_rate_never_dist, death_rate_combined_dist] = calculate_death_rate_prion_combined(data, startTime)



% get rid of rejected cell in further analysis
cells_to_reject = find(isnan(data.death_call));
for i=1:length(cells_to_reject)
    data.frameNo{cells_to_reject(i)} = 1000;
end

[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig] =prionLossDistribution(data,0,1, startTime);



neverTime= neverIndexNormal;
for i=1:length(neverIndexNormal)
   neverTime(i) = data.frameNo{neverIndexNormal(i)}(end) - startTime +1;
end
t=1:max([lostTimeBig lostTimeNormal keptTimeNormal keptTimeBig neverTime ]) ;

survival_never = nan(1,length(t));
survival_normal = survival_never;
survival_big = survival_never;
survival_mixed = survival_never;
number_normal = survival_never;
number_big = survival_never;
number_never = survival_never;

never_death = data.death_time(neverIndexNormal); %-startTime???
small_lost_death = data.death_time(lostIndexNormal);
small_kept_death = data.death_time(keptIndexNormal);
big_lost_death = data.death_time(lostIndexBig);
big_kept_death = data.death_time(keptIndexBig);

% Death rate per frame
death_rate_never = sum(never_death<inf) / sum(neverTime);
sum(never_death<inf) 
sum(neverTime)*8/60
death_rate_small = (sum(small_lost_death<lostTimeNormal) + sum(small_kept_death<inf)) / (sum(lostTimeNormal)+sum(keptTimeNormal));

death_rate_big = (sum(big_lost_death<lostTimeBig) + sum(big_kept_death<inf)) / (sum(lostTimeBig)+sum(keptTimeBig));
(sum(small_lost_death<lostTimeNormal) + sum(small_kept_death<inf) +sum(big_lost_death<lostTimeBig) + sum(big_kept_death<inf))
(sum(lostTimeBig)+sum(keptTimeBig) + sum(lostTimeNormal)+sum(keptTimeNormal))*8/60
death_rate_combined = (sum(small_lost_death<lostTimeNormal) + sum(small_kept_death<inf) +sum(big_lost_death<lostTimeBig) + sum(big_kept_death<inf)) /(sum(lostTimeBig)+sum(keptTimeBig) + sum(lostTimeNormal)+sum(keptTimeNormal));
N_BOOTSTRAP = 1000;
death_rate_never_dist = zeros(1,N_BOOTSTRAP);
death_rate_small_dist = death_rate_never_dist;
death_rate_big_dist = death_rate_never_dist;
death_rate_combined_dist = death_rate_never_dist;
for i =1:N_BOOTSTRAP
    random_indices = randi(length(never_death),[1 length(never_death)]);
    death_rate_never_dist(i) = sum(never_death(random_indices)<inf) / sum(neverTime(random_indices));

    random_indices = randi(length(lostTimeNormal),[1 length(lostTimeNormal)]);
    random_indices_kept = randi(length(small_kept_death),[1 length(small_kept_death)]);
    death_rate_small_dist(i) = (sum(small_lost_death(random_indices)<lostTimeNormal(random_indices)) + sum(small_kept_death(random_indices_kept)<inf)) / (sum(lostTimeNormal(random_indices))+sum(keptTimeNormal(random_indices_kept)));

    random_indices_big = randi(length(lostTimeBig),[1 length(lostTimeBig)]);
    random_indices_kept_big = randi(length(keptTimeBig),[1 length(keptTimeBig)]);
    death_rate_big_dist(i) = (sum(big_lost_death(random_indices_big)<lostTimeBig(random_indices_big)) + sum(big_kept_death(random_indices_kept_big)<inf)) / (sum(lostTimeBig(random_indices_big))+sum(keptTimeBig(random_indices_kept_big)));
    death_rate_combined_dist(i) = (sum(big_lost_death(random_indices_big)<lostTimeBig(random_indices_big)) + sum(big_kept_death(random_indices_kept_big)<inf) + sum(small_lost_death(random_indices)<lostTimeNormal(random_indices)) + sum(small_kept_death(random_indices_kept)<inf)) / (sum(lostTimeNormal(random_indices))+sum(keptTimeNormal(random_indices_kept))+sum(lostTimeBig(random_indices_big))+sum(keptTimeBig(random_indices_kept_big)));

end

end

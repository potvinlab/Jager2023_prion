function [data] = identify_death_time(data)

GROWTH_RATE_THRESHOLD = 0.2;

number_traces = length(data.normGrowthRate);
death_time = zeros(1,number_traces);
death_time(:) = inf;
death_call = zeros(1,number_traces); % 0 alive at the end, 1 dead, NaN reject trace


for i=1:number_traces
    growth_rate_trace = data.normGrowthRate{i};
    thresholded_growth = growth_rate_trace > GROWTH_RATE_THRESHOLD;

    %did the cell growth at all
    if find(thresholded_growth,1,'first')
        % Did the cell die
        if mean(data.normGrowthRate{i}(end-5:end)) < GROWTH_RATE_THRESHOLD
            [death_index] = strfind(thresholded_growth, [1 0]);
            if length(death_index)==1
                death_call(i) = 1;
                death_time(i) = death_index;
            else
                death_call(i) = NaN;
                death_time(i) = NaN;
            end
        else
            death_call(i) = 0;
        end
    else
        death_call(i) = NaN;
        death_time(i) = NaN;
    end


end

data.death_call = death_call;
data.death_time = death_time;
end

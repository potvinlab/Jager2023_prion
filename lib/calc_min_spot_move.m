function [min_move] = calc_min_spot_move(spotCenTrace)
%   [min_move] = calc_min_spot_move(spotCenTrace)
%   Calculate the minimum eucledian distance between all spots from frame i
%   and all spots from frame i+1
%   spotCenTrace: cell containing the centroid of the spots over time
%   min_move: vector of length(spotCenTrace)-1 containing the min possible
%   movements
min_move = zeros(1,length(spotCenTrace)-1);
for time=1:(length(spotCenTrace)-1)
    if ~isempty(spotCenTrace{time}) && ~isempty(spotCenTrace{time+1})
        min_move_temp = 100000;
        for spot1=1:length(spotCenTrace{time})/2
            for spot2=1:length(spotCenTrace{time+1})/2
                min_move_temp = min(min_move_temp, euc_dist(spotCenTrace{time}((spot1-1)*2+1:(spot1-1)*2+2),spotCenTrace{time+1}((spot2-1)*2+1:(spot2-1)*2+2)));
            end
        end
        min_move(time) = min_move_temp;
    else % if there is a frame with missing spots, assign a large displacement
        min_move(time) = 5;
                       
    end
    
end

    function dist = euc_dist(pos1,pos2)
        dist = sqrt( (pos1(1)-pos2(1)).^2 + (pos1(2)-pos2(2)).^2);
        
    end

end

function [old_agg_pos,displacement] = old_pole_agg_position(cell_spot_cen)

old_agg_pos =[];
for t =1:length(cell_spot_cen)
    if length(cell_spot_cen{t})>=2
        %skip time points with no spots
        [~, old_pole_spot_index] = min(cell_spot_cen{t}(2:2:length(cell_spot_cen{t})));
        old_agg_pos=cat(2,old_agg_pos, cell_spot_cen{t}((old_pole_spot_index-1)*2+1:(old_pole_spot_index)*2));
    end
end

displacement = sqrt(diff(old_agg_pos(2,:)).^2 + diff(old_agg_pos(2,:)).^2);



end

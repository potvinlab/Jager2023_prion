function [dataCell] = ensure_divTime_in_absolute_time(dataCell)

in_abs_time = 1;

try
    for i = 1:length(dataCell.frameNo)
        if dataCell.frameNo{i}(1)>1


            for div = dataCell.divTime{i}
                div_index = find(dataCell.frameNo{i}==div,1,"first");
                assert(~isempty(div_index),"Division could not be found at trace " + num2str(i));
                area_increase= dataCell.Area{i}(div_index-1)>dataCell.Area{i}(div_index);
                assert(area_increase, "Area increase at division at trace "+ num2str(i))
            end


        end


    end
    disp("divTime in absolute time.")
catch err_msg
    disp(err_msg.message)
    in_abs_time =0;
    disp("divTime not in absolute time. Adjusting divTime.")
    for i = 1:length(dataCell.frameNo)
        if dataCell.frameNo{i}(1)>1
            %only necessary to adjust these
            dataCell.divTime{i} = dataCell.divTime{i} + dataCell.frameNo{i}(1) -1; 

        end
    end
end
    
if ~isfield(dataCell, 'prog') && isfield(dataCell, 'sisterYFPTot')
   dataCell.prog.frameNo = dataCell.frameNo;
   dataCell.prog.YFPTot = dataCell.sisterYFPTot;
   dataCell.prog.Area = dataCell.sisterArea;
   % Could copy everything else here too, but only need these partitioning
   % errors
end

if ~isfield(dataCell,'cellPhase')

    [ cellCyclePhase,fieldCell,cellCyclePhaseCell ] = cellcycleVsFieldv3( dataCell,'YFPAvg' );
    dataCell.cellPhase = cellCyclePhaseCell;

end

end

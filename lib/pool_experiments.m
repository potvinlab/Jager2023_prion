function [pool_dataCell] = pool_experiments(dataCell_array,start_time_array)



for i = 1:length(dataCell_array)
    if isa(dataCell_array,"struct")
        dataCell = dataCell_array(i);
    else
        dataCell = dataCell_array{i};
    end
    startTime = start_time_array(i);
    
    for j = 1:length(dataCell.frameNo)
        dataCell.frameNo{j} = dataCell.frameNo{j} - startTime +1;
        
    end
    if i==1
        pool_dataCell = dataCell;
        field_cell = fieldnames(dataCell);
    else
        for j =1:length(field_cell)
            if isfield(dataCell, field_cell{j})
                
                if field_cell{j}=="prog"
                    prog_field = fieldnames(pool_dataCell.prog);
                    for k=1:length(prog_field)
                        if isfield(dataCell.prog, prog_field{k})
                            pool_dataCell.prog.(prog_field{k}) = [pool_dataCell.prog.(prog_field{k}) dataCell.prog.(prog_field{k})];
                        end
                    end
                else
                    pool_dataCell.(field_cell{j}) = [pool_dataCell.(field_cell{j}) dataCell.(field_cell{j})];
                end
            end
        end
    end

    

end

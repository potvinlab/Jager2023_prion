function [varargout] = bootstrap_it(n_times,argout_index,function_handle,varargin)

varargout=cell(1,nargout);

output_cell = cell(n_times, nargout);
min_length = ones(1,nargout)*10000;
for i=1:n_times
    out_cell = cell(1,max(argout_index));
    data = varargin{1};
    random_indices = randi(length(data.Area),[1 length(data.Area)]);
    randomized_data = randomize_datacell(data,random_indices); 
    [out_cell{:}] = function_handle(randomized_data,varargin{2:end});
    output_cell(i,:) = out_cell(argout_index);
    for j=1:nargout
        min_length(j)=min(min_length(j),length(output_cell{i,j}));
    end
end

for i=1:n_times
    for j=1:nargout
    output_cell{i,j}=output_cell{i,j}(1:min_length(j));

    end
end

for j=1:nargout
    output_mat = cell2mat(output_cell(:,j));
    varargout{j}=std(output_mat,[],1);
end



end

function [rand_data] = randomize_datacell(data,random_indices)
    n_cells = length(data.Area);
    fields = fieldnames(data);
    for i=1:length(fields) 
        if length(data.(fields{i})) == n_cells
            rand_data.(fields{i}) = data.(fields{i})(random_indices);
        elseif strcmp(fields{i},'prog')
            rand_data.prog=randomize_datacell(data.prog,random_indices);
        else
            rand_data.(fields{i})=data.(fields{i});

        end
    end
end

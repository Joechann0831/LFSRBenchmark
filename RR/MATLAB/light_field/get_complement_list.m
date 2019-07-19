function [lf_names_train, datasets_train] = get_complement_list(lf_names_test, ...
    dataset_foldername)

% Get the list of folders contained in the dataset
list = get_folder_list(dataset_foldername);
m = 1;
for i = 1:size(list,1)
    % Derive the dataset to consider
    dataset = list{i};
    
    % Determine the foldername
    foldername = [dataset_foldername,list{i},'/'];
    
    if strcmp(dataset,'EPFL')  
        % Determine the light fields list
        list_file = get_file_list(foldername,'mat');
    elseif strcmp(dataset,'INRIA')
        list_file = get_file_list(foldername,'MAT');
    elseif strcmp(dataset,'HCI') || strcmp(dataset,'STANFORD')
        list_file = get_folder_list(foldername);
    end
    
    for j = 1:size(list_file,1)
        if strcmp(dataset,'HCI')||strcmp(dataset,'STANFORD')
            lf_name = list_file{j};
        else
            % Extract the light field name
            lf_name = list_file{j}(1:end-4);
        end
            
        flag = false;
        for k = 1:size(lf_names_test,2)
            if strcmp(lf_names_test{k},lf_name)
                flag = true;
                break;
            end
        end
            
        if flag == false
            lf_names_train{m} = lf_name;
            datasets_train{m} = dataset;
            m = m + 1;
        end
    end
end





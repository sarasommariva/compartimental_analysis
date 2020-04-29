clc
clear
close all

%% Step 1. Define general parameters
folder_data = './data';
data_path = fullfile(folder_data, ...
    'BRAIN_TOPI_per_analisi_compartimentale.xlsx');

%% Step 2. Read data from file excel

% 2.1. Time
time = xlsread(data_path, 'B5:B31');

% 2.2. First group
aux_group_diab = xlsread(data_path, 'C5:BN32');
aux_group_control = xlsread(data_path, 'BQ5:DT32');


%% Step 3. Organize data into matlab structure
% 3.1. Size of each group
n_mice_group_diab = size(aux_group_diab, 2)/4;
n_mice_group_control = size(aux_group_control, 2)/4;

if n_mice_group_diab ~= floor(n_mice_group_diab)
    error('Something went wrong when reading data..')
end

if n_mice_group_control ~= floor(n_mice_group_control)
    error('Something went wrong when reading data..')
end


% 3.2. Data from each group
%       - First group
mice_group_diab(n_mice_group_diab) = struct();
for im = 1:n_mice_group_diab
    
    idx_end = zeros(1, 4);
    for ic = 1:4
        aux_idx_end = find(isnan(aux_group_diab(:, 4*(im-1)+ic)), 1) - 1;
        if isempty(aux_idx_end) 
            idx_end(ic) = size(aux_group_diab, 1); 
        else
            idx_end(ic) = aux_idx_end;
        end
    end
    
    mice_group_diab(im).name = sprintf('mouse_diabetes_%d', im);
    mice_group_diab(im).if = aux_group_diab(1:idx_end(1), 4*im-3);
    mice_group_diab(im).wb = aux_group_diab(1:idx_end(2), 4*im-2);
    mice_group_diab(im).ant = aux_group_diab(1:idx_end(3), 4*im-1);
    mice_group_diab(im).post = aux_group_diab(1:idx_end(4), 4*im);
    
end

%       - Second group
mice_group_control(n_mice_group_control) = struct();
for im = 1:n_mice_group_control

    idx_end = zeros(1, 4);
    for ic = 1:4
        aux_idx_end = find(isnan(aux_group_control(:, 4*(im-1)+ic)), 1) - 1;
        if isempty(aux_idx_end) 
            idx_end(ic) = size(aux_group_control, 1); 
        else
            idx_end(ic) = aux_idx_end;
        end
    end

    mice_group_control(im).name = sprintf('mouse_%d', im);
    mice_group_control(im).if = aux_group_control(1:idx_end(1), 4*im-3);
    mice_group_control(im).wb = aux_group_control(1:idx_end(2), 4*im-2);
    mice_group_control(im).ant = aux_group_control(1:idx_end(3), 4*im-1);
    mice_group_control(im).post = aux_group_control(1:idx_end(4), 4*im);
    
end

%% Step 4. Save
save(fullfile(folder_data, 'diab_vs_control'), ...
    'mice_group_diab', 'mice_group_control', 'time')

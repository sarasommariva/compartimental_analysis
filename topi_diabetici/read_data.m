clc
clear
close all

folder_data = './data';
data_path = fullfile(folder_data, ...
    'BRAIN_TOPI_per_analisi_compartimentale.xlsx');

%% Step 1. Read data from file excel

% 1.1. Time
time = xlsread(data_path, 'B5:B31');

% 1.2. First group
aux_group_1 = xlsread(data_path, 'C5:BN32');
aux_group_2 = xlsread(data_path, 'BQ5:DT32');


%% Step 2. Organize data into matlab structure
% 2.1. Size of each group
n_mice_group_1 = size(aux_group_1, 2)/4;
n_mice_group_2 = size(aux_group_2, 2)/4;

if n_mice_group_1 ~= floor(n_mice_group_1)
    error('Something went wrong when reading data..')
end

if n_mice_group_2 ~= floor(n_mice_group_2)
    error('Something went wrong when reading data..')
end


% 2.2. Data from each group
%       - First group
mice_group_1(n_mice_group_1) = struct();
for im = 1:n_mice_group_1
    
    idx_end_if = find(isnan(aux_group_1(:, 4*im-3)), 1) - 1;
    if isempty(idx_end_if); idx_end_if = size(aux_group_1, 1); end
    idx_end_wb = find(isnan(aux_group_1(:, 4*im-2)), 1) - 1;
    if isempty(idx_end_wb); idx_end_wb = size(aux_group_1, 1); end
    idx_end_ant = find(isnan(aux_group_1(:, 4*im-1)), 1) - 1;
    if isempty(idx_end_ant); idx_end_ant = size(aux_group_1, 1); end
    idx_end_post = find(isnan(aux_group_1(:, 4*im)), 1) - 1;
    if isempty(idx_end_post); idx_end_post = size(aux_group_1, 1); end

    
    mice_group_1(im).name = sprintf('mouse_%d', im);
    mice_group_1(im).if = aux_group_1(1:idx_end_if, 4*im-3);
    mice_group_1(im).wb = aux_group_1(1:idx_end_wb, 4*im-2);
    mice_group_1(im).ant = aux_group_1(1:idx_end_ant, 4*im-1);
    mice_group_1(im).post = aux_group_1(1:idx_end_post, 4*im);
    
end

%       - Second group
mice_group_2(n_mice_group_2) = struct();
for im = 1:n_mice_group_2
    
    idx_end_if = find(isnan(aux_group_2(:, 4*im-3)), 1) - 1;
    if isempty(idx_end_if); idx_end_if = size(aux_group_2, 1); end
    idx_end_wb = find(isnan(aux_group_2(:, 4*im-2)), 1) - 1;
    if isempty(idx_end_wb); idx_end_wb = size(aux_group_2, 1); end
    idx_end_ant = find(isnan(aux_group_2(:, 4*im-1)), 1) - 1;
    if isempty(idx_end_ant); idx_end_ant = size(aux_group_2, 1); end
    idx_end_post = find(isnan(aux_group_2(:, 4*im)), 1) - 1;
    if isempty(idx_end_post); idx_end_post = size(aux_group_2, 1); end

    
    mice_group_2(im).name = sprintf('mouse_%d', im);
    mice_group_2(im).if = aux_group_2(1:idx_end_if, 4*im-3);
    mice_group_2(im).wb = aux_group_2(1:idx_end_wb, 4*im-2);
    mice_group_2(im).ant = aux_group_2(1:idx_end_ant, 4*im-1);
    mice_group_2(im).post = aux_group_2(1:idx_end_post, 4*im);
    
end

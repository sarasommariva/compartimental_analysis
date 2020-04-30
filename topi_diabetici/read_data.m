clc
clear
close all

%% Step 1. Define general parameters
folder_data = './data';
data_path = fullfile(folder_data, ...
    'BRAIN_TOPI_per_analisi_compartimentale.xlsx');

folder_results = './results';

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
save(fullfile(folder_results, 'diab_vs_control'), ...
    'mice_group_diab', 'mice_group_control', 'time')

%% Step 5. Preprocessing data
mice_group_diab_or = mice_group_diab;
mice_group_control_or = mice_group_control;

%%      5.a. Remove negative data
%   - First group
for im = 1:n_mice_group_diab
    mice_group_diab(im).if(mice_group_diab(im).if<0) = 0;
    mice_group_diab(im).wb(mice_group_diab(im).wb<0) = 0;
    mice_group_diab(im).ant(mice_group_diab(im).ant<0) = 0;
    mice_group_diab(im).post(mice_group_diab(im).post<0) = 0;
end
%   - Second group
for im = 1:n_mice_group_control
    mice_group_control(im).if(mice_group_control(im).if<0) = 0;
    mice_group_control(im).wb(mice_group_control(im).wb<0) = 0;
    mice_group_control(im).ant(mice_group_control(im).ant<0) = 0;
    mice_group_control(im).post(mice_group_control(im).post<0) = 0;
end

%%     5.b. Save
save(fullfile(folder_results, 'diab_vs_control_preproc'), ...
    'mice_group_diab', 'mice_group_control', 'time')


%% Step 6. Check on preprocessing
%%   6.1. Number negative values set to zero
fields = {'if', 'wb', 'ant', 'post'};
for ii = 1:numel(fields)
    % Group 1
    n_ist_neg_diab.(fields{ii}) = zeros(n_mice_group_diab, 1);
    max_negval_diab.(fields{ii}) = zeros(n_mice_group_diab, 1);
    for im = 1:n_mice_group_diab
        idx_neg = find(mice_group_diab(im).(fields{ii})==0);
        n_ist_neg_diab.(fields{ii})(im) = numel(idx_neg);
        if ~isempty(idx_neg)
        max_negval_diab.(fields{ii})(im) = ...
            max(abs(mice_group_diab_or(im).(fields{ii})(idx_neg)));
        end
    end
    % Group 2
    n_ist_neg_control.(fields{ii}) = zeros(n_mice_group_control, 1);
    max_negval_control.(fields{ii}) = zeros(n_mice_group_control, 1);
    for im = 1:n_mice_group_control
        idx_neg = find(mice_group_control(im).(fields{ii})==0);
        n_ist_neg_control.(fields{ii})(im) = numel(idx_neg);
        if ~isempty(idx_neg)
        max_negval_control.(fields{ii})(im) = ...
            max(abs(mice_group_control_or(im).(fields{ii})(idx_neg)));
        end
    end
end

% Plot
fc_negval = figure('units', 'normalized', 'outerposition',[0 0 1 1]);
subplot(2, 2, 1)
hold on
for ii = 1:numel(fields)
    plot(1:n_mice_group_diab, n_ist_neg_diab.(fields{ii}), 'o-', 'Linewidth', 2, ...
        'Displayname', fields{ii})
end
grid on
lgd = legend('show');
set(lgd, 'Location', 'Best', 'Fontsize', 15)

subplot(2, 2, 2)
hold on
for ii = 1:numel(fields)
    plot(1:n_mice_group_diab, max_negval_diab.(fields{ii}), 'o-', 'Linewidth', 2, ...
        'Displayname', fields{ii})
end
grid on
lgd = legend('show');
set(lgd, 'Location', 'Best', 'Fontsize', 15)

subplot(2, 2, 3)
hold on
for ii = 1:numel(fields)
    plot(1:n_mice_group_control, n_ist_neg_control.(fields{ii}), 'o-', 'Linewidth', 2, ...
        'Displayname', fields{ii})
end
grid on
lgd = legend('show');
set(lgd, 'Location', 'Best', 'Fontsize', 15)

subplot(2, 2, 4)
hold on
for ii = 1:numel(fields)
    plot(1:n_mice_group_control, max_negval_control.(fields{ii}), 'o-', 'Linewidth', 2, ...
        'Displayname', fields{ii})
end
grid on
lgd = legend('show');
set(lgd, 'Location', 'Best', 'Fontsize', 15)


%% 6.b

%% 6.c. Data mouse by mouse
fm = figure('units', 'normalized', 'outerposition',[0 0 1 1]);
for im = 1:n_mice_group_diab+n_mice_group_control
    if im <= n_mice_group_diab
        mouse_or = mice_group_diab_or(im);
        mouse = mice_group_diab(im);
        title_im = sprintf('Diab mouse %d', im);
    else
        mouse_or = mice_group_control_or(im-n_mice_group_diab);
        mouse = mice_group_control(im-n_mice_group_diab);
        title_im = sprintf('Control mouse %d', im - n_mice_group_diab);
    end
    
    for ii = 1:numel(fields)
        subplot(2, 2, ii)
        hold off
        plot(mouse_or.(fields{ii}), 'k', 'Linewidth', 2, ...
            'Displayname', sprintf('%s or', fields{ii}))
        hold on
        plot(mouse.(fields{ii}), 'r--', 'Linewidth', 2, ...
            'Displayname', sprintf('%s preproc', fields{ii}))
        if ii == 1
            title(title_im)
        end
    end
    pause
    
end






clc
clear
close all

%% Step 1. Define general parameters
folder_data = './data';
data_path = fullfile(folder_data, ...
    'BRAIN_TOPI_per_analisi_compartimentale.xlsx');

folder_results = './results';
folder_figures = './figures';

%% Step 2. Read data from file excel

% 2.1. Time
time = readmatrix(data_path, 'Range','B5:B31');

% 2.2. First group
aux_group_diab = [readmatrix(data_path, 'Range', 'C5:BN31'), ... 
                  readmatrix(data_path, 'Range', 'S36:BF62'), ...
                  readmatrix(data_path, 'Range', 'BK36:BN62')];
aux_group_control = [readmatrix(data_path, 'Range', 'BQ5:DT31'), ...
                     readmatrix(data_path, 'Range', 'CG36:CN62'), ... 
                     readmatrix(data_path, 'Range', 'DA36:DT62')];
    % Note: data-point outside the defined time interval are rejected.

%% Step 3. Organize data into matlab structure
% 3.1. Convert time from sec to min
time = time/60;

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

%%      5.b. Manually remove problematic time-points
%   - Mouse 15 control, posterior cortex, seems to miss one time-point
mice_group_control(15).post = [0; mice_group_control(15).post(1:end-1)];

time_prob_diab = zeros(n_mice_group_diab, 1);
time_prob_control = zeros(n_mice_group_control, 1);

time_prob_diab(13) = 16;
time_prob_control(1) = 16;
time_prob_control(6) = 16;
time_prob_control(8) = 16;
time_prob_control(9) = 16;
time_prob_control(11) = 21;
time_prob_control(15) = 16;

fields = {'if', 'wb', 'ant', 'post'};
for im = 1:n_mice_group_diab
    time_prob = time_prob_diab(im);
    if time_prob > 0
        for ii = 1:numel(fields)
            aux_data = mice_group_diab(im).(fields{ii}); 
            aux_data(time_prob) = (aux_data(time_prob+1)+aux_data(time_prob-1))/2;
            mice_group_diab(im).(fields{ii}) = aux_data(1:idx_end);
        end
    end
end

for im = 1:n_mice_group_control
    time_prob = time_prob_control(im);
    if time_prob > 0
        for ii = 1:numel(fields)
            aux_data = mice_group_control(im).(fields{ii}); 
            aux_data(time_prob) = (aux_data(time_prob+1)+aux_data(time_prob-1))/2;
            mice_group_control(im).(fields{ii}) = aux_data(1:idx_end);
        end
    end
end

%%      5.c. Manually define last time point
idx_end_diab = numel(time)*ones(n_mice_group_diab, 1);
idx_end_control = numel(time)*ones(n_mice_group_control, 1);

idx_end_diab(2) = 25;
idx_end_diab(3) = 25;
idx_end_diab(4) = 26;
idx_end_diab(13) = 26;
idx_end_diab(14) = 24;
idx_end_diab(25) = 26;
idx_end_diab(27) = 26;

idx_end_control(2) = 24;
idx_end_control(6) = 26;
idx_end_control(7) = 25;
idx_end_control(14) = 25;

for im = 1:n_mice_group_diab
    for ii = 1:numel(fields)
        idx_end = min(idx_end_diab(im), numel(mice_group_diab(im).(fields{ii})));
        mice_group_diab(im).(fields{ii}) = ...
            mice_group_diab(im).(fields{ii})(1:idx_end);
    end
end

for im = 1:n_mice_group_control
    for ii = 1:numel(fields)
        idx_end = min(idx_end_control(im), numel(mice_group_control(im).(fields{ii})));
        mice_group_control(im).(fields{ii}) = ...
            mice_group_control(im).(fields{ii})(1:idx_end);
    end
end


%%     5.b. Save
save(fullfile(folder_results, 'diab_vs_control_preproc'), ...
    'mice_group_diab', 'mice_group_control', 'time')


%% Step 6. Check on preprocessing
%   I Check for:
%       (1) The number of negative values set to zero
%       (2) The number of analysed time-point for each condition
fields = {'if', 'wb', 'ant', 'post'};
for ii = 1:numel(fields)
    % Group 1
    n_ist_neg_diab.(fields{ii}) = zeros(n_mice_group_diab, 1);
    max_negval_diab.(fields{ii}) = zeros(n_mice_group_diab, 1);
    n_ist_diab.(fields{ii}) = zeros(n_mice_group_diab, 1);
    for im = 1:n_mice_group_diab
        idx_neg = find(mice_group_diab(im).(fields{ii})==0);
        n_ist_neg_diab.(fields{ii})(im) = numel(idx_neg);
        if ~isempty(idx_neg)
        max_negval_diab.(fields{ii})(im) = ...
            max(abs(mice_group_diab_or(im).(fields{ii})(idx_neg)));
        end
        n_ist_diab.(fields{ii})(im) = ...
            numel(mice_group_diab(im).(fields{ii}));
    end
    
    % Group 2
    n_ist_neg_control.(fields{ii}) = zeros(n_mice_group_control, 1);
    max_negval_control.(fields{ii}) = zeros(n_mice_group_control, 1);
    n_ist_control.(fields{ii}) = zeros(n_mice_group_control, 1);
    for im = 1:n_mice_group_control
        idx_neg = find(mice_group_control(im).(fields{ii})==0);
        n_ist_neg_control.(fields{ii})(im) = numel(idx_neg);
        if ~isempty(idx_neg)
        max_negval_control.(fields{ii})(im) = ...
            max(abs(mice_group_control_or(im).(fields{ii})(idx_neg)));
        end
        n_ist_control.(fields{ii})(im) = ...
            numel(mice_group_control(im).(fields{ii}));
    end
end

%% Plots

%% P1. Number of negative values set to zero
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

%% P2 . Number of time points for each condition
fc_nist = figure('units', 'normalized', 'outerposition',[0 0 0.5 1]);

subplot(2, 1, 1)
hold on
for ii = 1:numel(fields)
    plot(1:n_mice_group_diab, n_ist_diab.(fields{ii}), 'o-', 'Linewidth', 2, ...
        'Displayname', fields{ii})
end
grid on
lgd = legend('show');
set(lgd, 'Location', 'Best', 'Fontsize', 15)

subplot(2, 1, 2)
hold on
for ii = 1:numel(fields)
    plot(1:n_mice_group_control, n_ist_control.(fields{ii}), 'o-', 'Linewidth', 2, ...
        'Displayname', fields{ii})
end
grid on
lgd = legend('show');
set(lgd, 'Location', 'Best', 'Fontsize', 15)



%% 6.b

%% 6.c. Data mouse by mouse
fields = {'if', 'wb', 'ant', 'post'};
for im = 1:n_mice_group_diab+n_mice_group_control
    
    fm = figure('units', 'normalized', 'outerposition',[0 0 1 1]);
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
        plot(mouse_or.(fields{ii}), 'k', 'Linewidth', 4, ...
            'Displayname', sprintf('%s or', fields{ii}), 'Displayname', 'Original')
        hold on
        plot(mouse.(fields{ii}), 'g', 'Linewidth', 3, ...
            'Displayname', sprintf('%s preproc', fields{ii}), 'Displayname', 'Corrected')
        lg = legend('show');
        if strcmp(fields{ii}, 'if')
            ylabel('C_b(t) [kBq/mL]', 'Fontsize', 20)
            set(lg, 'Fontsize', 20, 'Location', 'NorthEast')
        else
            ylabel(sprintf('C_T(t) [kBq/mL] %s', fields{ii}), 'Fontsize', 20)
            set(lg, 'Fontsize', 20, 'Location', 'SouthEast')
        end 
        xlim([0, numel(time)])
        if ii == 3 || ii == 4
            xlabel('Time t [min]', 'Fontsize', 20)
        end
    end
    title(mouse.name, 'Fontsize', 20, 'interpreter', 'none')
    
    saveas(fm, fullfile(folder_figures, ...
                    sprintf('data_corrected_%s.png', mouse.name)))
                
    pause(0.1)
    
    close(fm)
    
end






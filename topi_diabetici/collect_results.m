clc
clear
close all

%% Step 1. Define general parameter
% 1.a. Path
folder_func = fullfile('..', 'functions');
folder_results = fullfile('.', 'results');
folder_figures = fullfile('.', 'figures');

path_data = fullfile(folder_results, 'diab_vs_control_preproc.mat');
aux_path_result = 'result_%s_temp.mat';

% 1.b. Parameters for the analysis
Vb = 0.04;
Vi = 0;
v = 0.17;

t_0 = 0;
C_0_skf = [0, 0];
C_0_bcm = [0, 0, 0];

fields = {'ant', 'post'};
name_k_skf = {'k1', 'k2', 'k3', 'k4'}; n_k_skf = numel(name_k_skf);
name_k_bcm = {'k1', 'k2', 'k3', 'k5', 'k6'}; n_k_bcm = numel(name_k_bcm);

%% Step 2. Load
%   2.a. Data
addpath(folder_func)
load(path_data, 'mice_group_diab', 'mice_group_control', 'time')
time_or = time'; clear time 

idx_an_diab = [1:numel(mice_group_diab)];
idx_an_control = [1:numel(mice_group_control)];

%% Step 3. Store results
%   3.1. Initialize
n_an_diab = numel(idx_an_diab);
n_an_control = numel(idx_an_control);

for ii = 1:numel(fields)
    skf_diab.(fields{ii}).mean = zeros(n_an_diab, n_k_skf);
    skf_diab.(fields{ii}).std = zeros(n_an_diab, n_k_skf);
    skf_diab.(fields{ii}).opt = zeros(n_an_diab, n_k_skf);
    skf_diab.(fields{ii}).relerr_opt = zeros(n_an_diab, 1);
        
    skf_control.(fields{ii}).mean = zeros(n_an_control, n_k_skf);
    skf_control.(fields{ii}).std = zeros(n_an_control, n_k_skf);
    skf_control.(fields{ii}).opt = zeros(n_an_control, n_k_skf);
    skf_control.(fields{ii}).relerr_opt = zeros(n_an_control, 1);
    
    bcm_diab.(fields{ii}).mean = zeros(n_an_diab, n_k_bcm);
    bcm_diab.(fields{ii}).std = zeros(n_an_diab, n_k_bcm);
    bcm_diab.(fields{ii}).opt = zeros(n_an_diab, n_k_bcm);
    bcm_control.(fields{ii}).relerr_opt = zeros(n_an_diab, 1);
        
    bcm_control.(fields{ii}).mean = zeros(n_an_control, n_k_bcm);
    bcm_control.(fields{ii}).std = zeros(n_an_control, n_k_bcm);
    bcm_control.(fields{ii}).opt = zeros(n_an_control, n_k_bcm);
    bcm_control.(fields{ii}).relerr_opt = zeros(n_an_control, 1);
end

%   3.2. Results
for im = 1:n_an_diab+n_an_control
    
    if im <= n_an_diab
        aux_im = idx_an_diab(im);
        mouse = mice_group_diab(aux_im);
    else
        aux_im = idx_an_control(im-n_an_diab);
        mouse = mice_group_control(aux_im);
    end
    
    load(fullfile(folder_results, sprintf(aux_path_result, mouse.name)), 'ris');

    %  - Input function 
    c_if = mouse.if;
    time = time_or(1:numel(c_if));
    c_if = @(tt)(interp1([0 time],[0 c_if'], tt,'linear',0));

    
    for ii = 1:numel(fields)
        
%%      Skf model
        for jj = 1:n_k_skf
            
            aux_name_k = sprintf('%sx_skf', name_k_skf{jj});
            
            if im <= n_an_diab
    %  - Mean and std of k over ripetition
            skf_diab.(fields{ii}).mean(im, jj) = ...
                mean(ris.(fields{ii}).(aux_name_k));
            skf_diab.(fields{ii}).std(im, jj) = ...
                std(ris.(fields{ii}).(aux_name_k));
    %   - Optimal k 
            [~, idx_skf] = min(ris.(fields{ii}).relerr_skf);
            skf_diab.(fields{ii}).opt(im, jj) = ...
                ris.(fields{ii}).(aux_name_k)(idx_skf);
            skf_diab.(fields{ii}).relerr_opt(im) = ...
                ris.(fields{ii}).relerr_skf(idx_skf);
            else
    %  - Mean and std of k over ripetition
            skf_control.(fields{ii}).mean(im-n_an_diab, jj) = ...
                mean(ris.(fields{ii}).(aux_name_k));
            skf_control.(fields{ii}).std(im-n_an_diab, jj) = ...
                std(ris.(fields{ii}).(aux_name_k));
    %   - Optimal k 
            [~, idx_skf] = min(ris.(fields{ii}).relerr_skf);
            skf_control.(fields{ii}).opt(im-n_an_diab, jj) = ...
                ris.(fields{ii}).(aux_name_k)(idx_skf);
            skf_control.(fields{ii}).relerr_opt(im-n_an_diab) = ...
                ris.(fields{ii}).relerr_skf(idx_skf);
            end
        end
        
%%      Bcm model
        cond_idx = find(ris.(fields{ii}).k5x_bcm>=0.1);
        [~, aux_idx_bcm] = min(ris.(fields{ii}).relerr_bcm(cond_idx));
        idx_bcm = cond_idx(aux_idx_bcm);

        for jj = 1:numel(name_k_bcm)
            aux_name_k = sprintf('%sx_bcm', name_k_bcm{jj});
            if im <= n_an_diab
    %  - Mean and std of k over ripetition
            bcm_diab.(fields{ii}).mean(im, jj) = ...
                mean(ris.(fields{ii}).(aux_name_k));
            bcm_diab.(fields{ii}).std(im, jj) = ...
                std(ris.(fields{ii}).(aux_name_k));
    %   - Optimal k 
            bcm_diab.(fields{ii}).opt(im, jj) = ...
                ris.(fields{ii}).(aux_name_k)(idx_bcm);
            bcm_diab.(fields{ii}).relerr_opt(im) = ...
                ris.(fields{ii}).relerr_bcm(idx_bcm);
            else
    %  - Mean and std of k over ripetition
            bcm_control.(fields{ii}).mean(im-n_an_diab, jj) = ...
                mean(ris.(fields{ii}).(aux_name_k));
            bcm_control.(fields{ii}).std(im-n_an_diab, jj) = ...
                std(ris.(fields{ii}).(aux_name_k));
    %   - Optimal k 
            bcm_control.(fields{ii}).opt(im-n_an_diab, jj) = ...
                ris.(fields{ii}).(aux_name_k)(idx_bcm);
            bcm_control.(fields{ii}).relerr_opt(im-n_an_diab) = ...
                ris.(fields{ii}).relerr_bcm(idx_bcm);
            end  
        end
        
    end
    
%%  Step 4. Check whether reconstructions are ok
    y_max = max(max(mouse.(fields{1})), max(mouse.(fields{2})));
    f_ris = figure('units', 'normalized', 'outerposition',[0 0 1 1]);
    for ii = 1:numel(fields)
        % True valie
        c_tot = mouse.(fields{ii});
        % Reconstructed (average k over ripetition)
        if im <= n_an_diab
            aux_k_s = skf_diab.(fields{ii}).mean(im, :);
            aux_k_b = bcm_diab.(fields{ii}).mean(im, :);
        else
            aux_k_s = skf_control.(fields{ii}).mean(im-n_an_diab, :);
            aux_k_b = bcm_control.(fields{ii}).mean(im-n_an_diab, :);
        end
        ct_rec_skf_mean = forward_Skf(c_if, Vb, Vi, time, t_0, C_0_skf, ...
                     aux_k_s(1), aux_k_s(2), aux_k_s(3), aux_k_s(4));
        ct_rec_bcm_mean = forward_BCM(c_if, Vb, Vi, v, time, t_0, C_0_bcm, ...
                     aux_k_b(1), aux_k_b(2), aux_k_b(3), aux_k_b(4), aux_k_b(5));
        
        % Reconstructed (optimel k over ripetition)
        if im <= n_an_diab
            aux_k_s = skf_diab.(fields{ii}).opt(im, :);
            aux_k_b = bcm_diab.(fields{ii}).opt(im, :);
        else
            aux_k_s = skf_control.(fields{ii}).opt(im-n_an_diab, :);
            aux_k_b = bcm_control.(fields{ii}).opt(im-n_an_diab, :);
        end
        [ct_rec_skf_opt, c_rec_skf_opt.(fields{ii})] = ...
                    forward_Skf(c_if, Vb, Vi, time, t_0, C_0_skf, ...
                     aux_k_s(1), aux_k_s(2), aux_k_s(3), aux_k_s(4));
        [ct_rec_bcm_opt, c_rec_bcm_opt.(fields{ii})] = ...
                    forward_BCM(c_if, Vb, Vi, v, time, t_0, C_0_bcm, ...
                     aux_k_b(1), aux_k_b(2), aux_k_b(3), aux_k_b(4), aux_k_b(5));
        
        subplot(2, 2, 2*ii-1)
        hold on
        plot(time, ris.(fields{ii}).c_totx_skf, 'Color', [255, 204, 203]/255, ...
            'Linewidth', 2, 'HandleVisibility','off')
        plot(time, c_tot, 'k', 'Linewidth', 2, 'Displayname', 'True')
        plot(time, ct_rec_skf_mean, 'r--', 'Linewidth', 2, ...
            'Displayname', 'Rec mean')
        plot(time, ct_rec_skf_opt, 'r', 'Linewidth', 2, ...
            'Displayname', 'Rec opt')
        ylim([0, y_max])
        ylabel(sprintf('C_T %s', fields{ii}), 'Fontsize', 18)
        if ii == 1
            title(sprintf('%s - 2 comp', mouse.name), 'Fontsize', 18)
        end
        lgd = legend('show');
        set(lgd, 'Fontsize', 18, 'Location', 'SouthEast')

        subplot(2, 2, 2*ii)
        hold on
        plot(time, ris.(fields{ii}).c_totx_bcm, 'Color', [208, 240, 192]/255, ...
            'Linewidth', 2, 'HandleVisibility','off')
        plot(time, c_tot, 'k', 'Linewidth', 2, 'Displayname', 'True')
        plot(time, ct_rec_bcm_mean, 'g--', 'Linewidth', 2, ...
            'Displayname', 'Rec mean')
        plot(time, ct_rec_bcm_opt, 'g', 'Linewidth', 2, ...
            'Displayname', 'Rec opt')
        ylim([0, y_max])
        if ii == 1
            title('3 comp', 'Fontsize', 18)
        end
        xlabel('Time [min]', 'Fontsize', 18)
        lgd = legend('show');
        set(lgd, 'Fontsize', 18, 'Location', 'SouthEast')
        
    end
    
    saveas(f_ris, fullfile(folder_figures, ...
                    sprintf('result_%s.png', mouse.name)))
                
%% Step 5. Plot concentrations
    f_conc = figure('units', 'normalized', 'outerposition',[0 0 1 1]);
    conc = {'C_f(t) [kBq/mL]', 'C_p(t) [kBq/mL]', 'C_r(t) [kBq/mL]'};
    % skf
    for ic = 1:2
        subplot(3, 2, 2*ic-1)
        hold on
        plot(time, c_rec_skf_opt.ant(ic, :), 'k', 'Linewidth', 2, ...
            'Displayname', 'Ant');
        plot(time, c_rec_skf_opt.post(ic, :), 'r', 'Linewidth', 2, ...
            'Displayname', 'Post');
        set(gca, 'Fontsize', 13)
        lgd = legend('show');
        set(lgd, 'Location', 'Best', 'Fontsize', 15)
        ylabel(conc{ic}, 'Fontsize', 18)
        if ic == 2
            xlabel('Time t [min]', 'Fontsize', 18)
        elseif ic == 1
            title(sprintf('%s - skf', mouse.name), 'Fontsize', 18, 'Interpreter', 'None')
        end
    end
    % bcm
    for ic = 1:3
        subplot(3, 2, 2*ic)
        hold on
        plot(time, c_rec_bcm_opt.ant(ic, :), 'k', 'Linewidth', 2, ...
            'Displayname', 'Ant');
        plot(time, c_rec_bcm_opt.post(ic, :), 'r', 'Linewidth', 2, ...
            'Displayname', 'Post');
        set(gca, 'Fontsize', 13)
        lgd = legend('show');
        set(lgd, 'Location', 'Best', 'Fontsize', 15)
        ylabel(conc{ic}, 'Fontsize', 18)
        if ic == 1
            title('bcm', 'Fontsize', 18)
        elseif ic == 3
            xlabel('Time t [min]', 'Fontsize', 18)
        end
    end
    
    saveas(f_conc, fullfile(folder_figures, ...
                   sprintf('result_conc_%s.png', mouse.name)))

    pause(0.1)
    close(f_ris)
    close(f_conc)
    
end

%% Step 5. Save
ris_all.name_k_skf = name_k_skf;
ris_all.name_k_bcm = name_k_bcm;
ris_all.skf_diab = skf_diab;
ris_all.skf_control = skf_control;
ris_all.bcm_diab = bcm_diab;
ris_all.bcm_control = bcm_control;

save(fullfile(folder_results, 'ris_all_mouse'), 'ris_all')
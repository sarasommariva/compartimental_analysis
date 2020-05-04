clc
clear
close all

%% Step 1. Define general parameter
% 1.a. Path
folder_results = fullfile('.', 'results');
folder_figures = fullfile('.', 'figures');

path_results = fullfile(folder_results, 'ris_all_mouse');

%% Step 2. Load
load(path_results, 'ris_all')

name_k_skf = ris_all.name_k_skf;
name_k_bcm = ris_all.name_k_bcm;

n_an_control = size(ris_all.skf_control.ant.opt, 1);
n_an_diab = size(ris_all.skf_diab.ant.opt, 1);

type_k = 'opt';

mean_k_skf_ant_diab = mean(ris_all.skf_diab.ant.(type_k), 1);
mean_k_skf_post_diab = mean(ris_all.skf_diab.post.(type_k), 1);
mean_k_skf_ant_control = mean(ris_all.skf_control.ant.(type_k), 1);
mean_k_skf_post_control = mean(ris_all.skf_control.post.(type_k), 1);

std_k_skf_ant_diab = std(ris_all.skf_diab.ant.(type_k), [], 1);
std_k_skf_post_diab = std(ris_all.skf_diab.post.(type_k), [], 1);
std_k_skf_ant_control = std(ris_all.skf_control.ant.(type_k), [], 1);
std_k_skf_post_control = std(ris_all.skf_control.post.(type_k), [], 1);

mean_k_bcm_ant_diab = mean(ris_all.bcm_diab.ant.(type_k), 1);
mean_k_bcm_post_diab = mean(ris_all.bcm_diab.post.(type_k), 1);
mean_k_bcm_ant_control = mean(ris_all.bcm_control.ant.(type_k), 1);
mean_k_bcm_post_control = mean(ris_all.bcm_control.post.(type_k), 1);

std_k_bcm_ant_diab = std(ris_all.bcm_diab.ant.(type_k), [], 1);
std_k_bcm_post_diab = std(ris_all.bcm_diab.post.(type_k), [], 1);
std_k_bcm_ant_control = std(ris_all.bcm_control.ant.(type_k), [], 1);
std_k_bcm_post_control = std(ris_all.bcm_control.post.(type_k), [], 1);


%% Step 3. SKF
for ik = 1:numel(name_k_skf)
    skf.(name_k_skf{ik}) = ...
        [mean_k_skf_ant_diab(ik), mean_k_skf_post_diab(ik); ...
        mean_k_skf_ant_control(ik), mean_k_skf_post_control(ik)];
    std_skf.(name_k_skf{ik}) = ...
        [std_k_skf_ant_diab(ik), std_k_skf_post_diab(ik); ...
        std_k_skf_ant_control(ik), std_k_skf_post_control(ik)];
end

f_skf = figure('units', 'normalized', 'outerposition',[0 0 1 1]);
for ik = 1:numel(name_k_skf)
    subplot(2, 2, ik)
    bar(skf.(name_k_skf{ik}), 'grouped');
    hold on
    
    set(gca,'XTickLabel',{'Diab','Control'}, 'Fontsize', 18);
    title(name_k_skf{ik}, 'Fontsize', 18)
    
    ngroups = 2; nbars = 2; groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, skf.(name_k_skf{ik})(:, i), std_skf.(name_k_skf{ik})(:, i), ...
            'k', 'linestyle', 'none');
    end
    
end

%% Step 4. BCM
for ik = 1:numel(name_k_bcm)
    bcm.(name_k_bcm{ik}) = ...
        [mean_k_bcm_ant_diab(ik), mean_k_bcm_post_diab(ik); ...
        mean_k_bcm_ant_control(ik), mean_k_bcm_post_control(ik)];
    std_bcm.(name_k_bcm{ik}) = ...
        [std_k_bcm_ant_diab(ik), std_k_bcm_post_diab(ik); ...
        std_k_bcm_ant_control(ik), std_k_bcm_post_control(ik)];
end

f_bcm = figure('units', 'normalized', 'outerposition',[0 0 1 1]);
for ik = 1:numel(name_k_bcm)
    subplot(3, 2, ik)
    bar(bcm.(name_k_bcm{ik}), 'grouped');
    hold on
    
    set(gca,'XTickLabel',{'Diab','Control'}, 'Fontsize', 18);
    title(name_k_bcm{ik}, 'Fontsize', 18)
    
%     ngroups = 2; nbars = 2; groupwidth = min(0.8, nbars/(nbars + 1.5));
%     for i = 1:nbars
%         x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%         errorbar(x, bcm.(name_k_bcm{ik})(:, i), std_bcm.(name_k_bcm{ik})(:, i), ...
%             'k', 'linestyle', 'none');
%     end
    
end

%% Step 5. Plot mouse by mouse
f_sing = figure('units', 'normalized', 'outerposition',[0 0 1 1]);

for jj = 1:numel(name_k_skf)
    subplot(numel(name_k_skf), 2, 2*jj-1)
    plot(1:n_an_diab, ris_all.skf_diab.ant.opt(:, jj), ...
            'Linewidth', 2, 'Displayname', 'Ant')
    hold on
    plot(1:n_an_diab, ris_all.skf_diab.post.opt(:, jj), ...
            'Linewidth', 2, 'Displayname', 'Post')
    
    subplot(numel(name_k_skf), 2, 2*jj)
    plot(1:n_an_control, ris_all.skf_control.ant.opt(:, jj), ...
            'Linewidth', 2, 'Displayname', 'Ant')
    hold on
    plot(1:n_an_control, ris_all.skf_control.post.opt(:, jj), ...
            'Linewidth', 2, 'Displayname', 'Post')
    
end

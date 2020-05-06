%% Domanda: 
% (1) Devo controllare le unità di misura
% (2) Fare anche la parte wb?
% (3) Ho sbagliato il nome dato ai campi di ris 
% (4) Forse dovrei salvare i parametri utilizzati (Vb, Vi, v, da qlc parte)

clc
clear
close all

%% Step 1. Define general parameter
% 1.a. Path
folder_func = fullfile('..', 'functions');
folder_results = fullfile('.', 'results');

path_data = fullfile(folder_results, 'diab_vs_control_preproc.mat');

% 1.b. Parameters for the analysis
num_rip = 20;
Vb = 0.04;
Vi = 0;
v = 0.17;

t_0 = 0;
C_0_skf = [0, 0];
C_0_bcm = [0, 0, 0];

fields = {'ant', 'post'};

%% Step 2. Load
addpath(folder_func)
load(path_data, 'mice_group_diab', 'mice_group_control', 'time')
time_or = time'; clear time

idx_an_diab = [1:numel(mice_group_diab)];
idx_an_control = [1:numel(mice_group_control)];

%% Step 3. Analysis mouse by mouse
n_an_diab = numel(idx_an_diab);
n_an_control = numel(idx_an_control);

for im = 1:n_an_diab+n_an_control
    if im <= n_an_diab
        aux_im = idx_an_diab(im);
        mouse = mice_group_diab(aux_im);
    else
        aux_im = idx_an_control(im-n_an_diab);
        mouse = mice_group_control(aux_im);
    end
    
    c_if = mouse.if;
    time = time_or(1:numel(c_if));
    c_if = @(tt)(interp1([0 time],[0 c_if'], tt,'linear',0));
    
    for ii = 1:numel(fields)

    fprintf('Analysing mouse %s - Condition %s \n', mouse.name, fields{ii})

    c_tot = mouse.(fields{ii}); % <--- Devo rendere sistematica anche la condizione

    ris.(fields{ii}).k1x_skf = zeros(num_rip, 1); 
    ris.(fields{ii}).k2x_skf = zeros(num_rip, 1); 
    ris.(fields{ii}).k3x_skf = zeros(num_rip, 1);
    ris.(fields{ii}).k4x_skf = zeros(num_rip, 1);
    ris.(fields{ii}).relerr_skf = zeros(num_rip, 1);
    ris.(fields{ii}).iter_skf = zeros(num_rip, 1);
    ris.(fields{ii}).c_totx_skf = zeros(numel(time), num_rip);

    ris.(fields{ii}).k1x_bcm = zeros(num_rip, 1);
    ris.(fields{ii}).k2x_bcm = zeros(num_rip, 1);
    ris.(fields{ii}).k3x_bcm = zeros(num_rip, 1);
    ris.(fields{ii}).k5x_bcm = zeros(num_rip, 1);
    ris.(fields{ii}).k6x_bcm = zeros(num_rip, 1);
    ris.(fields{ii}).relerr_bcm = zeros(num_rip, 1);
    ris.(fields{ii}).iter_bcm = zeros(num_rip, 1);
    ris.(fields{ii}).c_totx_bcm = zeros(numel(time), num_rip);

    for it = 1:num_rip
        fprintf('Run = %d \n', it)
        [ris.(fields{ii}).k1x_skf(it), ris.(fields{ii}).k2x_skf(it), ...
            ris.(fields{ii}).k3x_skf(it), ris.(fields{ii}).k4x_skf(it), ...
            ~, ris.(fields{ii}).c_totx_skf(:, it), ris.(fields{ii}).relerr_skf(it), ...
            ris.(fields{ii}).iter_skf(it)] = ...
            reconstruction_Skf(c_tot, c_if, Vb, Vi, time, t_0, C_0_skf);
        
        [ris.(fields{ii}).k1x_bcm(it), ris.(fields{ii}).k2x_bcm(it), ...
            ris.(fields{ii}).k3x_bcm(it), ris.(fields{ii}).k5x_bcm(it), ...
            ris.(fields{ii}).k6x_bcm(it), ~, ris.(fields{ii}).c_totx_bcm(:, it), ...
            ris.(fields{ii}).relerr_bcm(it), ris.(fields{ii}).iter_bcm(it)] = ...
            reconstruction_BCM(c_tot, c_if, Vb, Vi, v, time, t_0, C_0_bcm); 

    end

    clear c_tot

    end

    save(fullfile(folder_results, sprintf('result_%s_temp', mouse.name)), 'ris')

    clear c_if

end






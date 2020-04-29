%% Domanda: 
% (1) Devo fare qualche correzione dei dati:
%   - dati negativi?
%   - sorta di eliminazione dei picchi alla Mara? 
% (2) Devo controllare le unità di misura
% (3) Quello che fa Mara è lanciare 50 volte l'algoritmo e poi prendere
% come risultato la media dei risultati.

clc
clear
close all

%% Step 1. Define general parameter
% 1.a. Path
folder_func = fullfile('..', 'functions');
folder_results = fullfile('.', 'results');

path_data = fullfile(folder_results, 'diab_vs_control_preproc.mat');

% 1.b. Parameters for the analysis
num_rip = 50;
% num_rip = 1;
Vb = 0.04;
Vi = 0;
v = 0.17;

t_0 = 0;
C_0_skf = [0, 0];
C_0_bcm = [0, 0, 0];

% idx_an_diab = [5, 6, 7];
% idx_an_control = [3, 4, 10];
idx_an_diab = [5];

%% Step 2. Load
addpath(folder_func)
load(path_data, 'mice_group_diab', 'mice_group_control', 'time')
time = time'; 

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
    
    fprintf('Analysing mouse %s \n', mouse.name)
    
    c_if = mouse.if(1:numel(time));
    c_tot = mouse.post(1:numel(time)); % <--- Devo rendere sistematica anche la condizione
    
    c_if = @(tt)(interp1([0 time],[0 c_if'], tt,'linear',0));
    
    % Initializa
    ris.c_tot = c_tot;
    
    ris.k1x_skf = zeros(num_rip, 1); 
    ris.k2x_skf = zeros(num_rip, 1); 
    ris.k3x_skf = zeros(num_rip, 1);
    ris.k4x_skf = zeros(num_rip, 1);
    ris.relerr_skf = zeros(num_rip, 1);
    ris.iter_skf = zeros(num_rip, 1);
    ris.c_totx_skf = zeros(numel(time), num_rip);
    
    ris.k1x_bcm = zeros(num_rip, 1);
    ris.k2x_bcm = zeros(num_rip, 1);
    ris.k3x_bcm = zeros(num_rip, 1);
    ris.k5x_bcm = zeros(num_rip, 1);
    ris.k6x_bcm = zeros(num_rip, 1);
    ris.relerr_bcm = zeros(num_rip, 1);
    ris.iter_bcm = zeros(num_rip, 1);
    ris.c_totx_bcm = zeros(numel(time), num_rip);
    
    for it = 1:num_rip
        fprintf('Run = %d \n', it)
        [ris.k1x_skf(it), ris.k2x_skf(it), ris.k3x_skf(it), ris.k4x_skf(it), ...
            ~, ris.c_totx_skf(:, it), ris.relerr_skf(it), ris.iter_skf(it)] = ...
            reconstruction_Skf(c_tot, c_if, Vb, Vi, time, t_0, C_0_skf);
        
        [ris.k1x_bcm(it), ris.k2x_bcm(it), ris.k3x_bcm(it), ris.k5x_bcm(it), ris.k6x_bcm(it), ...
            ~, ris.c_totx_bcm(:, it), ris.relerr_bcm(it), ris.iter_bcm(it)] = ...
            reconstruction_BCM(c_tot, c_if, Vb, Vi, v, time, t_0, C_0_bcm); 

    end

    save(fullfile(folder_results, sprintf('result_%s', mouse.name)), 'ris')
    
    clear c_if c_tot
    
end




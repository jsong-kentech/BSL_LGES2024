%% 
% This function 
% (1) finds the best-fit parameters for a halfcells
% calls a EIS model function

% Version
% V1_3E_combined: EIS_function_V1_3E_combined: cathode-anode combined fitting
% more outputs: func, data, plots


clear; clc; close all


%% Configurations

% EIS data path
    path_folder = 'G:\Shared drives\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';
    %path_folder = 'G:\Shared drives\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 2';
    name_file_n = 'PEIS_C09_anode_cycle_soc90.csv';
    name_file_p = 'PEIS_C09_cathode_cycle_soc90.csv';

% SOC and T (for initial guess - they are functions of soc, T)
    soc = 0.9; % [1]
    T = 298.15; %[K]

% Fitting configuration
    type_weight = 1; % 0 for absolute error, 1 for relative error
    type_acf =3; % 1 for anode, 2 for cathode, 3 for full cell

% Optimization options
    options= optimset('display','iter','MaxIter',100,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');

% Parameters 
    %   1         2     3     4      5     6         7      8     9   10    11     12 
    % R_itsc_n, i0_n, Cdl_n, Ds_n, Av_n, R_itsc_p, i0_p, Cdl_p, Ds_p, Av_p, k_el, D_el
    
    bounds = [...
         0.5 2 % (1) R_n
         0.1 500; % (2) i0_n
         0.1 500; % (3) C_dl_n
         0.1 50; % (4) Ds_n
         0.02 10; % (5) Av_n
         0.5 2;  % (6) R_p
         0.1 50; % (7) i0_p
         0.1 50; % (8) C_dl_p
         0.1 50; % (9) Ds_p
         0.1 10; % (10) Av_p
         0.1 10; % (11) kappa_el
         0.01 10; % (12) D_el
         ]; 
    lb = bounds(:,1);
    ub = bounds(:,2);

    factors_ini = ones(size(bounds,1),1);
    

%% Structure Explained

    % Data and model func
    % for n - frequency points
    % w_data (n x 1) vector
    % [ Z_n_re Z_n_im Z_p_re Z_p_im] (n x 4) matrix

%% Load and Pre-processing Data

    % load EIS data
    data_n = load([path_folder filesep name_file_n]);
    data_p = load([path_folder filesep name_file_p]);
    f_data = [data_n(:,1) data_p(:,1)];
        if any(f_data(:,1)~=f_data(:,1))
            error('freq range of EIS_p and EIS_n are different')
        end
    z_data = [data_n(:,2) data_n(:,3) data_p(:,2) data_p(:,3)];
           % [ Z_n_re     Z_n_im      Z_p_re      Z_p_im]    (n x 4) matrix

    figure(1)
    plot(data_n(:,2),-data_n(:,3),'ro'); hold on; grid on
    plot(data_p(:,2),-data_p(:,3),'bo'); 
    
    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    % trim the high-frequency inductance part (positive imag data)
    ind_keep = data_n(:,3)<=0 | data_p(:,3)<=0;
    f_data = f_data(ind_keep);
    z_data = z_data(ind_keep,:);
    figure(1)
    plot(z_data(:,1),-z_data(:,2),'ro','MarkerFaceColor','red'); hold on; grid on
    plot(z_data(:,3),-z_data(:,4),'bo','MarkerFaceColor','blue'); 

    

%% Define weighting vector
    
    if type_weight == 1 % if minimizing the relative error
    weight_n = (z_data(:,1).^2 + z_data(:,2).^2).^(-0.5);
    weight_p = (z_data(:,3).^2 + z_data(:,4).^2).^(-0.5);
    weight_matrix = [weight_n weight_n weight_p weight_p];    
    elseif type_weight == 0 % if minimizing the absolute error
    weight_matrix = ones(size(z_data));
    end




%% FITTING
%   Call EIS model
   weighted_model = @(factors,f_data)BSL_func_EISmodel_V1_3E(f_data, factors,soc,T,type_acf)...
       .*weight_matrix;
   weighted_data = z_data.*weight_matrix;
%   fitting
   tic;
   [factors_hat, resnorm,residual,~,~,~,jacobian_hat] ...
        = lsqcurvefit(weighted_model,factors_ini,...
                      f_data,weighted_data, lb, ub, options);
   toc;


%% Plot Results
[z_model0, paras0] = BSL_func_EISmodel_V1_3E(f_data,factors_ini,soc,T,type_acf);
[z_model1, paras1] = BSL_func_EISmodel_V1_3E(f_data,factors_hat,soc,T,type_acf);

% Nyquist Plot
    % line colors
    cmat_jet = jet(16);
    cmat_n = [cmat_jet(16,:);cmat_jet(14,:);cmat_jet(11,:)];
    cmat_p = [cmat_jet(1,:);cmat_jet(3,:);cmat_jet(6,:)];

figure(2)
t = tiledlayout(2,2);
nexttile
plot(z_data(:,1),-z_data(:,2),'o','linewidth',1,'color',cmat_n(1,:)); hold on
plot(z_model0(:,1),-z_model0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
plot(z_model1(:,1),-z_model1(:,2),'o','linewidth',1,'color',cmat_n(2,:))
legend('Exp Data','Model Initial','Model Fit')
%legend('Exp Data','Model Fit')
daspect ([1 1 2])


    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    hold off
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

nexttile
plot(z_data(:,3),-z_data(:,4),'o','linewidth',1,'color',cmat_p(1,:)); hold on
plot(z_model0(:,3),-z_model0(:,4),'o','linewidth',1,'color',cmat_p(3,:))
plot(z_model1(:,3),-z_model1(:,4),'o','linewidth',1,'color',cmat_p(2,:))
legend('Exp Data','Model Initial','Model Fit')
%legend('Exp Data','Model Fit')
daspect ([1 1 2])

    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    hold off
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

% Zoom-in semicircle
nexttile
plot(z_data(:,1),-z_data(:,2),'o','linewidth',1,'color',cmat_n(1,:)); hold on
plot(z_model0(:,1),-z_model0(:,2),'o','linewidth',1,'color',cmat_n(3,:))
plot(z_model1(:,1),-z_model1(:,2),'o','linewidth',1,'color',cmat_n(2,:))
legend('Exp Data','Model Initial','Model Fit')
%legend('Exp Data','Model Fit')



    f_zoom_lb = 10; %[Hz] 
    idx_zoom = f_data>f_zoom_lb;
    axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    hold off
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')


nexttile
plot(z_data(:,3),-z_data(:,4),'o','linewidth',1,'color',cmat_p(1,:)); hold on
plot(z_model0(:,3),-z_model0(:,4),'o','linewidth',1,'color',cmat_p(3,:))
plot(z_model1(:,3),-z_model1(:,4),'o','linewidth',1,'color',cmat_p(2,:))
legend('Exp Data','Model Initial','Model Fit')
%legend('Exp Data','Model Fit')

    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    hold off
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'position',[100 100 1200 1000])
%% Result Summary
   
Result.factors_hat = factors_hat;
Result.paras_hat = paras1;
Result.z_model = z_model1;


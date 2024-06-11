%% 
% This function 
% (1) finds the best-fit parameters for a halfcells
% calls a EIS model function

clear; clc; close all


%% Configurations

% EIS data path
    path_folder = 'G:\Shared drives\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 1';
    %path_folder = 'G:\Shared drives\BSL-Data\LGES\LG raw data\12_6cm2_soc10_EIS # Sample 2';
    name_file = 'PEIS_C09_cathode_cycle_soc30.csv';

% SOC and T (for initial guess - they are functions of soc, T)
    soc = 0.3; % [1]
    T = 298.15; %[K]

% Simul configuration % LGES 05
    type_dist = 0; % 0 for DRT, 1 for DDT
    type_acf =2; % 1 for anode, 2 for cathode, 3 for full cell

% Optimization options
    options= optimset('display','iter','MaxIter',100,'MaxFunEvals',1e5,...
        'TolFun',1e-8,'TolX',1e-8,'FinDiffType','central');


% Parameters 
    factors_ini = [0.865241179222448	7.39543618232181	5.76847499611458	17.2448054354830	2.09466741915660	0.942836139030348	0.534603743428756];

%% Load and Pre-processing Data

    % load EIS data
    data = load([path_folder filesep name_file]);
    f_data = data(:,1);
    z_re_data = data(:,2);
    z_im_data = data(:,3);
    z_data = [z_re_data z_im_data];

    figure(1)
    plot(z_re_data,-z_im_data,'o'); hold on; grid on
    
    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')

    % trim the high-frequency inductance part
    f_data = f_data(z_im_data<=0);
    z_re_data = z_re_data(z_im_data<=0);
    z_im_data = z_im_data(z_im_data<=0);
    z_data = [z_re_data z_im_data];
    figure(1)
    plot(z_re_data,-z_im_data,'o'); hold on; grid on



%% Plot Results


% Nyquist Plot
figure(2)
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on
% Zoom-in semicircle
figure(3)
plot(z_data(:,1),-z_data(:,2),'ok','linewidth',1); hold on

% plot model with different std

std_vec = [0 0.1 0.3 0.5 1 2];
c_mat = jet(length(std_vec));


for i = 1:length(std_vec)

[z_model0, paras0] = BSL_func_EISmodel_V2_half_Dist(f_data,[factors_ini std_vec(i)],soc,T,type_acf,type_dist);

figure(2)
plot(z_model0(:,1),-z_model0(:,2),'-b','linewidth',1,'Color',c_mat(i,:))
%plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
%legend('Exp Data','Model Initial','Model Fit')
daspect ([1 1 2])

    axis_limit = 1.1*max(max(abs(z_data)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')


% Zoom-in semicircle
figure(3)
plot(z_model0(:,1),-z_model0(:,2),'-b','linewidth',1,'Color',c_mat(i,:))
%plot(z_model1(:,1),-z_model1(:,2),'or','linewidth',1)
%legend('Exp Data','Model Initial','Model Fit')


    f_zoom_lb = 10; %[Hz] 
    idx_zoom = f_data>f_zoom_lb;
    axis_limit = 1.1*max(max(abs(z_data(idx_zoom,:))));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')


end
%% Result Summary
   
%Result.factors_hat = factors_hat;
%Result.paras_hat = paras1;
%Result.z_model = z_model1;


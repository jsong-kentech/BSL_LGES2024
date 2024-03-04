%1. dUdc
clear; clc; close all

addpath('C:\Users\jsong\Documents\MATLAB\BSL_EIS\1_standalone\interpolations')


soc = linspace(-1,2,101);


x_1 = 0.8781; % anode stoic when soc =0
x_0 = 0.0216; % anode stoic when soc =1
y_0 = 0.9319; % cathode stoic when soc =0
y_1 = 0.3532; % cathode stoic when soc =1

x = x_0 + (x_1 - x_0)*soc; % anode stoic
y = y_0 + (y_1 - y_0)*soc; % cathode stoic


Ua = Ua_function_v2(x,0.5);
Uc = Uc_function_v2(y,0.5);

cmat = lines(3);

figure(1)
subplot(1,4,1)
plot(x,Ua,'color',cmat(2,:))
xlim([0 1])
ylim([0 1])

subplot(1,4,2)
plot(y,Uc,'color',cmat(1,:))
xlim([0 1])
ylim([3 4.3])

subplot(1,4,3)
yyaxis left
plot(soc,Uc,'color',cmat(1,:))
yyaxis right
plot(soc,Ua,'color',cmat(2,:))
xline(0)
xline(1)



subplot(1,4,4)
yyaxis left
plot(soc,-gradient(Uc,y),'color',cmat(1,:))

yyaxis right
plot(soc,-gradient(Ua,x),'color',cmat(2,:))
xline(0)
xline(1)



%%
ctc = 48786;
cs0c = y*ctc;
c_e_ref = 1000; 
F = 96487;                  % OK
c_e = 1120; 
alphaa = 0.5;
alphac = 0.5;
kc = kc_function(y,298.15,ctc); % {modified}

i0c = F*kc.*((c_e/c_e_ref).^alphaa).*((ctc-cs0c).^alphaa).*cs0c.^alphac;                    % {modified} c_e_ref


Dsc = Dsc_function(y,298.15);

figure(2)
subplot(1,2,1)
semilogy(y(~isnan(Uc)),i0c(~isnan(Uc)))
xlim([0 1])
ylim([10^-2 10])
subplot(1,2,2)
semilogy(y(~isnan(Uc)),Dsc(~isnan(Uc)))
xlim([0 1])



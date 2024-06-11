R = 1;
C = 1;

log10_w = linspace(-3,3,61);
w = 10.^log10_omega;


Y = R^-1 + 1i*w*C;
Z = Y.^-1;


figure(1)
plot(real(Z),-imag(Z)); hold on
daspect ([1 1 2])

Z_drt = zeros(size(w));
Y_drt = zeros(size(w));

for n = 1:length(w)

    Z_drt(n) = z_total_drt(w(n),R,C,1);
    Y_drt(n) = y_total_drt(w(n),R,C,1);

end

figure(1)
plot(real(Z_drt),-imag(Z_drt),'ob'); hold on
plot(real(Y_drt.^-1),-imag(Y_drt.^-1),'or'); hold on




function z_drt = z_total_drt(w,R,C,std)

[mu,sigma] = normal_para(1,std);

upper_func = @(t)lognpdf(t,mu,sigma).*z_local_drt(t,w,R,C);
lower_func = @(t)lognpdf(t,mu,sigma);

upper = integral(upper_func,10^-3,10^3);
lower = integral(lower_func,10^-3,10^3);

lower

z_drt = upper./lower;


end


function z_loc = z_local_drt(t,w,R,C)


Y = R^-1 + 1i*w*C*t;
z_loc = Y.^-1;

end


function y_drt = y_total_drt(w,R,C,std)

[mu,sigma] = normal_para(1,std);

upper_func = @(t)lognpdf(t,mu,sigma).*y_local_drt(t,w,R,C);
lower_func = @(t)lognpdf(t,mu,sigma);

upper = integral(upper_func,10^-3,10^3);
lower = integral(lower_func,10^-3,10^3);

lower

y_drt = upper./lower;


end





function y_loc = y_local_drt(t,w,R,C)


Y = R^-1 + 1i*w*C*t;
y_loc = Y;

end




function [mu,sig]=normal_para(mean,std) % LGES V2 2024 05
% converting parameters of observed lognormal distribution
% to parameters of underlying normal distribution
m=mean;
v=std^2;
mu=log(m^2/(v+m^2)^0.5);
sig=(log(v/m^2+1))^0.5;


end

% compute the errors of NQ
f = @(x) exp(2*x) + 2*sin(x) + x.^4 + log(abs(x)) + exp(-x.*x);
% x = linspace(0.01,2,100); figure();plot(x,f(x));

q = 3; % # quad pts
[xi, wi] = GH_quad(q);

f_df  = @(x) 2*exp(2*x) + 2*cos(x) + 4*x.^3 + 1./x - 2.*x.*exp(-x.*x);
f_d2f = @(x) 4*exp(2*x) - 2*sin(x) + 12*x.^2 - 1./(x.^2) - 2*exp(-x.*x) + 4*x.^2.*exp(-x.*x) ;
f_d3f = @(x) 8*exp(2*x) - 2*cos(x) + 24*x + 2/x^3 + 12*x*exp(-x*x) - 8*x^3*exp(-x*x);
f_d4f = @(x) 16*exp(2*x) + 2*sin(x) + 24 - 6/x^4 + 12*exp(-x*x) - 48*x*x*exp(-x*x) + 16*x^4*exp(-x*x);
f_d5f = @(x) 32*exp(2*x) + 2*cos(x)  + 24/x^5 - 120*x*exp(-x*x) + 160*x^3*exp(-x*x) - 32*x^5*exp(-x*x);
f_d6f = @(x) 64*exp(2*x) - 2*sin(x)  - 120/x^6 - 120*exp(-x*x) + 720*x^2*exp(-x*x) - 480*x^4*exp(-x*x) + 64*x^6*exp(-x*x);
x = 1.81; % pt of interest
fx = f(x); % compute in advance
v = 1; % 1 if 1-D

hvec = 1.25.^(-15:-4); % this varies
len = length(hvec);
approx_vals = zeros(6,len);
correct_vals = [
    f_df(x);
    f_d2f(x);
    f_d3f(x);
    f_d4f(x);
    f_df(x);
    f_d2f(x)
    ] *ones(1,len); % 6-by-len
for i = 1:len
    [approx_vals(1,i), approx_vals(2,i), approx_vals(3,i), approx_vals(4,i), ~, ~] = GH_Deriv(q, xi, wi, f, hvec(i), v, x, fx);
    % approximated through NQ
    approx_vals(5,i) = (f(x+hvec(i))-f(x-hvec(i)))/hvec(i)/2;
    approx_vals(6,i) = (f(x+hvec(i))-2*fx+f(x-hvec(i)))/hvec(i)^2;
    % approximated through FD
    % what if...
    correct_vals(1,i) = smoothing(f_df,x,hvec(i));
%     correct_vals(5,i) = correct_vals(1,i);
    correct_vals(2,i) = smoothing(f_d2f,x,hvec(i));
%     correct_vals(6,i) = correct_vals(2,i);
end
err = approx_vals - correct_vals; % compute the diff

figure();hold on;
disp = 2;
for j=1:disp
    plot(log(hvec),log(abs(err(j,:))), '-*');
end
legendnames = {'df','d2f','d3f','d4f','FD1', 'FD2'};
legend(legendnames(1:disp), 'Location', 'best');

function gfx = smoothing(f, x, mu)
rho = @(x) 1/sqrt(2*pi) * exp(-x.^2/2);
fun = @(t) f(x+mu*t).*rho(t);
gfx = integral(fun, -Inf, +Inf);
end
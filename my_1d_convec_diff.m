%% 1D Burgers' Equation 'du/dt + u(du/dx) = v(d2u/dx2)'
% convection and diffusion terms combined
% Discretization:
% Forward difference for du/dt
% Backward difference for first du/dx
% Central difference for d2u/dx2
%
% Discretize:
% (u(i,n+1) - u(i,n))/(delta t) + u(i,n)(u(i,n)-u(i+1,n))/(delta x) =...
%    v*(u(n,i+1)-2u(i,n)+u(i-1,n))/(delta x)^2
% Transpose (solve for only unkonwn)
% u(i,n+1) = u(i,n) ...
%    - u(i,n)*(delta t)/(delta x)*(u(i,n)-u(i-1,n)) ...
%    + v(delta t)/(delta x)^2 * (u(n,i+1)-2u(i,n)+u(i-1,n))
%

% Given initial condition:
% u = -2v*((d phi)/dx )/phi + 4
% phi = exp(-x^2/(4v)) + exp((-x-2pi)^2/4v)
clear
clc

v = 0.1;

syms('x')
phi = exp(-(x^2)/(4*v)) + exp(-(x-2*pi)^2/(4*v));
dp = diff(phi, x);

% and Boundary condition is u(0) = u(2*pi)
%    (function is periodic)
nx = 20;
x = linspace(0, 2*pi, nx);
i = x;
u = -2*v.*(subs(dp))./subs(phi)+4;
u = double(u);
p = plot(i, u);
axis([0,i(end), 0, max(u)])
dx = i(2)-i(1);

nt = 100;
t_max = 1.5;
n = linspace(0, t_max, nt);
dt = n(2)-n(1);

un = u;
p = plot(i, un);
axis([0,i(end), 0, max(u)])


for t = 1:(nt*t_max)
    % u(i,n+1) = u(i,n) ...
    %    - u(i,n)*(delta t)/(delta x)*(u(i,n)-u(i-1,n)) ...
    %    + v(delta t)/(delta x)^2 * (u(n,i+1)-2u(i,n)+u(i-1,n))
    un = [un(1),...
        un(2:end-1)...
        - un(2:end-1) .* dt/dx .* (un(2:end-1)-un(1:end-2))...
        + v * dt/dx^2 .*(un(3:end)-2*un(2:end-1)+un(1:end-2))...
        , un(end)];
    set(p,'YData',un);
    title(['time = ',num2str(t/nt)])
    drawnow()
end





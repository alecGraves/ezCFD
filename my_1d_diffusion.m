%% 1d diffusion (du/dt) = v(d2u/dx2)
% (heat eq. if u is the temp!)
%
% v is a constant (viscosity) => exact solutions are known
%
% Consider solutions of type u = Ue^(i*(kx-wt)), i=sqrt(-1)
%    Which represents a plane wave of amplitude U, wave number k, 
%    and w=2*pi*f and f is frequency
%
% Gives us iw = vk^2
% -> solution: u = Ue^(ikx)*e^(-tvk^2)
% --> a wave exponentially dampened in time by e^(-tvk^2)
%    where v is the diffusion coefficient. 
%    v has to be positive for physical dampening 
%    (neg. v would represent exponential growth,
%    like an explosion or model for concentration of wealth in wealthy)
%
% Pure diffusion represents an exponentially dampened wave
%
% The physics of diffusion is isotropic, so best finite-difference
%    formula for derivative approximation is the central difference.
%
% numerical methods:
% forward difference in time
% central difference in space
%
% Discretization:
% (u(i,n+1)-u(i,n))/(delta t)= ...
%     v*(u(n,i+1)-2u(i,n)+u(i-1,n))/(delta x)^2
% Transoposed form (solved for u(i, n+1))
% u(i,n+1) = u(i,n) + v(delta t)/(delta x)^2 ...
%            * u(n,i+1)-2u(i,n)+u(i-1,n)
%
clear();
clc();

rx = 44; % resolution
rt = 200;
t_max = 1;
i = linspace(0, 2, rx); % space
n = linspace(0, t_max, rt); % time

u = ones(1, rx);
u((0.5 <= i) & (i <= 1)) = 2;
% p = plot(i, u);

v = .1; % viscosity

% assign time intervals:
dt = n(2)-n(1);
dx = i(2)-i(1);
% v*dt/dx^2

un = u;
p = plot(i, un);
axis([0 2 0 2]);

for t = 1:(rt*t_max)
    % u(i,n+1) = u(i,n) + v(delta t)/(delta x)^2 ...
    %            * (u(n,i+1)-2u(i,n)+u(i-1,n))
    un = [un(1), un(2:end-1) + v*dt/dx^2 .* ...
        (un(3:end)-2*un(2:end-1)+un(1:end-2)), un(end)];
    set(p,'YData',un);
    title(['time = ',num2str(t/rt)])
    drawnow()
end




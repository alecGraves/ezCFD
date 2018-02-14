%% 1d linear convection example:
% du/dt + c(du/dx) = 0
% 
% c = constant transport velocity
% u(x, t=0) = u0(x)
% After time t, u(x) is just the initial profile
%    displaced with distance x = ct
%
% Solution:
% u(x,t) = u0(x-ct)
% -> a Wave propogation equation
% -> 'u' is a wave amplitude with a face propogation speed 'c'
% 
% This equation underlies a relationship between pure convetion 
%    and wave propogation.
%
clear
clc

% space-time discretization


% index of grid in x
r = 100; % x resolution
i = linspace(0, 2, r);
dx = i(2)-i(1);

% time index
t_max = 1;
n = linspace(0, t_max, r);
dt = n(2)-n(1);

% For this example, we will use 
%    forward difference in time (for du/dt)
%    and backward difference in space (for du/dx)
%
% (u(i,n+1) - u(i,n))/(delta t) + c*(u(i,n) - u(i-1,n))/(delta x) = 0
% given initial condition u(i,n), solve for next time u(i,n+1):
% u(i,n+1) = u(i,n) - c*(delta t)/(delta x)*(u(i,n)-u(i-1,n))
% 
% initial condition: u = 2 at 0.5 <= x <= 1
%                    u = 1 elsewhere
%                    u=1 at x=0, x=2 (boundary condition for sim)

c = 1; % wavespeed of 1

% Initial conditions
u = ones(1, r);
u((0.5 <= i) & (i <= 1)) = 2;


un = u;
p = plot(i, un);
axis([0 2 0 2])

for t = 1:(r*t_max)
    un = [un(1), un(2:end) - (c*dt/dx).*(un(2:end)-un(1:end-1))];
    set(p,'YData',un);
    title(['time = ',num2str(t/r)])
    drawnow()
end

% the end, my first fluid simulation

%% Inviscid Burgers: (du/dt) + u(du/dx) = 0
%    -> Can gernerate discontinuous solutions from
%    smooth initial conditions (similar to shocks in supersonic 
%    flows).


rate = 1/r;
un = u;
p = plot(i, un);
axis([0 2 0 2])

for t = 1:(r*t_max)
    tic()
    un = [un(1), un(2:end) - (un(2:end)*dt/dx).*(un(2:end)-un(1:end-1))];
    set(p,'YData',un);
    title(['time = ',num2str(t/r)])
    drawnow()
end

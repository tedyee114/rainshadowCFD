%% System parameters
dx = 1/40;
dy = 1/40;
x_min = 0; x_max = 5;
y_min = 0; y_max = 1;

%% Set obstruction geometry data
x_obs_front = .5*(x_max-x_min);
y_obs_top = 0.5*(y_max-y_min);

% Temporal Domain
dt = 0.2;
t = 0:dt:2;
Nt = length(t);

% Set the Reynolds number
Re = 1;

% Set the relaxation parameters
alpha_u = .2;
alpha_v = .2;
alpha_p = .0025;

% Set the boundary speeds
speeds.u_bot =  0; speeds.v_bot = 0;
speeds.u_top =  0; speeds.v_top = 0;
speeds.u_lef =  1; speeds.v_lef = 0; 
speeds.u_rig =  0; speeds.v_rig = 0; 

% Set the stop conditions
II_max = 300;
u_tol = 10^-6;
II = 1;
u_change = 1;
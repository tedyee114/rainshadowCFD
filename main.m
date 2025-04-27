clear; close all; clc;
tic

%% 
cell_visualizations_on_off =1;
blocked = 1;
rainshadow = 0;
noslipground = 0;

%%%%ideally a geometry gui would set the values currently in
%%%%system_parameters, then once submitted would generate, then gui would
%%%%setup methods then begin once the initialize and run button is clicked
system_parameters                                     %% Set system, geometry, and solver parameters
solver_gui
phi_in = zeros(N.y_p);
phi_wall = 1;
if blocked ==1
    set_cell_type_blocked
else
    set_cell_type
end
if rainshadow ==1
    %% Mountain Parameters
    rise_slope = .5;
    fall_slope = -1;
    x_mountain_start = floor(N.x_p*50/100);
    x_mountain_peak = floor(N.x_p*60/100);
    x_mountain_end = floor(N.x_p*120/100);
    ambienthumidity = 0;
    scalefactor     = 10000;                    %Default = 10000. Vertical scale factor in simulation length units to m, i.e., 1 simulation length unit is 10km

    %%Atmospheric Parameters
    T_surface       = 300;                      %Default = 300 Kelvin. Temperature at sea level
    lapse_rate      = 9.8/1000*scalefactor;     %Default = 9.8/1000 Kelvin/meter.  Lapse rate
    criticalmass    = .025;                      %Default = 3.2 for phi_wall=1,phi_in=0. Nonrealistic parameter controlling rain threshold
    rainremovalrate = .7;                       %default = 0.6 under 0.5 is 100% removal,  up to 1 is 50%. it skips over cells via round(rainremovalrate*rand)

    %% Set the heat transfer properties
    kappa = 26.3e-3/scalefactor; % 26.3e-3 for air at 300K
    phi_wall = 0; %amount of moisture coming from above atmosphere. not needed for rainshadow

    % Initialize rain gauge to track moisture cleared in each column
    rain_gauge = zeros(1, N.x_p);
    create_mountain
end
plot_cell_type
guess_initialization                                  %% Generate the initial guess for the system


for n = 1:Nt-1                                       %%%%%%%%%%%%%% TIME MARCHING
    disp(['TIME STEP NUMBER: ', num2str(n)])
    disp('time_step_initialization')
    time_step_initialization                          %% Set the variables to begin a new time step
    disp('simple_algorithm')
    simple_algorithm                                  %% Begin recursive calculation of u^n+1 and v^n+1 using A_u, A_v, and A_p
    disp('heat_transfer_step')
    heat_transfer_step                                %% Perform a heat transfer step
    disp('plot_store_results')
    plot_store_results                                %% Plot results and store the u^n+1, v^n+1, and p^n+1 fields
    if rainshadow ==1
        disp('clouds')
        clouds                                            %% Clouds and Rain
    end
end

figure(19)
hold on
x_mountain = [x_mountain_start, x_mountain_peak, x_mountain_end];
% Number of points on each side of mountain
n_left = x_mountain_peak - x_mountain_start + 1;
n_right = x_mountain_end - x_mountain_peak;

left_x = linspace(x_mountain_start, x_mountain_peak, n_left);
left_y = round(rise_slope * (left_x - x_mountain_start));
right_x = linspace(x_mountain_peak, x_mountain_end, n_right);
right_y = round(fall_slope * (right_x - x_mountain_peak));
x_values = [left_x, right_x];
y_mountain = [left_y, right_y];

x_values = [x_values, x_mountain_end];  % Add the end point for x
y_mountain = [y_mountain, 0];  % Add the bottom point (y = 0) for closing the shape
fill(x_values, y_mountain, 'k', 'FaceAlpha', 0.2);
%Quiver
xlocation = floor(2.75/dx);
x_quiver = ones(N.y_u, 1) * xlocation;
left_boundary_y = (1:N.y_u)';
quiver(x_quiver, left_boundary_y, u_guess(1:end, xlocation), v_guess(1:end-1, xlocation), 0.5, 'r', 'LineWidth', 1.5, 'AutoScale', 'off', 'AutoScaleFactor', 1000, 'MaxHeadSize', 2);
title(['Horizontal Velocity Profile at x = ', num2str(xlocation*dx), ' with No-Slip Mountain']);
xlabel('Velocity Magnitude (representative values)');
ylabel('Altitude (y-location)');
axis([0 N.x_u 0 N.y_u]);
hold off;

%%
toc
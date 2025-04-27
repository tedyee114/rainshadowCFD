%% Set cell size, assuming uniform mesh
dx = grids.x_u(2)-grids.x_u(1);
dy = grids.y_u(2)-grids.y_u(1);

%% Create the cell type matrices
u_type = zeros(N.y_u,N.x_u);
v_type = zeros(N.y_v,N.x_v);
p_type = zeros(N.y_p,N.x_p);

%% Set the u-cell types

% Cells with a south face on a no slip boundary
f_x = find( grids.x_u <= x_obs_front & grids.x_u > 0); % Points after the obstruction
u_type(1,f_x) = 4; % Points ON the top of the obstruction
f_x = find(grids.x_u >= x_obs_front & grids.x_u < grids.x_u(end));
f_y = find(grids.y_u > y_obs_top, 1, 'first');
u_type(f_y,f_x) = 4; % Points in the middle that NEIGHBOR the bottom boundary

u_type(:,1)         = 6;            %%Left/West Edge

% Cells with the east face above obstruction
f = find( grids.y_u > y_obs_top);
u_type(f,end) = 2; 

u_type(end,2:end-1) = 8;            %%Top/North Edge

% Cells with a cell center in or on the obstruction
f_y = find( grids.y_u < y_obs_top);
f_x = find( grids.x_u >= x_obs_front);
u_type(f_y,f_x) = -1; % Points in the middle that NEIGHBOR the top boundary



%% Set the v-cell type

% Cells ON the no slip boundary
f_x = find(grids.x_v < x_obs_front); % Points before the obstruction
v_type(1,f_x) = 4; 
f_x = find(grids.x_v > x_obs_front);
f_y = find( abs(grids.y_v-y_obs_top) < 10^-4 );
v_type(f_y,f_x) = 4; % Points ON the top of the obstruction

v_type(2:end-1,1)   = 6;            %%Left/West Edge
v_type(end,:)       = 8;            %%Top/North Edge

% Cells with the East face ON the outlet
f_y = find(grids.y_v > y_obs_top+dy/2 & grids.y_v < grids.y_v(end)); % Above the obstruction and below the top
v_type(f_y,end) = 2; 

% Cells with the East face on the no-slip obstruction
f_x = find(grids.x_v < x_obs_front , 1, 'last');
f_y = find( abs(grids.y_v - y_obs_top) < 10^-4);
v_type(2:f_y,f_x) = 2; 

% Cells with a cell center in the obstruction
f_y = find( grids.y_v <= y_obs_top);
f_x = find( grids.x_v >= x_obs_front);
v_type(f_y,f_x) = -1; % Points in the middle that NEIGHBOR the top boundary

%% Set the p-cell type
p_type(1,1) = 5;                   %%Southwest Corner
p_type(end,1) = 7;                 %%Northwest Corner
p_type(2:end-1,1) = 6;             %%Left/West Edge 

% Cells with boundaries on the South and East face
f_x = find(grids.x_p < x_obs_front,1,'Last');
p_type(1,f_x) = 3;
f_y = find(grids.y_p > y_obs_top,1,'First');
p_type(f_y,end) = 3;

p_type(end,end) = 1;               %%Northeast Corner

% Cells with boundaries on the East face
p_type(f_y+1:end-1,end) = 2;
p_type(2:f_y-1,f_x) = 2;

% Cell with boundaries on the South face
p_type(1,2:f_x-1) = 4;
p_type(f_y,f_x+1:end-1) = 4;

p_type(end,2:end-1) = 8;            %%Top/North Edge

p_type(1:f_y-1,f_x+1:end) = -1;     %%Obstruction Cells
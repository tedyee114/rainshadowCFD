%% Set cell size, assuming uniform mesh
dx = grids.x_u(2)-grids.x_u(1);
dy = grids.y_u(2)-grids.y_u(1);

%% Create the cell type matrices
u_type = zeros(N.y_u,N.x_u);
v_type = zeros(N.y_v,N.x_v);
p_type = zeros(N.y_p,N.x_p);

%% Set the u-cell types
u_type(:,end)       = 2;            %%Right/East Edge 
u_type(1,2:end-1)   = 4;            %%Bottom/South Edge
u_type(:,1)         = 6;            %%Left/West Edge
u_type(end,2:end-1) = 8;            %%Top/North Edge

%% Set the v-cell type
v_type(2:end-1,end) = 2;            %%Right/East Edge 
v_type(1,:)         = 4;            %%Bottom/South Edge
v_type(2:end-1,1)   = 6;            %%Left/West Edge 
v_type(end,:)       = 8;            %%Top/North Edge

%% Set the p-cell type
p_type(end,end) = 1;        %Northeast Corner
p_type(2:end-1,end) = 2;    %East Face
p_type(1,end) = 3;          %Southeast Corner
p_type(1,2:end-1) = 4;      %South Face
p_type(1,1) = 5;            %Southwest Corner
p_type(2:end-1,1) = 6;      %West Face
p_type(end,1) = 7;          %Northwest Corner
p_type(end,2:end-1) = 8;    %North Face

plot_cell_type
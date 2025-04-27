function [A_p, b_p, Ap_u, Ap_v] = A_p_creation(u_star, v_star, p_guess, A_u, A_v, ...
    dx, dy, p_type, A_p)

%% I have no clue why, but the if statements should be 
% p_type(i) == 5 || p_type(i) == 11
% p_type(i) == 7 || p_type(i) == 13
% p_type(i) == 6 || p_type(i) == 12
% p_type(i) == 3 || p_type(i) == 17
% p_type(i) == 1 || p_type(i) == 15
% p_type(i) == 2 || p_type(i) == 16
% p_type(i) == 4 || p_type(i) == 10
% p_type(i) == 8 || p_type(i) == 14
% but when these are used, the matrix becomes singular. It seems to work
% without though.


[Ny_u, Nx_u] = size(u_star);
[Ny_v, Nx_v] = size(v_star);
[Ny_p, Nx_p] = size(p_guess);

 %% Get the velocity at the midpoints of the Pressure cell

uw_star = u_star(:,1:end-1);
ue_star = u_star(:,2:end);

vn_star = v_star(2:end,:);
vs_star = v_star(1:end-1,:);

%% Set the value of aP for the u and v velocities

Ap_u = reshape(diag(A_u),Ny_u,Nx_u);    

Ap_v = reshape(diag(A_v),Ny_v,Nx_v);

%% Calculate the pressure correction matrix coefficients

Cw = 1./Ap_u(:,1:end-1)/dx^2;
Ce = 1./Ap_u(:,2:end)/dx^2;
Cn = 1./Ap_v(2:end,:)/dy^2;
Cs = 1./Ap_v(1:end-1,:)/dy^2;

%% Create the b vector

b_p = zeros(length(p_guess(:)),1);

%% Set the pressure correction coefficients

for i = 1:Ny_p*Nx_p      

    % Boundaries on the west and south face
    if p_type(i) == 5
        A_p(i,i) = Ce(i) + Cn(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i+1)    = -Cn(i);

    % Boundaries on the west and north face
    elseif p_type(i) == 7
        A_p(i,i) = Ce(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i-1)    = -Cs(i);

    % Boundaries on the west face
    elseif p_type(i) == 6        
        A_p(i,i) = Ce(i) + Cn(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);

    % Boundaries on the east and south face
    elseif p_type(i) == 3
        A_p(i,i) = Cw(i) + Cn(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-Ny_p) = -Cw(i);

    % Boundaries on the east and north face
    elseif p_type(i) == 1
        A_p(i,i) = Cw(i) + Cs(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw(i);

    % Boundaries on the east face
    elseif p_type(i) == 2      
        A_p(i,i) = Cw(i) + Cn(i) + Cs(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw(i);

    % Boundaries on the south face
    elseif p_type(i) == 4
        A_p(i,i) = Cw(i) + Ce(i) + Cn(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-Ny_p) = -Cw(i);            

    % Boundaries on the north face
    elseif p_type(i) == 8
        A_p(i,i) = Cw(i) + Ce(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw(i);

    elseif p_type(i) == -1
        A_p(i,i) = 1;     
        
    % Interior cells
    else
        A_p(i,i) = Cw(i) + Ce(i) + Cn(i) + Cs(i);
        A_p(i,i+Ny_p) = -Ce(i);
        A_p(i,i+1) = -Cn(i);
        A_p(i,i-1) = -Cs(i);
        A_p(i,i-Ny_p) = -Cw(i);            
    end            

    if p_type(i) == -1
        b_p(i) = 0;
    else
        b_p(i) = -(ue_star(i)-uw_star(i))/dx -(vn_star(i)-vs_star(i))/dy;
    end
    
    
end

%% Set the correction of the pressure at the middle of the domain to zero

ind = Ny_p*Nx_p; 

A_p(ind,:) = 0;
A_p(ind,ind) = 1;
b_p(ind) = 0;



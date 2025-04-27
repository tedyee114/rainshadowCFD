function [A_u, b_u] = A_u_creation_efficient(u_guess, v_guess, p_guess, dx, dy, dt, Re, u_type, u_prevTime,  A_u, grids, speeds, BC, noslipground);

[Ny_u, Nx_u] = size(u_guess);

u_P = u_guess;
u_E = [u_guess(:,2:end) zeros(Ny_u,1)];
u_W = [zeros(Ny_u,1) u_guess(:,1:end-1)];

v_ne = [v_guess(2:end,:) zeros(Ny_u,1)];
v_nw = [zeros(Ny_u,1) v_guess(2:end,:)];
v_se = [v_guess(1:end-1,:) zeros(Ny_u,1)];
v_sw = [zeros(Ny_u,1) v_guess(1:end-1,:)];

%%Calculation of the velocity at the midpoints (nonlinear terms)
ue = (u_P + u_E)/2;
uw = (u_P + u_W)/2;
vn = (v_ne + v_nw)/2;
vs = (v_se + v_sw)/2;

%%Calculate the pressures at the east and west midpoint
p_e = [p_guess zeros(Ny_u,1)];
p_w = [zeros(Ny_u,1) p_guess];

%%Coefficients for the A matrix
Fe_u = ue/2/dx; Fw_u = uw/2/dx; Fn_u = vn/2/dy; Fs_u = vs/2/dy;
De = 1/dx^2/Re; Dw = 1/dx^2/Re; Dn = 1/dy^2/Re; Ds = 1/dy^2/Re;        

%%Create the b vector
b_u = zeros(length(u_guess(:)),1);
A_u = sparse(Ny_u * Nx_u, Ny_u * Nx_u);

index=ones(1,4);    %%% indices for scrolling through nonuniform BC

%% Update Matrices by cell type
for i = 1:length(u_type(:))
    %%Right/East Edge
    if u_type(i) == 2 
            A_u(i,i) = 1;  
            A_u(i,i-1) = -1;
            b_u(i) =  BC.u_rig_value(index(1));

    
    %%Bottom/South Edge
    elseif u_type(i) == 4
        if BC.u_bot_type == "Dirichlet"
            A_u(i,i)      =  1;
            b_u(i) = BC.u_bot_value(index(2));
            index(2)=index(2)+1;
        elseif BC.u_bot_type == "Neumann"
            A_u(i,i) = 1;  
            A_u(i,i-1) = -1;
            b_u(i) =  BC.u_bot_value(index(2));
            index(2)=index(2)+1;
        else                   %%% for free-slip
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+Dw+Dn+2*Ds + 1/dt;            
            A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
            A_u(i,i+1)    =  Fn_u(i)-Dn;  
            % A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i)        = -(p_e(i)-p_w(i))/dx + 2*speeds.u_bot*Ds + u_prevTime(i)/dt;
        end

    %%Left/West Edge
    elseif u_type(i) == 6
        if BC.u_lef_type == "Dirichlet"
            A_u(i,i)      =  1;
            b_u(i) = BC.u_lef_value(index(3));
            index(3)=index(3)+1;
        elseif BC.u_lef_type == "Neumann"
            A_u(i,i) = 1;  
            A_u(i,i-1) = -1;
            b_u(i) =  BC.u_lef_value(index(3));
            index(3)=index(3)+1;
        else                   %%% for free-slip
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+2*Dw+Dn+Ds + 1/dt;            
            % A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
            A_u(i,i+1)    =  Fn_u(i)-Dn;  
            A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i)        = -(p_e(i)-p_w(i))/dx + 2*speeds.u_lef*Dw + u_prevTime(i)/dt;
        end

    %%Top/North Edge
    elseif u_type(i) == 8
        if BC.u_top_type == "Dirichlet"
            A_u(i,i)      =  1;
            b_u(i) = BC.u_top_value(index(4));
            index(4)=index(4)+1;
        elseif BC.u_top_type == "Neumann"
            A_u(i,i) = 1;  
            A_u(i,i-1) = -1;
            b_u(i) =  BC.u_top_value(index(4));
            index(4)=index(4)+1;
        else                   %%% for free-slip
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+Dw+2*Dn+Ds + 1/dt;
                        A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
            % A_u(i,i+1)    =  Fn_u(i)-Dn;  
            A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i)        = -(p_e(i)-p_w(i))/dx + 2*speeds.u_top*Dn + u_prevTime(i)/dt;
        end
    
    %%If the cell has a bottom face on the obstruction
    elseif u_type(i) == 10
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+Dw+Dn+2*Ds + 1/dt;
            A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
            A_u(i,i+1)    =  Fn_u(i)-Dn;  
            % A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt;
        end

    %%If the cell has a bottom and a west face on the obstruction
    elseif u_type(i) == 11
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+2*Dw+Dn+2*Ds + 1/dt;
                % A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
            A_u(i,i+1)    =  Fn_u(i)-Dn;  
                % A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt;
        end

    %%If the cell has a west face on the obstruction
    elseif u_type(i) == 12
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+2*Dw+Dn+Ds + 1/dt;
                % A_u(i,i-Ny_u) = -Fw_u(i)-Dw; 
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
            A_u(i,i+1)    =  Fn_u(i)-Dn;  
            A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt;
        end

    %%If the cell has a top and an west face on the obstruction
    elseif u_type(i) == 13
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+2*Dw+2*Dn+Ds + 1/dt;
                %A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
                % A_u(i,i+1)    =  Fn_u(i)-Dn;  
            A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt;
        end

    %%If the cell has a top face on the obstruction
    elseif u_type(i) == 14
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+Dw+2*Dn+Ds + 1/dt;
            A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
            A_u(i,i+Ny_u) =  Fe_u(i)-De;
                % A_u(i,i+1)    =  Fn_u(i)-Dn;  
            A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt;
        end

    %%If the cell has a top and an east face on the obstruction
    elseif u_type(i) == 15
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      = Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+2*De+Dw+2*Dn+Ds + 1/dt;
            A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
                % A_u(i,i+Ny_u) =  Fe_u(i)-De;
                % A_u(i,i+1)    =  Fn_u(i)-Dn;  
            A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt; 
        end
    
    %%If the cell has an east face on the obstruction
    elseif u_type(i) == 16
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+2*De+Dw+Dn+Ds + 1/dt;
            A_u(i,i+1)    =  Fn_u(i)-Dn;  
            A_u(i,i-1)    = -Fs_u(i)-Ds; 
            A_u(i,i-Ny_u) = -Fw_u(i)-Dw; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt; 
        end

    %%If the cell has a bottom and east face on the obstruction
    elseif u_type(i) == 17
        if noslipground == 1
            A_u(i,i)      =  1;
            b_u(i) = 0;
        else
            A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+2*De+Dw+Dn+2*Ds + 1/dt;
            A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
                % A_u(i,i+Ny_u) =  Fe_u(i)-De;
            A_u(i,i+1)    =  Fn_u(i)-Dn;  
                % A_u(i,i-1)    = -Fs_u(i)-Ds; 
            b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt; 
        end
            
    %%Obstruction Cells
    elseif u_type(i) == -1
        A_u(i,i) = 1;  
        b_u(i) = 0;      
        
    %%Interior Cells
    else
        A_u(i,i)      =  Fe_u(i)-Fw_u(i)+Fn_u(i)-Fs_u(i)+De+Dw+Dn+Ds + 1/dt;
        A_u(i,i-Ny_u) = -Fw_u(i)-Dw;             
        A_u(i,i+Ny_u) =  Fe_u(i)-De;
        A_u(i,i+1)    =  Fn_u(i)-Dn;  
        A_u(i,i-1)    = -Fs_u(i)-Ds; 
        b_u(i) = -(p_e(i)-p_w(i))/dx + u_prevTime(i)/dt;              
    end
end
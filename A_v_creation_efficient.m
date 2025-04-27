function [A_v, b_v] = A_v_creation_efficient(u_guess, v_guess, p_guess, dx, dy, dt, Re, v_type, v_prevTime, A_v, grids, speeds, BC, noslipground); 
       
[Ny_v, Nx_v] = size(v_guess);

%%Calculation of the velocity at the midpoints (nonlinear terms)
u_ne = [u_guess(:,2:end); zeros(1,Nx_v)];
u_se = [zeros(1,Nx_v); u_guess(:,2:end)];
u_nw = [u_guess(:,1:end-1); zeros(1,Nx_v)];
u_sw = [zeros(1,Nx_v); u_guess(:,1:end-1)];

v_P = v_guess;
v_N = [v_guess(2:end,:); zeros(1,Nx_v)];
v_S = [zeros(1,Nx_v); v_guess(1:end-1,:)];

ue = (u_ne+u_se)/2;
uw = (u_nw+u_sw)/2;
vn = (v_P+v_N)/2;
vs = (v_P+v_S)/2;


%%Calculate of pressure at the north and south midpoints
p_n = [zeros(1,Nx_v); p_guess(2:end,:); zeros(1,Nx_v)];
p_s = [zeros(1,Nx_v); p_guess(1:end-1,:); zeros(1,Nx_v)];

%%Coefficients for the A matrix
Fe_v = ue/2/dx; Fw_v = uw/2/dx; Fn_v = vn/2/dy; Fs_v = vs/2/dy;
De = 1/dx^2/Re; Dw = 1/dx^2/Re; Dn = 1/dy^2/Re; Ds = 1/dy^2/Re;

%%Create the b vector for the vertical velocity calculation
b_v = zeros(length(v_guess(:)),1);
A_v = sparse(Ny_v * Nx_v, Ny_v * Nx_v);

index=ones(1,4);    %%% indices for scrolling through nonuniform BC

%% Update Matrices by cell type
for i = 1:Ny_v*Nx_v
    %%Right/East Edge
    if v_type(i) == 2
        if BC.v_rig_type == "Dirichlet"
            A_v(i,i)      =  1;
            b_v(i) = BC.v_rig_value(index(1));
            index(1)=index(1)+1;
        elseif BC.v_rig_type == "Neumann"
            A_v(i,i) = 1;  
            A_v(i,i-1) = -1;
            b_v(i) =  BC.v_rig_value(index(1));
            index(1)=index(1)+1;
        else                   %%% for free-slip
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+2*De+Dw+Dn+Ds + 1/dt;            
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            % A_v(i,i+Ny_v) =  Fe_v(i)-De;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i)        = -(p_n(i)-p_s(i))/dy + 2*speeds.v_rig*De+ v_prevTime(i)/dt;   
        end


    %%Bottom/South Edge
    elseif v_type(i) == 4
        if BC.v_bot_type == "Dirichlet"
            A_v(i,i)      =  1;
            b_v(i) = BC.v_bot_value(index(2));
            index(2)=index(2)+1;
        elseif BC.v_bot_type == "Neumann"
            A_v(i,i) = 1;  
            A_v(i,i-1) = -1;
            b_v(i) =  BC.v_bot_value(index(2));
            index(2)=index(2)+1;
        else                   %%% for free-slip
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+Dw+Dn+2*Ds + 1/dt;            
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
            % A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i)        = -(p_n(i)-p_s(i))/dy + 2*speeds.v_bot*Ds + v_prevTime(i)/dt;
        end

    %%Left/West Edge
    elseif v_type(i) == 6
        if BC.v_lef_type == "Dirichlet"
            A_v(i,i)      =  1;
            b_v(i) = BC.v_lef_value(index(3));
            index(3)=index(3)+1;
        elseif BC.v_lef_type == "Neumann"
            A_v(i,i) = 1;  
            A_v(i,i-1) = -1;
            b_v(i) =  BC.v_lef_value(index(3));
            index(3)=index(3)+1;
        else                   %%% for free-slip
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+2*Dw+Dn+Ds + 1/dt;            
            % A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i)        = -(p_n(i)-p_s(i))/dy + 2*speeds.v_lef*Dw + v_prevTime(i)/dt;
        end

    %%Top/North Edge
    elseif v_type(i) == 8
        if BC.v_top_type == "Dirichlet"
            A_v(i,i)      =  1;
            b_v(i) = BC.v_top_value(index(4));
            index(4)=index(4)+1;
        elseif BC.v_top_type == "Neumann"
            A_v(i,i) = 1;  
            A_v(i,i-1) = -1;
            b_v(i) =  BC.v_top_value(index(4));
            index(4)=index(4)+1;
        else                   %%% for free-slip
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+Dw+2*Dn+Ds + 1/dt;
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
            % A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i)        = -(p_n(i)-p_s(i))/dy + 2*speeds.v_top*Dn + v_prevTime(i)/dt;
        end
    
    %%If the cell has a bottom face on the obstruction
    elseif v_type(i) == 10
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+Dw+Dn+2*Ds + 1/dt;
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
            % A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt;              
        end

    %%If the cell has a bottom and a west face on the obstruction
    elseif v_type(i) == 11
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+2*Dw+Dn+2*Ds + 1/dt;
                % A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
                % A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt;
        end

    %%If the cell has a west face on the obstruction
    elseif v_type(i) == 12
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+2*Dw+Dn+Ds + 1/dt;
                % A_v(i,i-Ny_v) = -Fw_v(i)-Dw; 
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt; 
        end

    %%If the cell has a top and an west face on the obstruction
    elseif v_type(i) == 13
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+2*Dw+2*Dn+Ds + 1/dt;
                %A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
                % A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt; 
        end

    %%If the cell has a top face on the obstruction
    elseif v_type(i) == 14
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+Dw+2*Dn+Ds + 1/dt;
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
            A_v(i,i+Ny_v) =  Fe_v(i)-De;
                % A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt;
        end

    %%If the cell has a top and an east face on the obstruction
    elseif v_type(i) == 15
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      = Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+2*De+Dw+2*Dn+Ds + 1/dt;
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
                % A_v(i,i+Ny_v) =  Fe_v(i)-De;
                % A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt;  
        end
        
    %%If the cell has an east face on the obstruction
    elseif v_type(i) == 16
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+2*De+Dw+Dn+Ds + 1/dt;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
            A_v(i,i-1)    = -Fs_v(i)-Ds; 
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt; 
        end

    %%If the cell has a bottom and east face on the obstruction
    elseif v_type(i) == 17
        if noslipground == 1
            A_v(i,i)      =  1;
            b_v(i) = 0;
        else
            A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+2*De+Dw+Dn+2*Ds + 1/dt;
            A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
                % A_v(i,i+Ny_v) =  Fe_v(i)-De;
            A_v(i,i+1)    =  Fn_v(i)-Dn;  
                % A_v(i,i-1)    = -Fs_v(i)-Ds; 
            b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt;  
        end
        
    %%Obstruction Cells
    elseif v_type(i) == -1
        A_v(i,i) = 1;  
        b_v(i) = 0;      
        
    %%Interior Cells
    else
        A_v(i,i)      =  Fe_v(i)-Fw_v(i)+Fn_v(i)-Fs_v(i)+De+Dw+Dn+Ds + 1/dt;
        A_v(i,i-Ny_v) = -Fw_v(i)-Dw;             
        A_v(i,i+Ny_v) =  Fe_v(i)-De;
        A_v(i,i+1)    =  Fn_v(i)-Dn;  
        A_v(i,i-1)    = -Fs_v(i)-Ds; 
        b_v(i) = -(p_n(i)-p_s(i))/dy + v_prevTime(i)/dt;              
    end
end
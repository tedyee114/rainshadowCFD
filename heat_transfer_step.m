%% Perform a step of the heat transfer analysis
phi_n = phi_np1;
A_phi = zeros(length(phi_n(:)));
b_phi = zeros(length(phi_n(:)),1);


%% Get the velocity at the midpoints of the Pressure/Temperature cell
uw = u_star(:,1:end-1);
ue = u_star(:,2:end);

vn = v_star(2:end,:);
vs = v_star(1:end-1,:);

%% Coefficients for the A matrix
Fe = ue/2/dx;   Fw = uw/2/dx;   Fn = vn/2/dy;   Fs = vs/2/dy;
De = 1/dx^2/Re; Dw = 1/dx^2/Re; Dn = 1/dy^2/Re; Ds = 1/dy^2/Re;

%% Create the A matrix for the phi system
index = 1;
for i = 1:N.y_p*N.x_p      
    

    if p_type(i) == 5           %%Southwest Corner
        A_phi(i,i) = Fe(i)+Fn(i)-2*Fs(i) + De+2*Dw+Dn + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i+1)    = Fn(i)-Dn;
        b_phi(i) = phi_n(i)/dt + 2*(Fw(i)+Dw)*phi_in(index);
        index = index + 1;

    elseif p_type(i) == 7       %%Northwest Corner
        A_phi(i,i) = Fe(i)-Fs(i) + De+2*Dw+2*Dn+Ds + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        b_phi(i) = phi_n(i)/dt + 2*(Fw(i)+Dw)*phi_in(index) + 2*(-Fn(i)+Dn)*phi_wall;
        index = index + 1;

    elseif p_type(i) == 6      %%Left/West Edge 
        A_phi(i,i) = Fe(i)+Fn(i)-Fs(i) + De+2*Dw+Dn+Ds + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        b_phi(i) = phi_n(i)/dt + 2*(Fw(i)+Dw)*phi_in(index);
        index = index + 1;

    elseif p_type(i) == 11           %%Southwest Corner
        A_phi(i,i) = Fe(i)+Fn(i)-2*Fs(i) + De+2*Dw+Dn + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i+1)    = Fn(i)-Dn;
        b_phi(i) = phi_n(i)/dt;

    elseif p_type(i) == 13       %%Northwest Corner
        A_phi(i,i) = Fe(i)-Fs(i) + De+2*Dw+2*Dn+Ds + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        b_phi(i) = phi_n(i)/dt;

    elseif p_type(i) == 12      %%Left/West Edge 
        A_phi(i,i) = Fe(i)+Fn(i)-Fs(i) + De+2*Dw+Dn+Ds + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        b_phi(i) = phi_n(i)/dt;
    
    elseif p_type(i) == 3 || p_type(i) == 17      %%Southeast Corner
        A_phi(i,i) = 2*Fe(i)-Fw(i)+Fn(i)-2*Fs(i) +Dw+Dn + 1/dt;        
        A_phi(i,i+1)    = Fn(i)-Dn;        
        A_phi(i,i-N.y_p) =-Fw(i)-Dw;
        b_phi(i) = phi_n(i)/dt;

    elseif p_type(i) == 1 || p_type(i) == 15      %%Northeast Corner
        A_phi(i,i) = 2*Fe(i)-Fw(i)-Fs(i) + Dw+Ds+2*Dn + 1/dt;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-N.y_p) =-Fw(i)-Dw;
        b_phi(i) = phi_n(i)/dt + 2*(-Fn(i)+Dn)*phi_wall;

    elseif p_type(i) == 2 || p_type(i) == 16      %%Right/East Edge
        A_phi(i,i) = 2*Fe(i)-Fw(i)+Fn(i)-Fs(i) + +Dw+Dn+Ds + 1/dt;        
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-N.y_p) =-Fw(i)-Dw;
        b_phi(i) = phi_n(i)/dt;

    elseif p_type(i) == 4 || p_type(i) == 10      %%Bottom/South Edge
        A_phi(i,i) = Fe(i)-Fw(i)+Fn(i)-2*Fs(i) + De+Dw+Dn + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-N.y_p) =-Fw(i)-Dw;
        b_phi(i) = phi_n(i)/dt;           

    elseif p_type(i) == 8 || p_type(i) == 14      %%Top/North Edge
        A_phi(i,i) = Fe(i)-Fw(i)   -Fs(i) + De+Dw+2*Dn+Ds + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-N.y_p) =-Fw(i)-Dw;
        b_phi(i) = phi_n(i)/dt + 2*(-Fn(i)+Dn)*phi_wall;

    elseif p_type(i) == -1    %%Obstruction Cells
        A_phi(i,i) = 1;  
        b_phi(i) = 0;      
        
    % Interior cells
    else
        A_phi(i,i) = Fe(i)-Fw(i)+Fn(i)-Fs(i) + De+Dw+Dn+Ds + 1/dt;
        A_phi(i,i+N.y_p) = Fe(i)-De;
        A_phi(i,i+1)    = Fn(i)-Dn;
        A_phi(i,i-1)    =-Fs(i)-Ds;
        A_phi(i,i-N.y_p) =-Fw(i)-Dw;
        b_phi(i) = phi_n(i)/dt;
    end             
    
end

phi_np1 = A_phi\b_phi;
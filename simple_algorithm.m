%% Perform simple algorithm to calculate u^n+1, v^n+1, and p^n+1

while II <= II_max && abs(u_change(end)) > u_tol
    u_guess_old = u_guess;
    v_guess_old = v_guess;
    p_guess_old = p_guess;    

    %%Calculation of u-star %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A_u, b_u] = A_u_creation_efficient(u_guess, v_guess, p_guess, dx, dy, dt, Re, u_type, u_prevTime,  A_u, grids, speeds, BC, noslipground);
    u_star = reshape(A_u\b_u,N.y_u,N.x_u);              %%Calculate the u star value
    u_star = u_guess + alpha_u*(u_star-u_guess);        %%Update the u_star value with relaxation

    %%Calculation of v-star %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A_v, b_v] = A_v_creation_efficient(u_guess, v_guess, p_guess, dx, dy, dt, Re, v_type, v_prevTime, A_v, grids, speeds, BC, noslipground); 
    v_star = reshape(A_v\b_v,N.y_v,N.x_v);              %%Calculate the v-star value
    v_star = v_guess + alpha_v*(v_star-v_guess);        %%Use relaxation to set the v-star value

    %%Calculation of the pressure correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A_p, b_p, Ap_u, Ap_v] = A_p_creation(u_star, v_star, p_guess, A_u, A_v, dx, dy, p_type, A_p);
    p_correction = reshape(A_p\b_p,N.y_p,N.x_p);        %%Calculate the pressure correction and apply it with under-relaxation


    %% Calculation of the velocity correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%Calculate the pressure correction on the east and west faces of u-cells
    pc_E = p_correction(:,2:end);
    pc_W = [p_correction(:,1:end-1)];

    %%Calculate the velocity correction for all cells IN the domain
    u_correction = 0*u_guess;
    u_correction(:,2:end-1) = (pc_W-pc_E)/dx./Ap_u(:,2:end-1);
    u_correction(u_type==-1) = 0;

    %%Calculate the pressure correction on the north and south faces of v-cells
    pc_N = [p_correction(2:end,:)];
    pc_S = [p_correction(1:end-1,:)];

    %%Calculate the velocity correction for all cells IN the domain
    v_correction = 0*v_guess;
    v_correction(2:end-1,:) = (pc_S-pc_N)/dy./Ap_v(2:end-1,:);

    % Set the correction for cells with a center ON the boundary
    v_correction(v_type==4) = 0;
    v_correction(v_type==8) = 0;
    v_correction(v_type==-1) = 0;
    
    %% Corect the velocity directly and pressure with under-relaxation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_guess = u_star + u_correction;
    v_guess = v_star + v_correction;
    p_guess = p_guess + alpha_p*p_correction;

    %%Impose flux conservation
    m_in = sum(u_guess(:,1)*dy);
    m_out = sum(u_guess(:,end-1)*dy);    %%%calculated at the second to last column

    u_guess(:,end) = m_in/m_out*u_guess(:,end-1);

    %% Calculate the continuity residual and changes in u,v, and p %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ue = u_guess(:,2:end);
    uw = u_guess(:,1:end-1);
    vn = v_guess(2:end,1:end);
    vs = v_guess(1:end-1,1:end);

    continuity_residual(II) = sqrt(sum(((ue(:)-uw(:))/dx+(vn(:)-vs(:))/dy).^2));

    u_change(II) = sqrt(sum((u_guess(:)-u_guess_old(:)).^2))/N.x_u/N.y_u;
    v_change(II) = sqrt(sum((v_guess(:)-v_guess_old(:)).^2))/N.x_v/N.y_v;
    p_change(II) = sqrt(sum((p_guess(:)-p_guess_old(:)).^2))/N.x_p/N.y_p;

    %% Plot the continuity residual %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     close all;
    if n == 1
        if mod(II,20) == 0
            ax4 = axes('Parent', panel_2, 'Position', [0.05, 0.95, 0.9, 0.3]);
            imagesc(ax4, flipud(u_guess)); colorbar; caxis([-1 1]*max(abs(u_guess(:))));
            title('Horizontal (u) Velocity')
            colorbar
            set(gca,'YDir','normal')
            
            ax5 = axes('Parent', panel_2, 'Position', [0.05, 0.5, 0.9, 0.3]);
            imagesc(ax5, flipud(v_guess)); colorbar; caxis([-1 1]*max(abs(v_guess(:))));
            title('Vertical (v) Velocity')
            colorbar
            set(gca,'YDir','normal')
            
            ax6 = axes('Parent', panel_2, 'Position', [0.05, 0.1, 0.9, 0.3]);
            imagesc(ax6, flipud(p_guess)); colorbar; caxis([-1 1]*max(abs(p_guess(:))));
            title('Pressure')
            set(gca,'YDir','normal')
            colorbar
            colormap jet
            drawnow
        end
    end
    II = II + 1;
    
end
% Left Side of Mountain
for i = x_mountain_start+1 : x_mountain_peak
    j = round(rise_slope * (i - x_mountain_start));
    if j <= N.y_u
        u_type(1:j, i) = -1;
    end
    if j <= N.y_v
        v_type(1:j, i) = -1;
    end
    if j <= N.y_p
        p_type(1:j, i) = -1;
    end
end

% Right Side Of Mountain
for i = x_mountain_peak+1 : x_mountain_end
    j = round(fall_slope * (i - x_mountain_peak) + j);  % continue from where rise ended
    if j >= 1 && j <= N.y_u
        u_type(1:j, i) = -1;
    end
    if j >= 1 && j <= N.y_v
        v_type(1:j, i) = -1;
    end
    if j >= 1 && j <= N.y_p
        p_type(1:j, i) = -1;
    end
end

%%%Boundary cells u-type
for i = 1:N.y_u-1
    for j = 2:N.x_u-1
        if u_type(i,j) == 0
            below = u_type(i-1,  j);
            above  = u_type(i+1, j);
            right = u_type(i,    j+1);
            left  = u_type(i,    j-1);            
            
            if above == -1
                u_type(i,j) = 14;
            elseif below == -1
                u_type(i,j) = 10;
            elseif right == -1
                u_type(i,j) = 16;
            elseif left == -1               %%left and right overwrite at corners
                u_type(i,j) = 12;
            
            end
        end
    end
end

%%%Boundary cells v-type
for i = 1:N.y_v-1
    for j = 2:N.x_v-1
        if v_type(i,j) == 0
            below = v_type(i-1,  j);
            above  = v_type(i+1, j);
            right = v_type(i,    j+1);
            left  = v_type(i,    j-1);            
            
            if below == -1 && left == -1
                v_type(i,j) = 11;
            elseif left == -1
                v_type(i,j) = 12;
            elseif above == -1 && left == -1
                v_type(i,j) = 13;    
            elseif above == -1
                v_type(i,j) = 14;
            elseif above == -1 && right == -1
                v_type(i,j) = 15;
            elseif right == -1
                v_type(i,j) = 16;
            elseif below == -1 && right == -1
                v_type(i,j) = 17;
            elseif below == -1
                v_type(i,j) = 10;
            end
        end
    end
end

%%%Boundary cells p-type
for i = 1:N.y_p-1
    for j = 2:N.x_p-1
        if p_type(i,j) == 0
            below = p_type(i-1,  j);
            above  = p_type(i+1, j);
            right = p_type(i,    j+1);
            left  = p_type(i,    j-1);            
            
            if below == -1 && left == -1
                p_type(i,j) = 11;
            elseif left == -1
                p_type(i,j) = 12;
            elseif above == -1 && left == -1
                p_type(i,j) = 13;    
            elseif above == -1
                p_type(i,j) = 14;
            elseif above == -1 && right == -1
                p_type(i,j) = 15;
            elseif right == -1
                p_type(i,j) = 16;
            elseif below == -1 && right == -1
                p_type(i,j) = 17;
            elseif below == -1
                p_type(i,j) = 10;
            end
        end
    end
end
% Rain Calculations
rainmap = zeros(length(cloudmap),1);
cloud_density_before_rain = zeros(length(cloudmap),1);


for i = 1:length(cloudmap)
    % Initialize total moisture
    total_moisture = 0;
    if cloudmap(i) == 1
        total_moisture = phi_np1(i);  
    end
    % Calculate row and column from 1D index
    row = ceil(i / N.y_p);
    col = mod(i - 1, N.y_p) + 1;
    
    % Check left neighbor
    if col > 1 % Not at left edge
        left_idx = i - 1;
        if cloudmap(left_idx) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(left_idx);
        end
    else % At left edge, use own value
        if cloudmap(i) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(i);
        end
    end
    
    % Check right neighbor
    if col < N.y_p % Not at right edge
        right_idx = i + 1;
        if cloudmap(right_idx) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(right_idx);
        end
    else % At right edge, use own value
        if cloudmap(i) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(i);
        end
    end
    
    % Check above neighbor
    if row > 1 % Not at top edge
        above_idx = i - N.y_p;
        if cloudmap(above_idx) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(above_idx);
        end
    else % At top edge, use own value
        if cloudmap(i) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(i);
        end
    end
    
    % Check below neighbor
    if row < N.x_p % Not at bottom edge
        below_idx = i + N.y_p;
        if cloudmap(below_idx) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(below_idx);
        end
    else % At bottom edge, use own value
        if cloudmap(i) == 1 % If it's a cloud
            total_moisture = total_moisture + phi_np1(i);
        end
    end
    
    cloud_density_before_rain(i) = total_moisture;
    
    % Check if critical mass is reached
    if total_moisture > criticalmass
        for k = 1:N.y_p
            row_idx = (row - 1) * N.y_p + k + round(rainremovalrate*rand);
            if row_idx <= length(cloudmap)                  % in case the random extra factor goes outside of bounds
                rainmap(row_idx) = 2; % Set each element in the row to 2
                                
                % Add to rain gauge for this column the same amount that was just cleared
                rain_gauge(row) = rain_gauge(row) + phi_np1(row_idx);
            end
        end
    end
end

phi_np1(rainmap==2) = 0; % Reset moisture to 0 where rain is occurring at this timeframe

% Reshape and combine maps
rainmap_reshaped = reshape(rainmap, N.y_p, N.x_p);
rainmapwithclouds = rainmap_reshaped; cloudmap_mask = (cloudmap == 1); rainmapwithclouds(cloudmap_mask) = cloudmap(cloudmap_mask);

% Draw the mountain on the rainmap
p_type_reshaped = reshape(p_type, N.y_p, N.x_p);
p_type_mask = (p_type_reshaped == -1);
rainmapwithclouds(p_type_mask) = -1;
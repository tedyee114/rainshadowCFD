%% Clouds
[~, altitude] = meshgrid(grids.x_p, grids.y_p);
TK = T_surface - lapse_rate * altitude;
TK = TK(:);
TC = TK - 273.15;
Volume = dx * dy*scalefactor;
partial_pressure = phi_np1 .* 461.5 .* TK / Volume;
satvap_pressure = 610.94 * exp(17.625 * TC ./ (TC + 243.04));   % August-Roche-Magnus equation
cloudmap = partial_pressure >= satvap_pressure;                 % Array of 0 if evaporation and 1 if condensation/cloud


%% Plotting before Rain
panel_3 = figure('Name',['Timestep ', num2str(n),', which is about ',num2str(n*dt*10), 'hrs']);
ax7 = tiledlayout(panel_3,2,2);


%% Rain
rain

%% Plotting After Rain
nexttile
imagesc(grids.x_p, grids.y_p, reshape(cloud_density_before_rain, N.y_p, N.x_p))
set(gca, 'YDir','normal')
title('Cloud Density Before Rain (Cumulative Area Average)')
colorbar

nexttile(ax7)
imagesc(grids.x_p, grids.y_p, rainmapwithclouds)
set(gca, 'YDir','normal')
title('rainmapwithclouds')
set(gca, 'YDir','normal')
raincolor = [0 0 1];
cloudcolor = [1 1 1];   
skycolor = [0 195/255 235/255];
mountaincolor = [125/255 162/255 126/255];
colormap(gca, [mountaincolor; skycolor; cloudcolor; raincolor])
set(gca, 'CLim', [-1 2]);
title(['Realistic Weather Rendering at Timestep ', num2str(n),', ',num2str(n*dt*10), 'hrs'])

nexttile(ax7)
imagesc(grids.x_p, grids.y_p, reshape(phi_np1, N.y_p, N.x_p))
set(gca, 'YDir','normal')
title('Moisture content after rain')
colorbar

nexttile(ax7)
plot(grids.x_p, rain_gauge, 'b-', 'LineWidth', 2)
title('Rain Gauge - Total Rain Fallen')
xlabel('x-location on Earth')
ylabel('Cumulative Rainfall')
grid on
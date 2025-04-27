% Create UI as a figure
base = uifigure('Name', 'ANTSSSSS Solver Setup Window', "Icon","tyicon-transparent.png", 'color', [.3,.3,.3]);


%% Inputs

fig = uigridlayout(base,[20,5]);

%%%Title
welcome = uilabel(fig);

welcome.Text = "<font style='font-weight: bold; color: orange; font-family: sans-serif; font-size: 20;'>Ʌ</font>" + ...
               "<font style='font-weight: bold; color: black; font-family: sans-serif; font-size: 20;'>NSYS Lite, aka ANTSSSSS. This window is for setting up your solver methods now that the geometry has been set.</font>";
welcome.Interpreter = "html";
welcome.Layout.Row = 1; 
welcome.Layout.Column = [1,4];


%%Left u
u_lef = uilabel(fig, 'Text', 'u_left');
u_lef_type_dropdown = uidropdown(fig, "Items", ["Dirichlet", "FreeSlip", "Neumann"]);
u_lef_value_box = create_textbox(fig, 'Value:', 'cos(8*pi*y)', [30, 60, 120, 22]);

u_lef.Layout.Row = 3; 
u_lef.Layout.Column = 1;
u_lef_type_dropdown.Layout.Row = 4; 
u_lef_type_dropdown.Layout.Column = 1;
u_lef_value_box.Layout.Row = 4; 
u_lef_value_box.Layout.Column = 2;

%%Right u
u_rig = uilabel(fig, 'Text', 'u_right');
u_rig_type_dropdown = uidropdown(fig, "Items", ["Neumann", "FreeSlip", "Dirichlet"]);
u_rig_value_box = create_textbox(fig, 'Value:', '0', [30, 60, 120, 22]);

u_rig.Layout.Row = 5; 
u_rig.Layout.Column = 1;
u_rig_type_dropdown.Layout.Row = 6; 
u_rig_type_dropdown.Layout.Column = 1;
u_rig_value_box.Layout.Row = 6; 
u_rig_value_box.Layout.Column = 2;

%Bottom u
u_bot = uilabel(fig, 'Text', 'u_bottom');
u_bot_type_dropdown = uidropdown(fig, "Items", [ "FreeSlip", "Dirichlet", "Neumann"]);
u_bot_value_box = create_textbox(fig, 'Value:', '0', [30, 60, 120, 22]);

u_bot.Layout.Row = 7; 
u_bot.Layout.Column = 1;
u_bot_type_dropdown.Layout.Row = 8; 
u_bot_type_dropdown.Layout.Column = 1;
u_bot_value_box.Layout.Row = 8; 
u_bot_value_box.Layout.Column = 2;

%Top u
u_top = uilabel(fig, 'Text', 'u_top');
u_top_type_dropdown = uidropdown(fig, "Items", [ "FreeSlip", "Dirichlet", "Neumann"]);
u_top_value_box = create_textbox(fig, 'Value:', '0', [30, 60, 120, 22]);

u_top.Layout.Row = 9; 
u_top.Layout.Column = 1;
u_top_type_dropdown.Layout.Row = 10; 
u_top_type_dropdown.Layout.Column = 1;
u_top_value_box.Layout.Row = 10; 
u_top_value_box.Layout.Column = 2;

%Left v
v_lef = uilabel(fig, 'Text', 'v_left');
v_lef_type_dropdown = uidropdown(fig, "Items", [ "FreeSlip", "Dirichlet", "Neumann"]);
v_lef_value_box = create_textbox(fig, 'Value:', '0', [30, 60, 120, 22]);

v_lef.Layout.Row = 12;  
v_lef.Layout.Column = 1;
v_lef_type_dropdown.Layout.Row = 13;  
v_lef_type_dropdown.Layout.Column = 1;
v_lef_value_box.Layout.Row = 13;  
v_lef_value_box.Layout.Column = 2;

% Right V
v_rig = uilabel(fig, 'Text', 'v_right');
v_rig_type_dropdown = uidropdown(fig, "Items", [ "FreeSlip", "Dirichlet", "Neumann"]);
v_rig_value_box = create_textbox(fig, 'Value:', '0', [30, 60, 120, 22]);

v_rig.Layout.Row = 14;  
v_rig.Layout.Column = 1;
v_rig_type_dropdown.Layout.Row = 15;  
v_rig_type_dropdown.Layout.Column = 1;
v_rig_value_box.Layout.Row = 15;  
v_rig_value_box.Layout.Column = 2;

% Bottom V
v_bot = uilabel(fig, 'Text', 'v_bottom');
v_bot_type_dropdown = uidropdown(fig, "Items", ["Dirichlet", "FreeSlip", "Neumann"]);
v_bot_value_box = create_textbox(fig, 'Value:', 'cos(8*pi*x)', [30, 60, 120, 22]);

v_bot.Layout.Row = 16;  
v_bot.Layout.Column = 1;
v_bot_type_dropdown.Layout.Row = 17;  
v_bot_type_dropdown.Layout.Column = 1;
v_bot_value_box.Layout.Row = 17;  
v_bot_value_box.Layout.Column = 2;

% Top V
v_top = uilabel(fig, 'Text', 'v_top');
v_top_type_dropdown = uidropdown(fig, "Items", ["Neumann", "FreeSlip", "Dirichlet"]);
v_top_value_box = create_textbox(fig, 'Value:', '0', [30, 60, 120, 22]);

v_top.Layout.Row = 18;  
v_top.Layout.Column = 1;
v_top_type_dropdown.Layout.Row = 19;  
v_top_type_dropdown.Layout.Column = 1;
v_top_value_box.Layout.Row = 19;  
v_top_value_box.Layout.Column = 2;

panel_1 = uipanel(fig);
panel_1.Layout.Row = [3, 18];
panel_1.Layout.Column = 3;

panel_2 = uipanel(fig);
panel_2.Layout.Row = [3, 18];
panel_2.Layout.Column = 4;

panel_3 = uipanel(fig);
panel_3.Layout.Row = [3, 18];
panel_3.Layout.Column = 5;




%% Create the domain, set the geometry, and display it
grid_generation
if blocked ==1
    set_cell_type_blocked
else
    set_cell_type
end
plot_cell_type


%% Initial parameter values
BC.u_lef_type = "Dirichlet";        %x=0
BC.u_lef_value = "cos(8*pi*y)";
BC.u_rig_type = "Neumann";         %x=1
BC.u_rig_value = 0;
BC.u_bot_type = "FreeSlip";
BC.u_bot_value = 0;
BC.u_top_type = "FreeSlip";
BC.u_top_value = 0;

BC.v_lef_type = "FreeSlip";
BC.v_lef_value = 0;
BC.v_rig_type = "FreeSlip";
BC.v_rig_value = 0;
BC.v_bot_type = "Dirichlet";        %y=0
BC.v_bot_value = "cos(8*pi*x)";
BC.v_top_type = "Neumann";          %y=1
BC.v_top_value = 0;


% Create the submit button
submit_btn = uibutton(fig, 'Text', 'Initialize and Calculate with these Settings', 'Position', [10, 220, 120, 30], 'ButtonPushedFcn', @(btn, event) update_plot(fig, grids, N,base));
rainshadow_button = uibutton(fig, 'Text', 'Run Preset Rainshadow Simulation Instead', 'Position', [60, 220, 170, 30], 'ButtonPushedFcn', @(btn, event) rainshadow_callback(N, base));

% Store the UI component handles in guidata
guidata(fig, struct('u_lef_type_dropdown', u_lef_type_dropdown, 'u_lef_value_box', u_lef_value_box, ...
                    'u_rig_type_dropdown', u_rig_type_dropdown, 'u_rig_value_box', u_rig_value_box, ...
                    'u_bot_type_dropdown', u_bot_type_dropdown, 'u_bot_value_box', u_bot_value_box, ...
                    'u_top_type_dropdown', u_top_type_dropdown, 'u_top_value_box', u_top_value_box, ...
                    'v_lef_type_dropdown', v_lef_type_dropdown, 'v_lef_value_box', v_lef_value_box, ...
                    'v_rig_type_dropdown', v_rig_type_dropdown, 'v_rig_value_box', v_rig_value_box, ...
                    'v_bot_type_dropdown', v_bot_type_dropdown, 'v_bot_value_box', v_bot_value_box, ...
                    'v_top_type_dropdown', v_top_type_dropdown, 'v_top_value_box', v_top_value_box));

uiwait(base);

function rainshadow_callback(N, base)
    [rainshadow, blocked, BC, phi_in, noslipground] = begin_rainshadow(N, base);
    assignin('base', 'rainshadow', rainshadow);
    assignin('base', 'blocked', blocked);
    assignin('base', 'BC', BC);
    assignin('base', 'phi_in', phi_in);
    assignin('base', 'noslipground', noslipground);
end

function [rainshadow,blocked, BC, phi_in, noslipground] = begin_rainshadow(N, base)
    rainshadow = 1;
    blocked = 0;
    noslipground = 0;
            BC.u_lef_type = "Dirichlet";            % angled incoming wind
            BC.u_lef_value = 10* ones(N.y_u, 1);
            BC.u_rig_type = "FreeSlip";    % outlet
            BC.u_rig_value = zeros(N.y_u,1);
            BC.u_bot_type = "Dirichlet";            % no slip at ground
            BC.u_bot_value = zeros(N.x_u,1);
            BC.u_top_type = "FreeSlip";    % edge of atmosphere
            BC.u_top_value = zeros(N.x_u,1);

            BC.v_lef_type = "Dirichlet";            % horizontal wind
            BC.v_lef_value = zeros(N.y_v,1);
            BC.v_rig_type = "FreeSlip";    % outlet
            BC.v_rig_value = zeros(N.y_v,1);
            BC.v_bot_type = "Dirichlet";            % no slip at ground
            BC.v_bot_value = zeros(N.x_v,1);  
            BC.v_top_type = "Dirichlet";    % edge of atmosphere
            BC.v_top_value = zeros(N.x_v,1);
    x = linspace(-10, 10, N.y_p+1);      % Adjust range for sharpness
    phi_in = .25*(1- (1 ./ (1 + exp(-10-x))));      % Sigmoid function, adjusted to have a lot at the bottom
    % plot(phi_in,x)                                                                                              %%plotting the wind moisture profile
    % title('Inlet Humidity Profile')
    % ylabel('Altitude (Representative Values)')
    % xlabel('Humidity, aka Moisture Content (\Phi)')

    uiresume(base);
end

%% The callback function
function BC = update_plot(fig, grids, N, base)
    try
        % Retrieve the stored value from guidata
        data = guidata(fig);
        disp(fieldnames(data));  % This will help in debugging if the fields exist

        % Retrieve and evaluate inputs for boundary conditions
        BC.u_lef_value = evaluate_input(data.u_lef_value_box.Value, grids, N, 'u_left');
        BC.u_rig_value = evaluate_input(data.u_rig_value_box.Value, grids, N, 'u_right');
        BC.u_bot_value = evaluate_input(data.u_bot_value_box.Value, grids, N, 'u_bottom');
        BC.u_top_value = evaluate_input(data.u_top_value_box.Value, grids, N, 'u_top');
        
        % Force all to be column vectors
        if size(BC.u_lef_value, 1) == 1, BC.u_lef_value = BC.u_lef_value'; end
        if size(BC.u_rig_value, 1) == 1, BC.u_rig_value = BC.u_rig_value'; end
        if size(BC.u_bot_value, 1) == 1, BC.u_bot_value = BC.u_bot_value'; end
        if size(BC.u_top_value, 1) == 1, BC.u_top_value = BC.u_top_value'; end
        
        BC.u_lef_type  = data.u_lef_type_dropdown.Value;
        BC.u_rig_type  = data.u_rig_type_dropdown.Value;
        BC.u_bot_type  = data.u_bot_type_dropdown.Value;
        BC.u_top_type  = data.u_top_type_dropdown.Value;

        
        BC.v_lef_value = evaluate_input(data.v_lef_value_box.Value, grids, N, 'v_left');
        BC.v_rig_value = evaluate_input(data.v_rig_value_box.Value, grids, N, 'v_right');
        BC.v_bot_value = evaluate_input(data.v_bot_value_box.Value, grids, N, 'v_bottom');
        BC.v_top_value = evaluate_input(data.v_top_value_box.Value, grids, N, 'v_top');
        
        % Force all to be column vectors
        if size(BC.v_lef_value, 1) == 1, BC.v_lef_value = BC.v_lef_value'; end
        if size(BC.v_rig_value, 1) == 1, BC.v_rig_value = BC.v_rig_value'; end
        if size(BC.v_bot_value, 1) == 1, BC.v_bot_value = BC.v_bot_value'; end
        if size(BC.v_top_value, 1) == 1, BC.v_top_value = BC.v_top_value'; end
        
        BC.v_lef_type  = data.v_lef_type_dropdown.Value;
        BC.v_rig_type  = data.v_rig_type_dropdown.Value;
        BC.v_bot_type  = data.v_bot_type_dropdown.Value;
        BC.v_top_type  = data.v_top_type_dropdown.Value;

        assignin('base', 'BC', BC);
        disp('✅ BC stored in base workspace');
        uiresume(base);

    catch ME
        disp('❌ Error in update_plot:');
        disp(getReport(ME, 'extended'));
    end
end

function value = evaluate_input(input_str, grids, N, boundary_location)
    try
        % Replace 'x' and 'y' with grid variables
        input_str = strrep(input_str, 'x', 'grids.x_p');
        input_str = strrep(input_str, 'y', 'grids.y_p');

        func = str2func(['@(grids) ' input_str]);
        value = func(grids);

        if isscalar(value)
            % Determine which grid to use based on the boundary_location first part (u_ or v_)
            if contains(boundary_location, 'u_') || boundary_location(1) == 'u'
                % For u-velocity boundaries
                if contains(boundary_location, 'bottom') || contains(boundary_location, 'top')
                    value = value * ones(N.x_u, 1); % For u top/bottom: N.x_u long
                else % 'left' or 'right'
                    value = value * ones(N.y_u, 1); % For u sides: N.y_u long
                end
            elseif contains(boundary_location, 'v_') || boundary_location(1) == 'v'
                % For v-velocity boundaries
                if contains(boundary_location, 'bottom') || contains(boundary_location, 'top')
                    value = value * ones(N.x_v, 1); % For v top/bottom: N.x_v long
                else % 'left' or 'right'
                    value = value * ones(N.y_v, 1); % For v sides: N.y_v long
                end
            else
                % Default to pressure grid if boundary type is not specified
                if contains(boundary_location, 'bottom') || contains(boundary_location, 'top')
                    value = value * ones(N.x_p, 1); % Default to pressure grid
                else % 'left' or 'right'
                    value = value * ones(N.y_p, 1); % Default to pressure grid
                end
            end
            
            % Ensure column vector
            if size(value, 1) == 1
                value = value';
            end
        end
    catch ME
        error('Error evaluating input "%s": %s', input_str, ME.message);
    end
end



function textbox = create_textbox(fig, name, initial_value, position)
    % Create a text box for the UI
    textbox = uieditfield(fig, 'text', 'Position', position);
    textbox.Value = initial_value;
    textbox.Tooltip = name;
end
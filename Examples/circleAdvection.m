% Level set method

clear;

% Mesh parameters

% Number of grid points
nx = 205; % SET NUMBER OF GRID POINTS
ny = nx;

% Extents of domain
x1 = -1; % SET EXTENTS OF DOMAIN
x2 = 1;
y1= -1;
y2= 1;

% SDF dt
SDF_dt = 0.005; % SET REINITIALIZATION TIME STEP
% LDM dt
dt = 0.005; % SET LEVEL SET METHOD TIME STEP
% Total number of LSM time steps
num_time_steps = 400; % SET NUMBER OF LEVEL SET METHOD TIMESTEPS
% Time steps interval to save
save_step_interval = 100; % SET INTERVAL FOR TIMESTEPS TO SAVE

% Initialize grid
[X, Y, dx, dy] = create_grid(nx, ny, x1, x2, y1, y2);

% DEFINE LEVEL SET FUNCTIONS TO DESCRIBE INTERFACE FROM BASIC CIRCLES,
% RECTANGLES, AND PLANES
% Initialize level set field
center = [0, 0];
radius = 0.5;
phi = create_circle(X, Y, center, radius); 

%phi2 = create_rectangle(X, Y, [-0.1, -0.5], [0.1, 0.3]);

% COMBINE THE ABOVE LEVEL SET FUNCTIONS USING BOOLEAN OPERATIONS
% Combine several initial level set fields using boolean operations
%phi = combine_phi(phi, phi2, 'difference');

% Volume of interior region
Vol_i = calc_volume(phi, dx, dy);

% Convert phi to a signed distance function
phi = convert_to_SDF(phi, dx, dy, SDF_dt);

% Create external velocity field
[U, V] = create_velocity_field(X, Y);

% DEFINE A VELOCITY FIELD
function [U, V] = create_velocity_field(X, Y)
    U = 0.2+0*X;
    V = 0.2+0*Y;
end

% SET RATE OF MOTION IN THE NORMAL DIRECTION
a = 0;

% SET CURVATION COEFFICIENT
b = 0;

save_plot_initial(phi, X, Y);

% Increment phi over all time steps
% Create figure
fig = figure('Units', 'Inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');
for t = 1:num_time_steps
    % Increment phi
    phi = increment_phi(phi, U, V, a, b, dx, dy, dt);   
    % Convert to SDF
    phi = convert_to_SDF(phi, dx, dy, SDF_dt);
    % Compute volume
    vol_hist(t) = calc_volume(phi, dx, dy); 
    % Plot phi
    contourf(X, Y, phi, [0 0]);
    title_string = ['t=' num2str(t*dt)];
    tt = title(title_string, 'Interpreter','latex');
    add_plot_elements();
    drawnow;
    if mod(t, save_step_interval)==0
        save_plot(fig, t);
    end
end
%%
plot_vol_hist(vol_hist, Vol_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions

function add_plot_elements()
% Axes labels
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    
    % Plot borders
    box on;
    set(gca, 'LineWidth', 2);
    
    % Set font size for axes
    commonFigureFontSize = 15;
    ax = gca; % Get current axes
    ax.FontSize = commonFigureFontSize;

    axis square
end

function save_plot(fig, t)
    % Create the figure name
    figName = ['fig' num2str(t) '.pdf'];
    
    % Save the figure as a high-resolution PNG
    exportgraphics(fig, figName);
end

function save_plot_initial(phi, X, Y)
fig = figure('Units', 'Inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');
    % Plot phi
    contourf(X, Y, phi, [0 0]);
    title_string = ['t=' num2str(0)];
    tt = title(title_string, 'Interpreter','latex');
    add_plot_elements();
    save_plot(fig, 0);
end

function plot_vol_hist(vol_hist, Vol_i)
    fig = figure('Units', 'Inches', 'Position', [1, 1, 6, 4], 'PaperPositionMode', 'auto');
    p = plot(vol_hist/Vol_i);
    p.LineWidth=2;
    xlabel('Time step, $n$', 'Interpreter', 'latex');
    ylabel('$V_n/V_i$', 'Interpreter', 'latex');
    
    % Plot borders  
    box on;
    set(gca, 'LineWidth', 2);
    
    % Set font size for axes
    commonFigureFontSize = 15;
    ax = gca; % Get current axes
    ax.FontSize = commonFigureFontSize;
    save_plot(fig, 1);
end

% Function to initialize the grid
function [X, Y, dx, dy] = create_grid(nx, ny, x1, x2, y1, y2)
    x = linspace(x1, x2, nx);
    y = linspace(y1, y2, ny);

    dx = x(2)-x(1);
    dy = y(2)-y(1);

    [X, Y] = ndgrid(x, y);
end

% Function to initialize a circle
function phi = create_circle(X, Y, center, radius)
    phi = sqrt((X-center(1)).^2 + (Y-center(2)).^2) - radius;
end

% Function to initialize a rectangle
function phi = create_rectangle(X, Y, a, b)
    phi1 = create_plane(X, Y, a, [-1, 0]);
    phi2 = create_plane(X, Y, a, [0, -1]);
    phi3 = create_plane(X, Y, b, [1, 0]);
    phi4 = create_plane(X, Y, b, [0, 1]);
    phi5 = combine_phi(phi1, phi2, 'intersection');
    phi6 = combine_phi(phi3, phi4, 'intersection');
    phi = combine_phi(phi5, phi6, 'intersection');
end

% Function to initialize a plane
function phi = create_plane(X, Y, a, b)
    phi = (X-a(1)).*b(1) + (Y-a(2)).*b(2);
end

% Function to visualize phi as a contour
function plot_phi(X, Y, phi)
    figure;
    p = contourf(X, Y, phi, [0 0]);
end

% Function to combine two level set functions
function phi = combine_phi(phi1, phi2, operation)
    if strcmp(operation, 'union')
        phi = min(phi1, phi2);
    elseif strcmp(operation, 'intersection')
        phi = max(phi1, phi2);
    elseif strcmp(operation, 'difference')
        phi = max(phi1, -phi2);
    elseif strcmp(operation, 'complement')
        phi = -phi1;
    else
        disp('Invalid boolean operation!');
    end
end

% Obtain the forward differnce first derivative of phi in x
function phi_xp = Dxp(phi, dx)
    phi_xp = 0*phi;
    n = size(phi);
    rdx = 1/dx;
    for i=1:n(1)-1
        for j=1:n(2)
            phi_xp(i, j) = rdx*(phi(i+1, j) - phi(i, j));
        end
    end
end

% Obtain the backward differnce first derivative of phi in x
function phi_xm = Dxm(phi, dx)
    phi_xm = 0*phi;
    n = size(phi);
    rdx = 1/dx;
    for i=2:n(1)
        for j=1:n(2)
            phi_xm(i, j) = rdx*(phi(i, j) - phi(i-1, j));
        end
    end
end

% Obtain the central differnce first derivative of phi in x
function phi_xc = Dxc(phi, dx)
    phi_xc = 0*phi;
    n = size(phi);
    rdx = 0.5/dx;
    for i=2:n(1)-1
        for j=1:n(2)
            phi_xc(i, j) = rdx*(phi(i+1, j) - phi(i-1, j));
        end
    end
end

% Obtain the central difference second derivative of phi in x
function phi_xx = Dxx(phi, dx)
    phi_xx = 0*phi;
    n = size(phi);
    rdx = 1/dx/dx;
    for i=2:n(1)-1
        for j=1:n(2)
            phi_xx(i, j) = rdx*(phi(i+1, j) - 2*phi(i, j) + phi(i-1, j));
        end
    end
end

% Obtain the central difference first derivative of phi in y
function phi_yc = Dyc(phi, dy)
    phi_yc = 0*phi;
    n = size(phi);
    rdy = 0.5/dy;
    for i=1:n(1)
        for j=2:n(2)-1
            phi_yc(i, j) = rdy*(phi(i, j+1) - phi(i, j-1));
        end
    end
end

% Obtain the central difference second derivative of phi in y
function phi_yy = Dyy(phi, dy)
    phi_yy = 0*phi;
    n = size(phi);
    rdy = 1/dy/dy;
    for i=1:n(1)
        for j=2:n(2)-1
            phi_yy(i, j) = rdy*(phi(i, j+1) -2*phi(i, j) + phi(i, j-1));
        end
    end
end

% Obtain the central difference second derivative of phi in x and y
function phi_xy = Dxy(phi, dx, dy)
    phi_xy = 0*phi;
    n = size(phi);
    rdxdy = 0.25/dx/dy;
    for i=2:n(1)-1
        for j=2:n(2)-1
            phi_xy(i, j) = rdxdy*(phi(i+1, j+1) + phi(i-1, j-1) - phi(i+1, j-1) - phi(i-1, j+1));
        end
    end
end

% Obtain the forward differnce first derivative of phi in y
function phi_yp = Dyp(phi, dy)
    phi_yp = 0*phi;
    n = size(phi);
    rdy = 1/dy;
    for i=1:n(1)
        for j=1:n(2)-1
            phi_yp(i, j) = rdy*(phi(i, j+1) - phi(i, j));
        end
    end
end

% Obtain the backward differnce first derivative of phi in y
function phi_ym = Dym(phi, dy)
    phi_ym = 0*phi;
    n = size(phi);
    rdy = 1/dy;
    for i=1:n(1)
        for j=2:n(2)
            phi_ym(i, j) = rdy*(phi(i, j) - phi(i, j-1));
        end
    end
end

% Function to convert phi to a SDF
function phi = convert_to_SDF(phi, dx, dy, SDF_dt)

    % Numerically smeared sign function of phi
    S = phi./sqrt(phi.^2 + dx.^2);
   
    % Time step
    dt = SDF_dt;

    % Iterations to solve the reinitialization equation
    for i = 1:1000

        % Normal gradient of phi
        norm_phi = normal_gradient_phi(phi, S, dx, dy);

        % Update phi
        phi = phi - dt * S .* (norm_phi - 1);

        % Compute error
        err = norm(normal_gradient_phi(phi, S, dx, dy) - 1) / numel(phi);
        if (err < 1e-1)
            return
        end
    end
end

% Calculate the normal gradient of phi using its derivatives and the
% smeared sign function
function norm_phi = normal_gradient_phi(phi, S, dx, dy)

    % Forward differnce derivative of phi in x
    phi_xp = Dxp(phi, dx);
    % Backward differnce derivative of phi in x
    phi_xm = Dxm(phi, dx);
    % Forward differnce derivative of phi in y
    phi_yp = Dyp(phi, dy);
    % Backward differnce derivative of phi in y
    phi_ym = Dym(phi, dy);

    phi_x2 = 0*S;
    phi_y2 = 0*S;
    n = size(S);
    for i=1:n(1)
        for j=1:n(2)
            if(S(i, j)>=0)
                phi_x2(i, j) = max(max(phi_xm(i, j), 0)^2, min(phi_xp(i, j), 0)^2);
                phi_y2(i, j) = max(max(phi_ym(i, j), 0)^2, min(phi_yp(i, j), 0)^2);
            else
                phi_x2(i, j) = max(min(phi_xm(i, j), 0)^2, max(phi_xp(i, j), 0)^2);
                phi_y2(i, j) = max(min(phi_ym(i, j), 0)^2, max(phi_yp(i, j), 0)^2);
            end
        end
    end
    norm_phi = sqrt(phi_x2 + phi_y2);
end

% Increment phi
function phi = increment_phi(phi, U, V, a, b, dx, dy, dt)
    n = size(phi);
    % Forward difference derivative of phi in x
    phi_xp = Dxp(phi, dx);
    % Backward difference derivative of phi in x
    phi_xm = Dxm(phi, dx);
    % Forward difference derivative of phi in y
    phi_yp = Dyp(phi, dy);
    % Backward difference derivative of phi in y
    phi_ym = Dym(phi, dy);

    % Derivatives for motion under curvature
    phi_xc = Dxc(phi, dx);
    phi_yc = Dyc(phi, dy);

    phi_xx = Dxx(phi, dx);
    phi_yy = Dyy(phi, dy);
    phi_xy = Dxy(phi, dx, dy);

    px = 0*phi;
    py = 0*phi;
    H_v = 0*phi;

    for i=1:n(1)
        for j=1:n(2)
            if U(i,j)>=0
                px(i, j) = phi_xm(i, j);
            else
                px(i, j) = phi_xp(i, j); 
            end
            
            if V(i,j)>=0
                py(i, j) = phi_ym(i, j);
            else
                py(i, j) = phi_yp(i, j); 
            end
            
            H_v(i, j) = U(i, j) * px(i, j) + V(i, j) * py(i, j);
        end
    end

    % Compute the Hamiltonian for motion in the normal direction.
    H_normal = a .* normal_gradient_phi(phi, a, dx, dy);
    
    % Curvature calculation
    numr = phi_xc.^2 .* phi_yy - 2 * phi_xc .* phi_yc .* phi_xy + phi_yc.^2 .* phi_xx;
    denr = phi_xc.^2 + phi_yc.^2;
    denr(find(denr==0))=1;
    kappa_grad_phi = numr./denr;
    H_k = - b*kappa_grad_phi;

    kappa = numr./(denr.^(3/2));

    phi = phi - dt * H_v;
    phi = phi - dt * H_normal;
    phi = phi - dt * H_k;
end

% Function ot calculate the volume of phi
function vol = calc_volume(phi, dx, dy)
    n = size(phi);
    vol = 0;
    eps = 1e-3;
    for i=1:n(1)
        for j=1:n(2)
            vol = vol + (1-Hs(phi(i, j), eps))*dx*dy;
        end
    end
end

% Smeared heaviside function
function H = Hs(x, eps)
    if x<-eps
        H = 0;
    elseif x>eps
        H = 1;
    else
        H = 0.5*(1+x/eps + 1/pi * sin(pi*x/eps));
    end
end

% Function ot calculate the surface area of phi
function sa = calc_surface(phi, dx, dy)
    n = size(phi);
    sa = 0;
    eps = 1e-1;
    for i=1:n(1)
        for j=1:n(2)
            sa = sa + Dm(phi(i, j), eps)*dx*dy;
        end
    end
end

% Mollified Delta function
function D = Dm(x, eps)
    if abs(x) < eps
        D = (0.5/eps)*(1 + cos(pi*x/eps));
    else
        D = 0;
    end
end


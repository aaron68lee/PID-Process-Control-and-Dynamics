clc; 
clear all; 
close all;

%% Problem 2
% Utilize fsolve to solve for the steady state values, 3 different cases of
x_steady_states = zeros(4, 3);
steadystate_function = @calculate_steadystate;

for i = 1:3
    n = i - 2;
    initial_guess = [n n n n];
    x_steady_states(:, i) = fsolve(steadystate_function, initial_guess);
end

%% Problem 3
% Set constants
Da1 = 1e6; 
Da2 = 1e7; 
Da1p = 1.5 * 1e6; 
Da2p = 0; 
x30 = 0.025;
x40 = 0.025; 
E1 = 1.0066; 
E2 = 1.0532; 
U = 8; 
e1 = 12.5; 
e2 = 40; 
e3 = 1;

% Calculate Jacobian matrix for each steady state
for i = 1:3
    x = x_steady_states(:, i);
    J = calculate_jacobian(x, Da1, Da2, Da1p, Da2p, E1, E2, U, e1, e2, e3);
    eigenvalues(:, i) = eig(J);

    if any(eigenvalues(:, i) > 0)
        unstable_eigenvalues = eigenvalues(:, i);
        x_unstable = x;
    end
end

%% Problem 5
% Calculate Jacobian for unstable steady state
J_unstable = calculate_jacobian(x_unstable, Da1, Da2, Da1p, Da2p, E1, E2, U, e1, e2, e3);

%% Problem 6
% Plot step response
Gp = zpk(-1.2696, [64.3083, -0.972544, -601.204], 4000);
Kc = 350;
taui = 0.1;
Gc = tf([Kc * taui, Kc], [taui, 0]);
y_pi = (Gp * Gc) / (1 + Gp * Gc);
step(y_pi);

%% Problem 7
% Set constants
d = 0.001;
y_sp = 0.005;
x0 = zeros(5, 1);
xs = [x_unstable' 0];
tset = [0 0.1];
tdist = [0 0.3];

% Linear 0.005
[t, x_linear_005] = ode45(@(t, x) linear_system(t, x, J_unstable, d, y_sp, Kc, taui), tset, x0);
e_linear_005 = y_sp - x_linear_005(:, 3);
u_linear_005 = Kc * e_linear_005 + Kc / taui * x_linear_005(:, 5);

% Nonlinear 0.005
[t, x_nonlinear_005] = ode45(@(t, x) nonlinear_system(t, x, d, y_sp, Kc, taui, xs), tset, x0);
e_nonlinear_005 = y_sp - x_nonlinear_005(:, 3);
u_nonlinear_005 = Kc * e_nonlinear_005 + Kc / taui * x_nonlinear_005(:, 5);

% Plotting
plot_nonlinear_system_response(t, u_linear_005, x_linear_005(:, 3), '(magnitude 0.005)');
plot_nonlinear_system_response(t, u_nonlinear_005, x_nonlinear_005(:, 3), '(magnitude 0.005)');

% Set the set point to 0.01
y_sp = 0.01;

% Linear 0.01
[t, x_linear_01] = ode45(@(t, x) linear_system(t, x, J_unstable, d, y_sp, Kc, taui), tset, x0);
e_linear_01 = y_sp - x_linear_01(:, 3);
u_linear_01 = Kc * e_linear_01 + Kc / taui * x_linear_01(:, 5);

% Nonlinear 0.01
[t, x_nonlinear_01] = ode45(@(t, x) nonlinear_system(t, x, d, y_sp, Kc, taui, xs), tset, x0);
e_nonlinear_01 = y_sp - x_nonlinear_01(:, 3);
u_nonlinear_01 = Kc * e_nonlinear_01 + Kc / taui * x_nonlinear_01(:, 5);

% Plotting
plot_nonlinear_system_response(t, u_linear_01, x_linear_01(:, 3), '(magnitude 0.01)');
plot_nonlinear_system_response(t, u_nonlinear_01, x_nonlinear_01(:, 3), '(magnitude 0.01)');

%% Problem 8
d = 0.001; % Set the disturbance to 0.001
y_sp = 0;

% Plot linear, d = 0.001
[t, x_linear_dist_001] = ode45(@(t, x) linear_system(t, x, J_unstable, d, y_sp, Kc, taui), tdist, x0);
e_linear_dist_001 = y_sp - x_linear_dist_001(:, 3);
u_linear_dist_001 = Kc * e_linear_dist_001 + Kc / taui * x_linear_dist_001(:, 5);

% Plot nonlinear, d = 0.001
[t, x_nonlinear_dist_001] = ode45(@(t, x) nonlinear_system(t, x, d, y_sp, Kc, taui, xs), tdist, x0);
e_nonlinear_dist_001 = y_sp - x_nonlinear_dist_001(:, 3);
u_nonlinear_dist_001 = Kc * e_nonlinear_dist_001 + Kc / taui * x_nonlinear_dist_001(:, 5);

% Plotting
plot_nonlinear_system_response(t, u_linear_dist_001, x_linear_dist_001(:, 3), '(Disturbance 0.001)');
plot_nonlinear_system_response(t, u_nonlinear_dist_001, x_nonlinear_dist_001(:, 3), '(Disturbance 0.001)');

d = 0.002; % Set the disturbance to 0.002

% Plot linear, d = 0.002
[t, x_linear_dist_002] = ode45(@(t, x) linear_system(t, x, J_unstable, d, y_sp, Kc, taui), tdist, x0);
e_linear_dist_002 = y_sp - x_linear_dist_002(:, 3);
u_linear_dist_002 = Kc * e_linear_dist_002 + Kc / taui * x_linear_dist_002(:, 5);

% Plot nonlinear, d = 0.002
[t, x_nonlinear_dist_002] = ode45(@(t, x) nonlinear_system(t, x, d, y_sp, Kc, taui, xs), tdist, x0);
e_nonlinear_dist_002 = y_sp - x_nonlinear_dist_002(:, 3);
u_nonlinear_dist_002 = Kc * e_nonlinear_dist_002 + Kc / taui * x_nonlinear_dist_002(:, 5);

% Plotting
plot_nonlinear_system_response(t, u_linear_dist_002, x_linear_dist_002(:, 3), '(Disturbance 0.002)');
plot_nonlinear_system_response(t, u_nonlinear_dist_002, x_nonlinear_dist_002(:, 3), '(Disturbance 0.002)');

%% Problem 9
y_sp = 0.005; % Set point change to 0.005
d = 0;

% Nonlinear set point change
[t, x_nonlinear_setpoint] = ode45(@(t, x) nonlinear_system_new(t, x, d, y_sp, Kc, taui, xs), tset, x0);
e_nonlinear_setpoint = y_sp - x_nonlinear_setpoint(:, 3);
u_nonlinear_setpoint = Kc * e_nonlinear_setpoint + Kc / taui * x_nonlinear_setpoint(:, 5);

% Plotting
plot_nonlinear_system_response(t, u_nonlinear_setpoint, x_nonlinear_setpoint(:, 3), '(Set Point Change 10^5)');

y_sp = 0;
d = 0.002; % Set the disturbance to 0.002

% Nonlinear disturbance
[t, x_nonlinear_disturbance] = ode45(@(t, x) nonlinear_system_new(t, x, d, y_sp, Kc, taui, xs), tdist, x0);
e_nonlinear_disturbance = y_sp - x_nonlinear_disturbance(:, 3);
u_nonlinear_disturbance = Kc * e_nonlinear_disturbance + Kc / taui * x_nonlinear_disturbance(:, 5);

% Plotting
plot_nonlinear_system_response(t, u_nonlinear_disturbance, x_nonlinear_disturbance(:, 3), '(Disturbance 10^5)');

%% Functions %%

function plot_nonlinear_system_response(time, input, output, title_text)
    % Plot nonlinear system response
    
    figure();
    
    subplot(2, 1, 1);
    plot(time, input);
    title(['Nonlinear Input Response ' title_text]);
    xlabel('Time');
    ylabel('Manipulated Input, u');
    
    subplot(2, 1, 2);
    plot(time, output);
    title(['Nonlinear Output Response ' title_text]);
    xlabel('Time');
    ylabel('Controlled Output, y');
end

function x_ss = calculate_steadystate()
    % Calculate steady state values with different initial guesses
    x_ss = zeros(4, 3);
    fun = @steadystate;
    for i = 1:3
        n = i - 2;
        x0 = [n n n n];
        x_ss(:, i) = fsolve(fun, x0);
    end
end

function J = calculate_jacobian(x_ss)
    % Calculate Jacobian matrix for each steady state
    J = zeros(4, 4, 3);
    for i = 1:3
        x = x_ss(:, i);
        J(:, :, i) = calculate_jacobian_matrix(x);
    end
end

function J = calculate_jacobian_matrix(x)

    % Calculate Jacobian matrix for a given steady state
    J = zeros(4, 4);

    Da1 = 1e6; Da2 = 1e7; Da1p = 1.5 * 1e6; Da2p = 0;
    E1 = 1.0066; E2 = 1.0532; U = 8; e1 = 12.5; e2 = 40; e3 = 1;

    % Fill in the Jacobian matrix elements
    J(1, 1) = -1 - Da1 * exp(-E1 / x(3));
    J(1, 2) = 0;
    J(1, 3) = -Da1 * E1 / x(3)^2 * exp(-E1 / x(3)) * x(1);
    J(1, 4) = 0;

    J(2, 1) = Da1 * exp(-E1 / x(3));
    J(2, 2) = -1 - Da2 * exp(-E2 / x(3));
    J(2, 3) = Da1 * E1 / x(3)^2 * exp(-E1 / x(3)) * x(1) - Da2 * E2 / x(3)^2 * exp(-E2 / x(3)) * x(2);
    J(2, 4) = 0;

    J(3, 1) = Da1p * exp(-E1 / x(3));
    J(3, 2) = Da2p * exp(-E2 / x(3));
    J(3, 3) = -1 + Da1p * E1 / x(3)^2 * exp(-E1 / x(3)) * x(1) - Da2p * E2 / x(3)^2 * exp(-E2 / x(3)) * x(2) - U;
    J(3, 4) = U;

    J(4, 1) = 0;
    J(4, 2) = 0;
    J(4, 3) = U * e1 * e3;
    J(4, 4) = -e1 * e2 - U * e1 * e3;
end

function dxdt = linear_system(t, x, J, d, y_sp, Kc, taui)
    % Linear system dynamics
    dxdt = zeros(5, 1);
    dxdt(1) = J(1, 1)*x(1) + J(1, 3)*x(3);
    dxdt(2) = J(2, 1)*x(1) + J(2, 2)*x(2) + J(2, 3)*x(3);
    dxdt(3) = J(3, 1)*x(1) + J(3, 3)*x(3) + J(3, 4)*x(4) + d;
    dxdt(4) = J(4, 3)*x(3) + J(4, 4)*x(4) + 12.5 * 40 * (Kc * (y_sp - x(3)) + Kc/taui * x(5));
    dxdt(5) = y_sp - x(3);
end

function dxdt = nonlinear_system(t, x, d, y_sp, Kc, taui, xs)
    % Nonlinear system dynamics
    Da1 = 1e6; Da2 = 1e7; Da1p = 1.5 * 1e6; Da2p = 0; x30 = 0.025; x40 = 0.025;
    E1 = 1.0066; E2 = 1.0532; U = 8; e1 = 12.5; e2 = 40; e3 = 1;

    dxdt = zeros(5, 1);

    dxdt(1) = -x(1) - Da1 * exp(-E1 / (x(3) + xs(3))) * (x(1) + xs(1)) + Da1 * exp(-E1 / xs(3)) * xs(1);
    dxdt(2) = -x(2) + Da1 * exp(-E1 / (x(3) + xs(3))) * (x(1) + xs(1)) - Da1 * exp(-E1 / xs(3)) * xs(1) ...
            - Da2 * exp(-E2 / (x(3) + xs(3))) * (x(2) + xs(2)) + Da2 * exp(-E2 / xs(3)) * xs(2);
    dxdt(3) = d - x(3) + Da1p * exp(-E1 / (x(3) + xs(3))) * (x(1) + xs(1)) - Da1p * exp(-E1 / xs(3)) * xs(1) ...
            + Da2p * exp(-E2 / (x(3) + xs(3))) * (x(2) + xs(2)) - Da2p * exp(-E2 / xs(3)) * xs(2) - U * (x(3) - x(4));
    dxdt(4) = e1 * e2 * (Kc * (y_sp - x(3)) + Kc / taui * x(5) - x(4)) + U * e1 * e3 * (x(3) - x(4));
    dxdt(5) = y_sp - x(3);
end

function dxdt = nonlinear_system_new(t, x, d, y_sp, Kc, taui, xs)
    % New nonlinear system dynamics
    Da1 = 1e6; Da2 = 1e7; Da1p = 1.5 * 1e6; Da2p = 1e5; x30 = 0.025;

    dxdt = zeros(5, 1);

    dxdt(1) = -x(1) - Da1 * exp(-E1 / (x(3) + xs(3))) * (x(1) + xs(1)) + Da1 * exp(-E1 / xs(3)) * xs(1);
    dxdt(2) = -x(2) + Da1 * exp(-E1 / (x(3) + xs(3))) * (x(1) + xs(1)) - Da1 * exp(-E1 / xs(3)) * xs(1) ...
            - Da2 * exp(-E2 / (x(3) + xs(3))) * (x(2) + xs(2)) + Da2 * exp(-E2 / xs(3)) * xs(2);
    dxdt(3) = d - x(3) + Da1p * exp(-E1 / (x(3) + xs(3))) * (x(1) + xs(1)) - Da1p * exp(-E1 / xs(3)) * xs(1) ...
            + Da2p * exp(-E2 / (x(3) + xs(3))) * (x(2) + xs(2)) - Da2p * exp(-E2 / xs(3)) * xs(2) - U * (x(3) - x(4));
    dxdt(4) = e1 * e2 * (Kc * (y_sp - x(3)) + Kc / taui * x(5) - x(4)) + U * e1 * e3 * (x(3) - x(4));
    dxdt(5) = y_sp - x(3);
end


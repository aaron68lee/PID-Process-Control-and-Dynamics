clc; 
clear all; 
close all;

%% Problem 2
% Utilize fsolve to solve for the steady state values, 3 different cases of 
% initial guesses 
x_ss = zeros(4,3); 
fun = @steadystate; 

for i = 1:3 
    n = i - 2; 
    x0 = [n n n n]; 
    x_ss(:,i) = fsolve(fun, x0); 
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

% Jacobian matrix for each steady state 
for i = 1:3 
    x = x_ss(:,i); 
    J(1,1) = -1 - Da1 * exp(-E1 / x(3)); 
    J(1,2) = 0; 
    J(1,3) = -Da1 * E1 / x(3)^2 * exp(-E1 / x(3)) * x(1); 
    J(1,4) = 0; 
    J(2,1) = Da1 * exp(-E1 / x(3)); 
    J(2,2) = -1 - Da2 * exp(-E2 / x(3)); 
    J(2,3) = Da1 * E1 / x(3)^2 * exp(-E1 / x(3)) * x(1) - Da2 * E2 / x(3)^2 * exp(-E2 / x(3)) * x(2); 
    J(2,4) = 0; 
    J(3,1) = Da1p * exp(-E1 / x(3)); 
    J(3,2) = Da2p * exp(-E2 / x(3)); 
    J(3,3) = -1 + Da1p * E1 / x(3)^2 * exp(-E1 / x(3)) * x(1) - Da2p * E2 / x(3)^2 * exp(-E2 / x(3)) * x(2) - U; 
    J(3,4) = U; 
    J(4,1) = 0; 
    J(4,2) = 0; 
    J(4,3) = U * e1 * e3; 
    J(4,4) = -e1 * e2 - U * e1 * e3;

    % Calculate eigenvalues 
    eigen(:,i) = eig(J);

    % Grab the unsteady values, positive eigenvalues 
    if any(eigen(:,i) > 0) 
        unstable = eigen(:,i); 
        x_unstable = x; 
    end 
end

%% Problem 5
% Calculate Jacobian for unstable steady state 
J_unstable(1,1) = -1 - Da1 * exp(-E1 / x_unstable(3)); 
J_unstable(1,2) = 0; 
J_unstable(1,3) = -Da1 * E1 / x_unstable(3)^2 * exp(-E1 / x_unstable(3)) * x_unstable(1); 
J_unstable(1,4) = 0; 
J_unstable(2,1) = Da1 * exp(-E1 / x_unstable(3)); 
J_unstable(2,2) = -1 - Da2 * exp(-E2 / x_unstable(3)); 
J_unstable(2,3) = Da1 * E1 / x_unstable(3)^2 * exp(-E1 / x_unstable(3)) * x_unstable(1) - Da2 * E2 / x_unstable(3)^2 * exp(-E2 / x_unstable(3)) * x_unstable(2); 
J_unstable(2,4) = 0; 
J_unstable(3,1) = Da1p * exp(-E1 / x_unstable(3)); 
J_unstable(3,2) = Da2p * exp(-E2 / x_unstable(3));

J_unstable(3,3) = -1 + Da1p * E1 / x_unstable(3)^2 * exp(- E1 / x_unstable(3)) * x_unstable(1) - Da2p * E2 / x_unstable(3)^2 * exp(- E2 / x_unstable(3)) * x_unstable(2) - U; 
J_unstable(3,4) = U; 
J_unstable(4,1) = 0; 
J_unstable(4,2) = 0; 
J_unstable(4,3) = U * e1 * e3; 
J_unstable(4,4) = -e1 * e2 - U * e1 * e3; 

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
d = 0; % set the set point change to 0.005 
y_sp = 0.005; 
x0 = zeros(5,1); 
xs = [rot90(x_unstable) 0]; 
tset = [0 0.1]; 
tdist = [0 0.3];

% Linear 0.005 
% Solves the differential equations 
[t, x] = ode45(@(t,x)Linearization(x,J_unstable,d,y_sp,Kc,taui), tset, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Linear Input Response (magnitude 0.005)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Linear Output Response (magnitude 0.005)') 
xlabel('Time') 
ylabel('Controlled Output, y')

% Nonlinear 0.005
% Solves the differential equations 
[t, x] = ode45(@(t,x)NonLinearization(x,d,y_sp,Kc,taui,xs), tset, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Nonlinear Input Response (magnitude 0.005)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Nonlinear Output Response (magnitude 0.005)') 
xlabel('Time') 
ylabel('Controlled Output, y')

% Set the set point to 0.01 
y_sp = 0.01; 

% Linear 0.01
% Solves the differential equations 
[t, x] = ode45(@(t,x)Linearization(x,J_unstable,d,y_sp,Kc,taui), tset, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Linear Input Response (Magnitude 0.01) ') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Linear Output Response (Magnitude 0.01) ') 
xlabel('Time') 
ylabel('Controlled Output, y')

% Nonlinear 0.01
% Solves the differential equations 
[t, x] = ode45(@(t,x)NonLinearization(x,d,y_sp,Kc,taui,xs), tset, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Nonlinear Input Response (Magnitude 0.01)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Nonlinear Output Response (Magnitude 0.01)') 
xlabel('Time') 
ylabel('Controlled Output, y')

%% Problem 8 
d = 0.001; % set the disturbance to 0.001 
y_sp = 0;

% Plot linear, d = 0.001
% Solves the differential equations 
[t, x] = ode45(@(t,x)Linearization(x,J_unstable,d,y_sp,Kc,taui), tdist, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Linear Input Response (Disturbance 0.001)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Linear Output Response (Disturbance 0.001)') 
xlabel('Time') 
ylabel('Controlled Output, y')

% Plot nonlinear, d = 0.001
% Solves the differential equations 
[t, x] = ode45(@(t,x)NonLinear(x,d,y_sp,Kc,taui,xs), tdist, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Nonlinear Input Response (Disturbance 0.001)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Nonlinear Output Response (Disturbance 0.001)') 
xlabel('Time') 
ylabel('Controlled Output, y')

d = 0.002; % set the disturbance to 0.002

% Plot linear, d = 0.002
% Solves the differential equations 
[t, x] = ode45(@(t,x)Linear(x,J_unstable,d,y_sp,Kc,taui), tdist, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Linear Input Response (Disturbance 0.002)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Linear Output Response (Disturbance 0.002)') 
xlabel('Time') 
ylabel('Controlled Output, y')

% Plot nonlinear, d = 0.002 
[t, x] = ode45(@(t,x)NonLinear(x,d,y_sp,Kc,taui,xs), tdist, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Nonlinear Input Response (Disturbance 0.002)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Nonlinear Output Response (Disturbance 0.002)') 
xlabel('Time') 
ylabel('Controlled Output, y')

%% Problem 9 
y_sp = 0.005; % first case with the set point change to 0.005 
d = 0;

% Solves the differential equations 
[t, x] = ode45(@(t,x)NonLinearNew(x,d,y_sp,Kc,taui,xs), tset, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Nonlinear Set Point Change Input Response (10^5)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Nonlinear Set Point Change Output Response (10^5)') 
xlabel('Time') 
ylabel('Controlled Output, y')

y_sp = 0; 
d = 0.002; % set case with disturbance to 0.002 

% Solves the differential equations 
[t, x] = ode45(@(t,x)NonLinearNew(x,d,y_sp,Kc,taui,xs), tdist, x0); 
e = y_sp - x(:,3); 
u = Kc * e + Kc / taui * x(:,5);

% Plotting 
figure(); 
subplot(2,1,1); 
plot(t, u); 
title('Nonlinear Disturbance Input Response (10^5)') 
xlabel('Time') 
ylabel('Manipulated Input, u')

subplot(2,1,2); 
plot(t, x(:,3)); 
title('Nonlinear Disturbance Output Response (10^5)') 
xlabel('Time') 
ylabel('Controlled Output, y')

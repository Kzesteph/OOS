% Full MATLAB script for simulating and optimizing a Buck converter

% global variables
global Vin L C R_bat Vbat_nominal C_bat SOC_initial Vout_ref t_sim dt

% Parameters of the Buck Converter and Battery
Vin = 12;                % Input voltage (V)
L = 1e-3;              % Inductor (H)
C = 100e-6;             % Capacitor (F)
R_bat = 1;              % Equivalent battery resistance (Ohms)
Vbat_nominal = 1.5;     % Nominal battery voltage (V)
C_bat = 20e-3 * 3600;  % Battery capacity (Ampere-seconds)(extra low capacity just for the test)
SOC_initial = 0.9;      % Initial state of charge (SOC, between 0 and 1)
Vout_ref = 1.5;         % Desired output voltage (V)
t_sim = 2;              % Simulation time (s)
dt = 1e-3;              % Time step (s)

% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Initial guess for PID gains
initial_gains = [0.1, 100, 0];  % [Kp, Ki, Kd]

% Lower and upper bounds for the PID gains
lb = [0, 0, 0];  % Gains cannot be negative
ub = [1, 1000, 1];  % Reasonable upper bounds

% Optimization with constraints
[optimal_gains, optimal_cost] = fmincon(@objectiveFunction, initial_gains, [], [], [], [], lb, ub, @constraintsFunction, options);

% Display the optimal gains and the corresponding cost
disp('Optimal PID Gains:');
disp(['Kp = ', num2str(optimal_gains(1))]);
disp(['Ki = ', num2str(optimal_gains(2))]);
disp(['Kd = ', num2str(optimal_gains(3))]);
disp(['Optimal Cost (ISE) = ', num2str(optimal_cost)]);

% Simulation with Optimal Parameters
Kp_opt = optimal_gains(1);
Ki_opt = optimal_gains(2);
Kd_opt = optimal_gains(3);

% Initial Conditions
iL0 = 0;              
vC0 = 0;              
SOC0 = SOC_initial;   
x0 = [iL0; vC0; SOC0];  

% Time vector for simulation
time = 0:dt:t_sim;

% Variables initialization to store simulation results
iL = zeros(size(time));
vC = zeros(size(time));
SOC = zeros(size(time));
D = zeros(size(time));

% PID control variables
error_integral = 0;
prev_error = 0;

% Simulation loop with optimal PID gains
for k = 1:length(time)
    x = x0;

    % Duty cycle calculation
    [D(k), error_integral, prev_error] = PIDController(x(2), Vout_ref, Kp_opt, Ki_opt, Kd_opt, dt, error_integral, prev_error);

    % Computomputation of the next state including SOC using ODE solver
    [~, x_next] = ode45(@(t, x) buckConverterBatteryODE(t, x, Vin, L, C, x(3), R_bat, Vbat_nominal, C_bat, D(k)), [0 dt], x);

    % Store results
    iL(k) = x_next(end, 1);
    vC(k) = x_next(end, 2);
    SOC(k) = x_next(end, 3);

    % Update initial conditions for the next iteration
    x0 = x_next(end, :);
end

% Plot the results
figure;
subplot(4, 1, 1);
plot(time, iL);
title('Inductor Current (Optimal PID)');
xlabel('Time (s)');
ylabel('Current (A)');

subplot(4, 1, 2);
plot(time, vC);
title('Output Voltage (Optimal PID)');
xlabel('Time (s)');
ylabel('Voltage (V)');

subplot(4, 1, 3);
plot(time, SOC);
title('State of Charge (SOC) (Optimal PID)');
xlabel('Time (s)');
ylabel('SOC (0 to 1)');

subplot(4, 1, 4);
plot(time, D);
title('Duty Cycle (Optimal PID)');
xlabel('Time (s)');
ylabel('Duty Cycle');

% Objective Function to Minimize ISE
function cost = objectiveFunction(pid_gains)
    global Vin L C R_bat Vbat_nominal C_bat SOC_initial Vout_ref t_sim dt

    % Unpack PID gains
    Kp = pid_gains(1);
    Ki = pid_gains(2);
    Kd = pid_gains(3);
    
    % Initial Conditions
    iL0 = 0;
    vC0 = 0;
    SOC0 = SOC_initial;
    x0 = [iL0; vC0; SOC0];  

    % Time vector for simulation
    time = 0:dt:t_sim;

    % Error variable initialization
    total_error = 0;

    % PID control variables
    error_integral = 0;
    prev_error = 0;

    % Simulation loop
    for k = 1:length(time)
        x = x0;

        % Duty cycle Calculation
        [D, error_integral, prev_error] = PIDController(x(2), Vout_ref, Kp, Ki, Kd, dt, error_integral, prev_error);

        % Computation of the next state including SOC using ODE solver
        [~, x_next] = ode45(@(t, x) buckConverterBatteryODE(t, x, Vin, L, C, x(3), R_bat, Vbat_nominal, C_bat, D), [0 dt], x);

        % error Computation
        error = Vout_ref - x_next(end, 2);
        total_error = total_error + error^2;

        % Initial conditions for the next iteration
        x0 = x_next(end, :);
    end

    % The cost is the integral of the squared error (ISE)
    cost = total_error;
end

% PID Controller Function
function [D, error_integral, prev_error] = PIDController(vC, Vout_ref, Kp, Ki, Kd, dt, error_integral, prev_error)
    error = Vout_ref - vC;
    error_integral = error_integral + error * dt;
    error_derivative = (error - prev_error) / dt;
    
    % PID control
    D = Kp * error + Ki * error_integral + Kd * error_derivative;

    % Duty cycle limited to [0, 1]
    D = max(0, min(1, D));

    prev_error = error;
end

% ODE Function for Buck Converter with Battery
function dx = buckConverterBatteryODE(t, x, Vin, L, C, SOC, R_bat, Vbat_nominal, C_bat, D)
    iL = x(1);        % Inductor current
    vC = x(2);        % Capacitor voltage

    % Terminal voltage of battery based on SOC
    Vbat = Vbat_nominal * (1 - SOC) + R_bat * iL;

    % Differential equations
    diL_dt = (Vin * D - vC) / L;
    dvC_dt = (iL - vC / R_bat) / C;

    % SOC changes based on the inductor current (charging current)
    dSOC_dt = iL / C_bat;

    % Output vector (the derivatives)
    dx = [diL_dt; dvC_dt; dSOC_dt];
end

% Constraints Function
function [c, ceq] = constraintsFunction(pid_gains)
    global Vin L C R_bat Vbat_nominal C_bat SOC_initial Vout_ref t_sim dt

    % Unpack PID gains
    Kp = pid_gains(1);
    Ki = pid_gains(2);
    Kd = pid_gains(3);

    % Initial Conditions
    iL0 = 0;
    vC0 = 0;
    SOC0 = SOC_initial;
    x0 = [iL0; vC0; SOC0];  

    % Time vector for simulation
    time = 0:dt:t_sim;

    % Maximum allowable inductor current
    I_L_max = 10;  % Example constraint, 10 A maximum

    % Initialize variables for storing results
    max_iL = 0;
    max_D = 0;

    % PID control variables
    error_integral = 0;
    prev_error = 0;

    % Simulation loop
    for k = 1:length(time)
        x = x0;

        % Duty cycle Calculation
        [D, error_integral, prev_error] = PIDController(x(2), Vout_ref, Kp, Ki, Kd, dt, error_integral, prev_error);

        % Computation of the next state including SOC
        [~, x_next] = ode45(@(t, x) buckConverterBatteryODE(t, x, Vin, L, C, x(3), R_bat, Vbat_nominal, C_bat, D), [0 dt], x);

        % Tracking fo the maximum inductor current and duty cycle
        max_iL = max(max_iL, x_next(end, 1));
        max_D = max(max_D, D);

        % initial conditions for the next iteration
        x0 = x_next(end, :);
    end

    % Inequality constraints (c <= 0 for feasible solutions)
    c = [max_iL - I_L_max;  % Inductor current should not exceed I_L_max
         max_D - 1];        % Duty cycle should not exceed 1 (100%)
     
    % Equality constraint (ceq == 0 for feasible solutions)
    ceq = x0(3) - 0.95;  % Final SOC should be 0.95 (95% charged)
end


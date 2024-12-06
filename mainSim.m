%% mainSim script

% Created and mantained by
% Azevedo Da Silva Esteban
% Domenichelli Eleonora
% Donati Filippo 
% Gavidia Pantoja Maria Paulina

% We should consider to move the data generation into separare
% sub-functions as this will greatly improve readability


%% Falgs

flag.realsensors = 1; % Set to 1 to use real sensor models, with mesuring errors
flag.sensverb = 0; % Set to 1 to augment the number of data collected from sensors

%% Environment Data

% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Earth

env.Earth.n = 2*pi / (360*24*60*60); % [1x1] rad/2 - Earth mean Velocity
env.Earth.R = astroConstants(23); % [1x1] Km - Radius of the Earth
env.Earth.i = deg2rad(23.45); % [1x1] rad - Earth Rotation Axis Inclination

% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Sun

env.Sun.R = astroConstants(2); % [1x1] Km - Earth to Sun distance

%% Satellite Orbit Data

orb.a = 42000000; % [1x1] m - Semi-major axis 
orb.e = 0.1; % [1x1] - Eccentricity
orb.i = 0; % [1x1] rad - Inclination
orb.n = sqrt(6.6743E-11 * 5.972E24/(R^3)); % [1x1] rad/s - Mean orbital Velocity
orb.T = 2*pi / orb.n; % [1x1] s - Orbital Period

%% Satellite data

sat.Iv = [0.06;0.08;0.04]; % [3x1] Kgm^2 - Pincipal Inertial Moments as a Vector
sat.I = diag(sat.Iv); % [3x3] Kgm^2 - Inertia Matrix
sat.invI = inv(sat.I); % [3x3] Kgm^2 - Inverse Inertia Matrix

%% Sensor Data

% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Sun Sensor

% Generic Sensor Data
sens.ss.ADC.bit = 10;
sens.ss.fov = deg2rad(60); % Sensor FOV in radiants
sens.ss.freq = 50; % Sampling Frequency in Hz
sens.ss.accuracy = deg2rad(0.5); % Accuracy of sun sensor in radiants
sens.ss.precision = deg2rad(0.1); % Precision of the sensor in radiats
sens.ss.ADC.quanta = (2 * sens.ss.fov) / ((2)^(sens.ss.ADC.bit)); % ADC Quanta
% Note: current notation assume that ADC quanta on Voltage is linearly
% correlated with quantization of the angle. This isn't true, but good
% enough approximation

% Face specific Data

% +--------+
% | Face 1 |
% +--------+

sens.ss.S0.n_b = [1, 0, 0]';
sens.ss.S0.miss = randn(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.ss.S0.missalign = DCM(sens.ss.S0.miss); % Missalignement matrix
sens.ss.S0.A_ssb = eye(3); % Rotation matrix from body frame to Sensor frame
sens.ss.S0.bias.alpha = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S0.bias.beta = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S0.seed.beta = 0; % Random noise seed
sens.ss.S0.seed.alpha = 1; % Random Noise seed

% +--------+
% | Face 2 |
% +--------+


% +--------+
% | Face 3 |
% +--------+


% +--------+
% | Face 4 |
% +--------+

% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -

% Gyroscope

sens.gyro.freq = 2000; % Max value
sens.gyro.ADC.quanta = deg2rad(0.22)/3600;
sens.gyro.saturation = deg2rad(400);
sens.gyro.sf = randn(3,1) .* (500e-6);
sens.gyro.scaleFactor = eye(3) + diag(sens.gyro.sf);
sens.gyro.miss = randn(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.gyro.missalign = DCM(sens.gyro.miss); % Missalignement matrix
sens.gyro.orth = [[1, 1e-3, 1e-3]; [0, 1, 1e-3]; [0, 1, 1e-3]];
sens.gyro.bias = randn(3,1) .* (deg2rad(4)/3600);
sens.gyro.arw = deg2rad(0.15/60*sqrt(1/sens.gyro.freq)); % clearly wrong, but whatever 
sens.gyro.biasInstability = deg2rad(0.3)/3600; % NEED TO BE CONVERTED
sens.gyro.gsensing = eye(3);

% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -

%% Intial Conditions

IC.w0 = [1e-6; 1e-6; orb.n]; % [3x1] rad/s - Initial Angular rates
IC.angles = [0, 0.0000001, 0]; % [1x3] rad - Initial Euler Angles wrt ECI
IC.theta = 0; % [1x1] rad - Initial true anomaly of the Spacecraft

%% Simulation Options

simul.t0 = 0;
simul.tf = 3*orb.T;

%% Simulation Start

out = sim("RealWorld.slx", "StartTime", "simul.t0", "StopTime", "simul.tf");

%% Post-Processing

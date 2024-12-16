%% mainSim script

% Created and mantained by
% Azevedo Da Silva Esteban
% Domenichelli Eleonora
% Donati Filippo 
% Gavidia Pantoja Maria Paulina

% We should consider to move the data generation into separate
% sub-functions as this will greatly improve readability


%% Flags

flag.realsensors = 1; % Set to 1 to use real sensor models, with mesuring errors
flag.sensverb = 0; % Set to 1 to augment the number of data collected from sensors
flag.gg = 1; % Set to 1 to use gravity gradient perturbation
flag.srp = 1; % Set to 1 to consider srp perturbation
flag.mag = 1; % Set to 1 to consider magnetic field perturbation

%% Environment Data

env.c = astroConstants(5)*1000; % [1x1] m/s - Speed of light
env.G = astroConstants(1); % [1x1] km^3/(kg*s^2) - Universal gravity constant
% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Earth

env.Earth.n = 2*pi / (360*24*60*60); % [1x1] rad/2 - Earth mean Velocity
env.Earth.R = astroConstants(23); % [1x1] Km - Radius of the Earth
env.Earth.i = deg2rad(23.45); % [1x1] rad - Earth Rotation Axis Inclination
env.Earth.mu = astroConstants(13);
env.Earth.mass = env.Earth.mu/env.G; % [1x1] kg - Earth mass
env.Earth.omega = 7.2921159 * 1e-5; % [1x1] rad/s - Earth angular velocity
env.Earth.magInclination = deg2rad(11.5); % [1x1] rad - Magnetic field inclination
env.mag.dgrf2020 = [-29404.8; -1450.9; 4652.5].*1e-9; % [3x1] Teslas - Magnetic field costants
% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Sun

env.Sun.R = astroConstants(2); % [1x1] Km - Earth to Sun distance
env.Sun.Fe = 1358; % [1x1] W/m^2 - Solar radiation intensity

env.dis.P = env.Sun.Fe/env.c; % [1x1] kg/(m*s^2) - Average pressure due to radiation

%% Satellite Orbit Data

orb.a = 20000; % [1x1] m - Semi-major axis 
orb.e = 0; % [1x1] - Eccentricity
orb.i = 101; % [1x1] rad - Inclination
orb.n = sqrt(astroConstants(13)/(orb.a^3)); % [1x1] rad/s - Mean orbital Velocity
orb.T = 2*pi / orb.n; % [1x1] s - Orbital Period

orb.W = -2.0647E14 * orb.a^(-7/2) * cos(orb.i) * time; % [1x1] degrees/day - RAAN (for e = 0)

% IN GENERAL : orb.W = -3/2 * J2 * (env.Earth.R/orb.a)^2 * 1/(1 - orb.e^2) * sqrt(env.Earth.mu/orb.a^3) * cos(orb.i) * time; 


%% Satellite data

sat.Iv = [0.06;0.08;0.04]; % [3x1] Kgm^2 - Pincipal Inertial Moments as a Vector
sat.I = diag(sat.Iv); % [3x3] Kgm^2 - Inertia Matrix
sat.invI = inv(sat.I); % [3x3] Kgm^2 - Inverse Inertia Matrix
sat.n_b = [1, 0, 0; 0, 1, 0; -1, 0, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1]'; % [3x6] - Vector normal to the each surface of satellite in body frame
sat.A = [6*1e-2*ones(4,1); 4*1e-2*ones(2,1)]';% [1x6] m^2 - Surfaces (WRONG NUMBERS!!!)
sat.rho_s = 0.5*ones(6,1); % [6x1] - Surfaces' diffuse reflection coefficients (WRONG NUMBERS!!!)
sat.rho_d = 0.1*ones(1, 6); % [1x6] - Surfaces' specular reflection coefficients (WRONG NUMBERS!!!)
sat.r_CM = [10, 0, 15; 0, 10, 15; -10, 0, -15; 0, -10, -15; 0, 0, 30; 0, 0, 0]*1e-2; % [6x3] m - Distance from centre panel to CM (WRONG NUMBERS!!!)
sat.jB = [0.01; 0.05; 0.01]; %!!!!!!!

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

% Magentometer
sens.mag.f = 10;% [Hz]
sens.mag.SNR = 10^(70/20); % dB


%% Actuator

% Control Moment Gyroscope

% Generic Data 


% Gyro specific Data

% +--------+
% | Gyro 1 |
% +--------+



% +--------+
% | Gyro 2 |
% +--------+



% +--------+
% | Gyro 3 |
% +--------+



% +--------+
% | Gyro 4 |
% +--------+



%% Intial Conditions

IC.w0 = [1e-6; 1e-6; orb.n]; % [3x1] rad/s - Initial Angular rates
IC.angles = [0, 0.0000001, 0]; % [1x3] rad - Initial Euler Angles wrt ECI
IC.theta = 0; % [1x1] rad - Initial true anomaly of the Spacecraft

%% Simulation Options

simul.t0 = 0;
simul.tf = 3*orb.T;

%% Simulation Start

out = sim("CubeSat.slx", "StartTime", "simul.t0", "StopTime", "simul.tf");

%% Post-Processing

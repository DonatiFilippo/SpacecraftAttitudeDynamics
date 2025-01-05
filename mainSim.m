%% mainSim script

% Created and mantained by
% Azevedo Da Silva Esteban
% Domenichelli Eleonora
% Donati Filippo 
% Gavidia Pantoja Maria Paulina

% We should consider to move the data generation into separate
% sub-functions as this will greatly improve readability

%% Operational Mode

% Control and navigation modes 
% 0 - No control
% 1 - Sun Pointing / Slew
% 2 - Detumbling
op.mode = 1;

%% Flags

flag.realsensors = 1; % Set to 1 to use real sensor models, with mesuring errors
flag.realact = 1; % Set to 1 to use real sensor model
flag.gg = 1; % Set to 1 to use gravity gradient perturbation
flag.srp = 1; % Set to 1 to consider srp perturbation
flag.mag = 1; % Set to 1 to consider magnetic field perturbation
flag.aero = 1; % Set to 1 to consider atmospheric drag perturbation
flag.forcemode = 1; % Set to 1 to force the system into a particular control mode

%% Environment Data

env.c = astroConstants(5)*1000; % [1x1] m/s - Speed of light
env.G = astroConstants(1); % [1x1] km^3/(kg*s^2) - Universal gravity constant
env.atm.h0 = 1000;
env.atm.H = 268;
env.atm.rho0 = 3.019*10^(-15); 
% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Earth

env.Earth.n = 2*pi / (365*24*60*60); % [1x1] rad/2 - Earth mean Velocity
env.Earth.R = astroConstants(23); % [1x1] Km - Radius of the Earth
env.Earth.i = deg2rad(23.43928111); % [1x1] rad - Earth Rotation Axis Inclination
env.Earth.mu = astroConstants(13); % [1x1] km^3/s^2 - Earth gravity constat 
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

orb.a = env.Earth.R + 887.6790; % [1x1] Km - Semi-major axis 
orb.e = 0; % [1x1] - Eccentricity
orb.i = deg2rad(98.9897); % [1x1] rad - Inclination
orb.n = sqrt(astroConstants(13)/(orb.a^3)); % [1x1] rad/s - Mean orbital Velocity
orb.T = 2*pi / orb.n; % [1x1] s - Orbital Period
orb.OMprec = env.Earth.n; % [1x1] deg/d - RAAN precession for SSO

%% Satellite data

sat.Iv = [0.04327;0.095068;0.120327]; % [3x1] Kgm^2 - Pincipal Inertial Moments as a Vector
sat.I = diag(sat.Iv); % [3x3] Kgm^2 - Inertia Matrix
sat.invI = inv(sat.I); % [3x3] Kgm^2 - Inverse Inertia Matrix
sat.cd = 2.44; % [1x1] - Drag coefficient
sat.n_b = [1, 0, 0; 0, 1, 0; -1, 0, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1; 0, 0, 1; 0, 0, -1; 0, 0, 1; 0, 0, -1]'; % [3x10] - Vector normal to the each surface of satellite in body frame
sat.A = [26894, 82716, 26894, 82716, 43554, 43554, 23100, 23100, 23100, 23100]*1e-6;% [1x10] m^2 Surfaces
sat.rho_s = [0.8*ones(6,1); 0.35*ones(4,1)]; % [10x1] - Surfaces' diffuse reflection coefficients 
sat.rho_d = [0.45*ones(6,1);0.4*ones(4,1)]; % [10x1] - Surfaces' specular reflection coefficients
sat.r_CM = [113, 0, 1; 0, 59.5, 1; -113, 0, 1; 0, -59.5, 1; 0, 0, 184; 0, 0, -182; 218, 0, -182; 218, 0, -182; -218, 0, -182; -218, 0, -182]*1e-3; % [6x3] m - Distance from centre panel to CM 
sat.jB = [0.27; 0.53; 0.78]*0.01; % [3x1] Am^2 - Dipole moment of the spacecraft

%% Sensor Data

% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Sun Sensor

% Generic Sensor Data
sens.ss.ADC.bit = 10;
sens.ss.fov = deg2rad(60); % Sensor FOV in radiants
sens.ss.freq = 50; % Sampling Frequency in Hz
sens.ss.accuracy = deg2rad(0.5); % Accuracy of sun sensor in radiants
sens.ss.precision = deg2rad(0.1); % Precision of the sensor in radiants
sens.ss.ADC.quanta = (2 * sens.ss.fov) / ((2)^(sens.ss.ADC.bit)); % ADC Quanta
% Note: current notation assume that ADC quanta on Voltage is linearly
% correlated with quantization of the angle. This isn't true, but good
% enough approximation

% Face specific Data

% +--------+
% | Fine 0 |
% +--------+

sens.ss.S0.n_b = [0, 0, 1]'; % To change when new matrix K is made
sens.ss.S0.miss = sign(randn(3,1)).*rand(3,1) * deg2rad(1e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.ss.S0.missalign = DCM(sens.ss.S0.miss); % Missalignement matrix
sens.ss.S0.A_ssb = [[0, 0, -1]; [0, 1, 0]; [1, 0, 0]]; % Rotation matrix from body frame to Sensor frame
sens.ss.S0.bias.alpha = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S0.bias.beta = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S0.seed.beta = 0; % Random noise seed
sens.ss.S0.seed.alpha = 1; % Random Noise seed


% +--------+
% | Face 1 |
% +--------+

sens.ss.S1.n_b = [1, 0, 0]';
sens.ss.S1.miss = sign(randn(3,1)).*rand(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.ss.S1.missalign = DCM(sens.ss.S1.miss); % Missalignement matrix
sens.ss.S1.A_ssb = eye(3); % Rotation matrix from body frame to Sensor frame
sens.ss.S1.bias.alpha = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S1.bias.beta = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S1.seed.beta = 2; % Random noise seed
sens.ss.S1.seed.alpha = 3; % Random Noise seed

% +--------+
% | Face 2 |
% +--------+

sens.ss.S2.n_b = [0, 1, 0]';
sens.ss.S2.miss = sign(randn(3,1)).*rand(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.ss.S2.missalign = DCM(sens.ss.S2.miss); % Missalignement matrix
sens.ss.S2.A_ssb = [[0, 1, 0]; [-1, 0, 0]; [0, 0, 1]]; % Rotation matrix from body frame to Sensor frame
sens.ss.S2.bias.alpha = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S2.bias.beta = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S2.seed.beta = 4; % Random noise seed
sens.ss.S2.seed.alpha = 5; % Random Noise seed

% +--------+
% | Face 3 |
% +--------+

sens.ss.S3.n_b = [-1, 0, 0]';
sens.ss.S3.miss = sign(randn(3,1)).*rand(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.ss.S3.missalign = DCM(sens.ss.S3.miss); % Missalignement matrix
sens.ss.S3.A_ssb = [[-1, 0, 0]; [0, -1, 0]; [0, 0, 1]]; % Rotation matrix from body frame to Sensor frame
sens.ss.S3.bias.alpha = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S3.bias.beta = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S3.seed.beta = 6; % Random noise seed
sens.ss.S3.seed.alpha = 7; % Random Noise seed

% +--------+
% | Face 4 |
% +--------+

sens.ss.S4.n_b = [0, -1, 0]';
sens.ss.S4.miss = sign(randn(3,1)).*rand(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.ss.S4.missalign = DCM(sens.ss.S4.miss); % Missalignement matrix
sens.ss.S4.A_ssb = [[0, -1, 0]; [1, 0, 0]; [0, 0, 1]]; % Rotation matrix from body frame to Sensor frame
sens.ss.S4.bias.alpha = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S4.bias.beta = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S4.seed.beta = 8; % Random noise seed
sens.ss.S4.seed.alpha = 9; % Random Noise seed


% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + -

% Magentometer
% Hp: magnetometer assumed oriented as the principal inertia frame
% datas recovered from the data sheet
% put everything in nT 1 G = 1e5 nT

f = 10; % [Hz]
sens.mag.SNR = 10^(70/20); % dB 

% Misalignment matrix computation using small angles approximation
ang = deg2rad(sign(randn(3,1)).*rand(3,1)); % rad
sens.mag.A_mis = DCM(ang); 

% Non-orthogonality matrix computation, considering the cross-axis
% sensitivity (cxs) from the data sheet
cxs = 0.002;
rng('default');
s_xy = cxs*(2*rand-1);
s_xz = cxs*(2*rand-1);
s_yx = 0;
s_yz = cxs*(2*rand-1);
s_zx = 0;
s_zy = 0;
sens.mag.A_nonorth = [1, s_xy, s_xz; 
                      s_yx, 1, s_yz;
                      s_zx, s_zy, 1];


sens.mag.Ts = 1/(2*f);
sens.mag.Quant = 7*1e-7;
sens.mag.sat = [4, -4] .* 1e-4;


%% Actuator

% Control Moment Gyroscope

% Generic Data 
act.cmg.h = [1.4;1.4;1.4;1.4]*1e-3; % [1x1] Nm^2s - Angular momentum of each gyro
act.cmg.sat = 9e-3; % [1x1] N - Max Torque that can be produced by the cmg
act.cmg.beta = 54.73; % [1x1] deg - Skew angle nominal 
act.cmg.maxgimbrate = 10; % [1x1] rad/s - Max gimble rate
act.cmg.quanta = deg2rad(0.0879); % [1x1] deg - Encoder resolution
act.pwm.resolution = 0.010702; % [1x1] rad - PWM generator maximum 
act.cmg.err= sign(randn(4,1)).*rand(4,1); % [4x1] deg - Missalignment error of the gyros
act.cmg.betareal = act.cmg.beta*ones(4,1) + act.cmg.err; % [4x1] deg - Real betas

%% OBC Settings
obc.freq = 25;

%% Navigation

 nav.alpha1 = 1/sens.ss.accuracy^2;
 nav.alpha2 = sens.mag.SNR / (1 + norm(sens.mag.A_nonorth, 'fro'));

%nav.alpha1 = 0.75;
%nav.alpha2 = 0.25;

% State-Observer
nav.ws = 0.05;
ws = nav.ws;

A = [[0,(sat.Iv(2)-sat.Iv(3))/sat.Iv(1)*ws,  0]; [(sat.Iv(3)-sat.Iv(1))/sat.Iv(2)*ws, 0, 0]; [0, 0, 0]];
nav.B2= sat.invI;

C = [1, 0 , 0;
    0, 1, 0;
    0, 0, 1]; 

R = diag([1; 1; 3]);
Q = diag([4.5e-2; 4.5e-2; 2e-3]);

%R = diag([1; 3; 2.5]);
%Q = diag([4.5e-2; 4.5e-2; 2e-6]);

%Q = diag([5; 5; 1e-3]);
%R = diag([0.5; 1; 0.8]);

[K2,S,P] = lqr(A',C,Q,R);

nav.L = K2';
nav.A2 = A - nav.L*C;



%% Guidance

% Sun poiting 
guid.wsg = ws;
% wsg = 0.01;
%% Control
wsc = ws;
%wsc = 60 * orb.n;

Ac = [0, wsc, 0, -1, 0;
    -wsc, 0, 1, 0, 0;
    0, 0, 0, (sat.Iv(2)-sat.Iv(3))/sat.Iv(1)*wsc, 0;
    0, 0, (sat.Iv(3)-sat.Iv(1))/sat.Iv(2)*wsc, 0, 0;
    0, 0, 0, 0, 0];

Bc = [0, 0, 0;
    0, 0, 0;
    1/sat.Iv(1), 0, 0;
    0, 1/sat.Iv(2), 0;
    0, 0, sat.Iv(3)];

% s_x_max = 3;
% s_y_max = 3;
% w_x_max = 10^-2;
% w_y_max = 10^-2;
% w_z_max = 10^-3;

s_x_max = 3e-1;
s_y_max = 3e-1;
w_x_max = 6e-1;
w_y_max = 6e-1;
w_z_max = 4e-3;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2])*1.115;


u_x_max = 9e-3;
u_y_max = 9e-3;
u_z_max = 9e-3;
Rc = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2])*7e-1;


% Another option
% Qc = diag([7e-5; 7e-5; 1.2; 1.2; 0.35]);
% Rc = diag([1; 1; 1]);

[control.K,~,~] = lqr(Ac,Bc,Qc,Rc);

% Detumbling

control.kdet = [[-0.0272, 0, 0]; [0, -0.0272, 0]; [0, 0, -0.0272]];

%% Intial Conditions

IC.w0 = [0; 0; 1e-5]; % [3x1] rad/s - Initial Angular rates
IC.angles = [0.05, 0.05, 0.05]; % [1x3] rad - Initial Euler Angles wrt ECI
IC.theta = 0; % [1x1] rad - Initial true anomaly of the Spacecraft
IC.OM = deg2rad(-90); % Intial RAAN

%% Simulation Options

simul.t0 = 0;
simul.tf = orb.T;

%% Simulation Start

out = sim("CubeSat.slx", "StartTime", "simul.t0", "StopTime", "simul.tf");

%% Post-Processing
% Perturbations representation along the orbit

% SRP
figure
plot(out.tout, out.T_srp(:,1),'r-', out.tout, out.T_srp(:,2), 'g-', out.tout, out.T_srp(:,3), 'b-', 'LineWidth', 2)
grid on
legend('Tx_s_r_p', 'Ty_s_r_p', 'Tz_s_r_p')
xlabel('Time [s]')
ylabel('Torque [Nm]')
title('Solar Radiation Pressure torque')
xlim([0, orb.T])

% GG
figure
plot(out.tout, out.T_gg(:,1),'r-', out.tout, out.T_gg(:,2), 'g-', out.tout, out.T_gg(:,3), 'b-', 'LineWidth', 2)
grid on
legend('Tx_g_g', 'Ty_g_g', 'Tz_g_g')
xlabel('Time [s]')
ylabel('Torque [Nm]')
title('Gravity Gradient torque')
xlim([0, orb.T])

% MAG
figure
plot(out.tout, out.T_mag(:,1),'r-', out.tout, out.T_mag(:,2), 'g-', out.tout, out.T_mag(:,3), 'b-', 'LineWidth', 2)
grid on
legend('Tx_m_a_g', 'Ty_m_a_g', 'Tz_m_a_g')
xlabel('Time [s]')
ylabel('Torque [Nm]')
title('Magnetic field torque')
xlim([0, orb.T])

clear;
clc;
close all;

%% DATAS

% DYNAMICS
J = diag([0.0823, 0.0633, 0.0317]); %hypotized inertia matrix generated from ChatGPT known our chosen s/c
invJ = inv(J);

% ENVIRONMENT
env.orbit.i = deg2rad(97.36);
env.orbit.e = 0.0005568;
env.orbit.a = 6844;
env.orbit.T = 93.91*60;
env.orbit.n = 2*pi/env.orbit.T;
env.Earth.omega = 7.2921159 * 1e-5;
env.Earth.magInclination = deg2rad(11.5);
env.Earth.R = 6371000;
env.mag.dgrf2020 = [-29404.8; -1450.9; 4652.5].*1e-9;

% MAGNETOMETER
% Hp: magnetometer assumed oriented as the principal inertia frame
% datas recovered from the data sheet
% put everything in nT 1 G = 1e5 nT

f = 10; % [Hz]
sens.mag.SNR = 10^(70/20); % dB 

% Misalignment matrix computation using small angles approximation
ang = deg2rad([0.5, -0.3, 0.2]); % rad
sens.mag.A_mis = [1, ang(3), -ang(2);
                 -ang(3), 1, ang(1);
                 ang(2), -ang(1), 1]; % misalignment matrix, small angles approximation

% Non-orthogonality matrix computation, considering the cross-axis
% sensitivity (cxs) from the data sheet
cxs = 0.002;
rng('default');
s_xy = cxs*(2*rand-1);
s_xz = cxs*(2*rand-1);
s_yx = 0;
s_yz = cxs*(2*rand-1);
s_zx = 0;
s_zy = cxs*(2*rand-1);
sens.mag.A_nonorth = [1, s_xy, s_xz; 
                      s_yx, 1, s_yz;
                      s_zx, s_zy, 1];

sens.mag.flag = 1;
sens.mag.Ts = 1/(2*f);
sens.mag.Quant = 7*1e-7;
sens.mag.sat = [4, -4] .* 1e-4;


% INITIAL CONDITIONs
IC.w = [0.45, 0.52, 0.55];
IC.A = eye(3);
IC.theta = 0;

%% Simulation
t0 = 0;
tf = 100;

%% Data

%Falgs
flag.realsensors = 1;
flag.sensverb = 0; %activate debug option 

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
sens.ss.S0.n_b = [1, 0, 0]';
sens.ss.S0.miss = randn(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.ss.S0.missalign = sens.mag.A_mis; % Missalignement matrix
sens.ss.S0.A_ssb = eye(3); % Rotation matrix from body frame to Sensor frame
sens.ss.S0.bias.alpha = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S0.bias.beta = randn * (sens.ss.accuracy/3); % Accuracy of sun sensor in radiants
sens.ss.S0.seed.beta = 0; % Random noise seed
sens.ss.S0.seed.alpha = 1; % Random Noise seed

% ORBIT DATA
theta0 = 0;
e = 0.0005568;
a = 6844;
i = deg2rad(97.36);
G = 6.6743E-11;
Mt = 5.972E24;
T = 93.91*60;
n = 2*pi/T;

% KINEMATICS DATA
Ix = 0.0823;
Iy = 0.0633;
Iz = 0.0317;
wx0 = 0.45;
wy0 = 0.52;
wz0 = 0.55;

% KINEMATCIS QUANTITIES
J = diag([0.0823, 0.0633, 0.0317]);
w0 = [wx0; wy0; wz0];
ABN0 = diag([1, 1, 1]);
wLN = [0; 0; n];

% MAGNETIC TORQUE DATA
phi = deg2rad(11.5);
g01 = -29404.8E-9;
g11 = -1450.0E-9;
h11 = 4652.5E-9;
H0 = sqrt(g01^2 + g11^2 + h11^2);
wE = 7.2921159 * 1e-5;
RE = 6371000;
jB = [0.01; 0.05; 0.01];

% Sun Direction
Ts = 86400*365.25;
ns = 2*pi/Ts;
epsi = deg2rad(23.45);
es = 0.0167;
as = 149.5978707e11;

% Simulation 
t0 = 0;
tf = T;
out = sim("Q_method.slx", "StartTime", "t0", "StopTime", "tf", "FixedStep", "0.1");
%% LAB 8
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
env.Earth.R = astroConstants(23);
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
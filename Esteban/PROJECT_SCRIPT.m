clc
close all
clear all

% ORBIT DATA
theta0 = 0;
e = 0.1;
a = 42000000;
i = 0.5;
G = 6.6743E-11;
Mt = 5.972E24;
n = sqrt(G*Mt/a^3);
T = 2*pi/n;

% KINEMATICS DATA
Ix = 0.06;
Iy = 0.08;
Iz = 0.04;
wx0 = 1E-6;
wy0 = 1E-6;
wz0 = n;

% KINEMATCIS QUANTITIES
J = diag([Ix, Iy, Iz]);
w0 = [wx0; wy0; wz0];
ABN0 = diag([1, 1, 1]);
wLN = [0; 0; n];

% MAGNETIC TORQUE DATA
phi = deg2rad(11.5);
g01 = -29404.8E-9;
g11 = -1450.0E-9;
h11 = 4652.5E-9;
H0 = sqrt(g01^2 + g11^2 + h11^2);
wE = n;
RE = 6378000;
jB = [0.01; 0.05; 0.01];

% Simulation 
t0 = 0;
tf = T;
out = sim("PROJECT_SIMULINK.slx", "StartTime", "t0", "StopTime", "tf", "FixedStep", "0.1");
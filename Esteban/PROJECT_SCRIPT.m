clc
close all
clear all

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
wE = n;
RE = 6378000;
jB = [0.01; 0.05; 0.01];

% Simulation 
t0 = 0;
tf = 10;
out = sim("test.slx", "StartTime", "t0", "StopTime", "tf", "FixedStep", "0.1");
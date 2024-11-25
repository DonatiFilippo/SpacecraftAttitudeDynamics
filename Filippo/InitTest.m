clc 
close all
clear all

%% Data
Iv = [0.06;0.08;0.04];
R = 42000000;
orb.n = sqrt(6.6743E-11 * 5.972E24/(R^3));
env.Earth.n = 2*pi / (360*24*60*60);
T = 2*pi / orb.n;
orb.e = 0.1;
orb.i = 0;
flag.gg = 1;
orb.a = R;
env.Earth.i = deg2rad(23.45);
env.Earth.R = astroConstants(23);
env.mass = 5.972E24;
env.G = 6.6743E-11;
env.Sun.R = astroConstants(2);
I = diag(Iv);
invI = inv(I);
w0 = [1e-6; 1e-6; orb.n];
ABN0 = diag([1, 1, 1]);
IC.theta = 0;
IC.angles = [0, 0.0000001, 0]; %The system is smart and choose Euler rappresentation based on singularities
%% Sim config
t0 = 0;
tf = T;
%% Simulation Init
out = sim("testModel.slx", "StartTime", "t0", "StopTime", "tf");
clc 
close all
clear all

%% Data
Iv = [0.07;0.055;0.025];

I = diag(Iv);
invI = inv(I);
w0 = [0; 0; 0];
A = [[1, 0, 0]; [0, sqrt(2)/2, sqrt(2)/2]; [0, -sqrt(2)/2, sqrt(2)/2]];
IC.angles = [[1, 0 , 0]; [0, 1, 0]; [0, 0, 1]]; %The system is smart and choose Euler rappresentation based on singularities
env.R = 408000;
env.G = 6.6743e-11;
env.mass = 5.9722e24;
%% Sim config
t0 = 0;
tf = 10000;
%% Simulation Init
out = sim("Lab6.slx", "StartTime", "t0", "StopTime", "tf");

%% Plot
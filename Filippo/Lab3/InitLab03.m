clc 
close all
clear all

%% Data
Iv = [0.07;0.055;0.025];

I = diag(Iv);
invI = inv(I);
w0 = [2*pi; 0.5; 1];
%% Sim config
t0 = 0;
tf = 5000;
%% Simulation Init
out = sim("Lab3.slx", "StartTime", "t0", "StopTime", "tf");

%% Verification
plot(out.omega)
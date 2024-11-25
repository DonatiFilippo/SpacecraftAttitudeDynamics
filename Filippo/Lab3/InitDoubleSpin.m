clc 
close all
clear all

%% Data
Iv = [0.07;0.025;0.055];

I = diag(Iv);
invI = inv(I);
w0 = [1e-6; 0.005; 0.2];
wr0 = 2*pi;
Jv = [0; 0; 0.005];

J = diag(Jv);
invJ = inv(J);
torque = 0;
%% Sim config
t0 = 0;
tf = 5000;
%% Simulation Init
out = sim("DoubleSpin.slx", "StartTime", "t0", "StopTime", "tf");

%% Verification
plot(out.omega)
clear
close all
clear all

%% Data

mu = 1.327e11;
r0 = [2.029e8; -1.475e5; -1.1395e7];
v0 = [3.021; 24.139; 10.964];

%% Simulation Setting
t0 = 0;
tf = 864000;
out = sim("TwoBodySystem.slx","StartTime","t0","StopTime","tf");

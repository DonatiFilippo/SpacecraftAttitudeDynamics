clear
close all
clear all

%% Data

m = 2;
b = 5;
k = 6;
omega = 0.5;
x0 = 1;
v0 = 0;

%% Simulation Setting
t0 = 0;
tf = 200;
out = sim("Damper.slx","StartTime","t0","StopTime","tf");

%% Plot

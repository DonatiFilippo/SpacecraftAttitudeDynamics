clear
close all
clear all

%% Data

g = 9.81;
l = 5;
theta0 = 0;
theta_dot0 = pi/4;

%% Simulation Setting
t0 = 0;
tf = 200;
out = sim("Lab2.slx","StartTime","t0","StopTime","tf");

%% Plot

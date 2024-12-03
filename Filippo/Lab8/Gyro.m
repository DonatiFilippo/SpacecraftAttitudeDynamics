clc 
close all
clear all

%% Data

%Falgs
flag.realsensors = 1;
flag.sensverb = 0; %activate debug option

sens.gyro.freq = 2000; % Max value
sens.gyro.ADC.quanta = deg2rad(0.22)/3600;
sens.gyro.saturation = deg2rad(400);
sens.gyro.sf = randn(3,1) .* (500e-6);
sens.gyro.scaleFactor = eye(3) + diag(sens.gyro.sf);
sens.gyro.miss = randn(3,1) * deg2rad(3e-2); % Missaligment angles
% I have choosen a random value of 1'' as missaligment 1sigma
sens.gyro.missalign = DCM(sens.gyro.miss); % Missalignement matrix
sens.gyro.orth = [[1, 1e-3, 1e-3]; [0, 1, 1e-3]; [0, 1, 1e-3]];
sens.gyro.bias = randn(3,1) .* (deg2rad(4)/3600);
sens.gyro.arw = deg2rad(0.15/60*sqrt(1/sens.gyro.freq)); % clearly wrong, but whatever 
sens.gyro.biasInstability = deg2rad(0.3)/3600; % NEED TO BE CONVERTED
sens.gyro.gsensing = eye(3);
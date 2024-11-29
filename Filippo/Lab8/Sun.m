clc 
close all
clear all

%% Data

% Generic Sensor Data
sens.ss.fov = deg2rad(60); % Sensor FOV in radiants
sens.ss.freq = 50; % Sampling Frequency in Hz
sens.ss.precision = deg2rad(0.1); % Precision of the sensor in radiats
sens.ss.ADC.quanta = 0.001; % ADC Quanta
% Note: current notation assume that ADC quanta on Voltage is linearly
% correlated with quantization of the angle. This isn't true, but good
% enough approximation

% Face specific Data
sens.ss.S0.missalign = eye(3); % Missalignement matrix
sens.ss.S0.A_ssb = eye(3) % Rotation matrix from body frame to Sensor frame
sens.ss.S0.bias.alpha = 0; % Accuracy of sun sensor in radiants
sens.ss.S0.bias.beta = 0; % Accuracy of sun sensor in radiants
sens.ss.S0.seed.beta = 0; % Random noise seed
sens.ss.S0.seed.alpha = 1; % Random Noise seed
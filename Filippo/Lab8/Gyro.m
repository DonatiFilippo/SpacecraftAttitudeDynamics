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
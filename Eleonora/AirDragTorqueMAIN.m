%% AIR DRAG
clear;
clc;
close all;

%% mainSim script

% Created and mantained by
% Azevedo Da Silva Esteban
% Domenichelli Eleonora
% Donati Filippo 
% Gavidia Pantoja Maria Paulina

% We should consider to move the data generation into separate
% sub-functions as this will greatly improve readability


%% Flags

flag.realsensors = 1; % Set to 1 to use real sensor models, with mesuring errors
flag.sensverb = 0; % Set to 1 to augment the number of data collected from sensors
flag.gg = 1; % Set to 1 to use gravity gradient perturbation
flag.srp = 1; % Set to 1 to consider srp perturbation
flag.mag = 1; % Set to 1 to consider magnetic field perturbation

%% Environment Data

env.c = astroConstants(5)*1000; % [1x1] m/s - Speed of light
env.G = astroConstants(1); % [1x1] km^3/(kg*s^2) - Universal gravity constant
env.mag.dgrf2020 = [-29404.8; -1450.9; 4652.5].*1e-9; % [3x1] Teslas - Magnetic field costants
env.atm.H = 60.828;
env.atm.h0 = 450;
env.atm.rho0 = 1.585 * 1e-12;
% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Earth

env.Earth.n = 2*pi / (360*24*60*60); % [1x1] rad/2 - Earth mean Velocity
env.Earth.R = astroConstants(23); % [1x1] Km - Radius of the Earth
env.Earth.i = deg2rad(23.45); % [1x1] rad - Earth Rotation Axis Inclination
env.Earth.mass = astroConstants(13)/env.G; % [1x1] kg - Earth mass
env.Earth.omega = 7.2921159 * 1e-5; % [1x1] rad/s - Earth angular velocity
env.Earth.magInclination = deg2rad(11.5); % [1x1] rad - Magnetic field inclination

% - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + %

% Sun

env.Sun.R = astroConstants(2); % [1x1] Km - Earth to Sun distance
env.Sun.Fe = 1358; % [1x1] W/m^2 - Solar radiation intensity

env.dis.P = env.Sun.Fe/env.c; % [1x1] kg/(m*s^2) - Average pressure due to radiation

%% Satellite Orbit Data

orb.a = env.Earth.R + 500; % [1x1] km - Semi-major axis 
orb.e = 0; % [1x1] - Eccentricity
orb.i = 0; % [1x1] rad - Inclination
orb.n = sqrt(astroConstants(13)/(orb.a^3)); % [1x1] rad/s - Mean orbital Velocity
orb.T = 2*pi / orb.n; % [1x1] s - Orbital Period

%% Satellite data

sat.Iv = [0.07; 0.055; 0.025]; % [3x1] Kgm^2 - Pincipal Inertial Moments as a Vector
sat.I = diag(sat.Iv); % [3x3] Kgm^2 - Inertia Matrix
sat.invI = inv(sat.I); % [3x3] Kgm^2 - Inverse Inertia Matrix
sat.n_b = [1, 0, 0; 0, 1, 0; -1, 0, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1]'; % [3x6] - Vector normal to the each surface of satellite in body frame
sat.A = [6*1e-2*ones(4,1); 4*1e-2*ones(2,1)]';% [1x6] m^2 - Surfaces (WRONG NUMBERS!!!)
sat.rho_s = 0.5*ones(6,1); % [6x1] - Surfaces' diffuse reflection coefficients (WRONG NUMBERS!!!)
sat.rho_d = 0.1*ones(1, 6); % [1x6] - Surfaces' specular reflection coefficients (WRONG NUMBERS!!!)
sat.r_CM = [
    0.10,  0,  0;  
    0,  0.10,  0;  
   -0.10,  0,  0;  
    0, -0.10,  0;  
    0,  0,  0.15;  
    0,  0, -0.15]; % [6x3] m - Distance from centre panel to CM (WRONG NUMBERS!!!)
sat.jB = [0.01; 0.05; 0.01]; %!!!!!!!
sat.Cd = 1.942;

%% Intial Conditions

IC.w0 = [0.45 0.52 0.55]; % [3x1] rad/s - Initial Angular rates
IC.angles = [0, 0.0000001, 0]; % [1x3] rad - Initial Euler Angles wrt ECI
IC.theta = 0; % [1x1] rad - Initial true anomaly of the Spacecraft

%% Simulation Options

simul.t0 = 0;
simul.tf = orb.T;

%% Simulation Start

out = sim("AirDragTorque.slx", "StartTime", "simul.t0", "StopTime", "simul.tf");



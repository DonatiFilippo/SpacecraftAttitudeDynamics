clc
clear all

%% Data

Ix = 0.07;
Iy = 0.055;
Iz = 0.025;

J = diag([Ix Iy Iz]);
W0 = [0.45 0.52 0.55];

G = astroConstants(1) ;
M = astroConstants(13)/G;
R = astroConstants(23) ;
n = sqrt(G*M/(R^3));
T = 365*24*60*60;
nx= 2*pi()/T;
r_sun = astroConstants(2);
Ep =  deg2rad(23.45);

c = astroConstants(5)

[a e i theta0] = deal(9500, 0.3, 0, 0);

N_B = [ 1,  0,  0; 
        0,  1,  0; 
       -1,  0,  0; 
       0,  -1,  0; 
       0,   0,  1; 
       0,   0, -1;
       1,   0,  0; 
       -1,  0,  0;  
       1,   0,  0;  
       -1,  0,  0];

rho_s = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5; 0.1; 0.1; 0.1; 0.1];
rho_d = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
Ar = [6e-2; 6e-2; 6e-2; 6e-2; 4e-2; 4e-2; 12e-2; 12e-2; 12e-2; 12e-2];
r_Fi = [
    0.10,  0,  0;  
    0,  0.10,  0;  
   -0.10,  0,  0;  
    0, -0.10,  0;  
    0,  0,  0.15;  
    0,  0, -0.15;  
    0,  0.45,  0;  
    0,  0.45,  0;  
    0,  -0.45,  0;  
    0, -0.45, 0
    ];



%% Simulation

t0 = 0;
tf = 100;
out = sim("Lab7_Kep.slx", "StartTime", "t0", "StopTime", "tf", "FixedStep", "0.1");
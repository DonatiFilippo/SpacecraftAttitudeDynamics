sat.Iv = [0.04327;0.095068;0.120327];

env.Earth.R = astroConstants(23);
orb.a = env.Earth.R + 887.6790; % [1x1] Km - Semi-major axis 
orb.e = 0; % [1x1] - Eccentricity
orb.i = deg2rad(103.39); % [1x1] rad - Inclination
orb.n = sqrt(astroConstants(13)/(orb.a^3)); % [1x1] rad/s - Mean orbital Velocity
orb.T = 2*pi / orb.n; % [1x1] s - Orbital Period

ws = 0.02;

A = [0, ws, 0, -1, 0;
    -ws, 0, 1, 0, 0;
    0, 0, 0, (sat.Iv(2)-sat.Iv(3))/sat.Iv(1)*ws, 0;
    0, 0, (sat.Iv(3)-sat.Iv(1))/sat.Iv(2)*ws, 0, 0;
    0, 0, 0, 0, 0];

B = [0, 0, 0;
    0, 0, 0;
    1/sat.Iv(1), 0, 0;
    0, 1/sat.Iv(2), 0;
    0, 0, sat.Iv(3)];

%Q = diag([7e-5; 7e-5; 1.2; 1.2; 7e-5]);
%R = diag([0.9; 2; 1]);

Q = diag([7e-5; 7e-5; 0.8; 0.8; 7e-7]);
R = diag([0.009; 0.02; 0.01]);

[K,S,P] = lqr(A,B,Q,R);
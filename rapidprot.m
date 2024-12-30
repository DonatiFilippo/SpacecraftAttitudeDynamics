sat.Iv = [0.04327;0.095068;0.120327];

env.Earth.R = astroConstants(23);
orb.a = env.Earth.R + 700; % [1x1] Km - Semi-major axis 
orb.e = 0; % [1x1] - Eccentricity
orb.i = deg2rad(103.39); % [1x1] rad - Inclination
orb.n = sqrt(astroConstants(13)/(orb.a^3)); % [1x1] rad/s - Mean orbital Velocity
orb.T = 2*pi / orb.n; % [1x1] s - Orbital Period

ws = 100*orb.n;

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

Q = diag([7e-5; 7e-5; 1.2; 1.2; 7e-5]);
R = diag([0.9; 2; 1]);

[K,S,P] = lqr(A,B,Q,R);

%% State observer
A = [[0,(sat.Iv(2)-sat.Iv(3))/sat.Iv(1)*ws,  0]; [(sat.Iv(3)-sat.Iv(1))/sat.Iv(2)*ws, 0, 0]; [0, 0, 0]];
C = eye(3,3);
R = diag([0.9; 2; 1]);
Q = diag([1.2; 1.2; 7e-5]);

[K2,S,P] = lqr(A',C,Q,R);

L = K2';
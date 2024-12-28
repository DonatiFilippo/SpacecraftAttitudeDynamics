%% LQR CONTROL

sat.Iv = [0.04327;0.095068;0.120327];

env.Earth.R = astroConstants(23);
orb.a = env.Earth.R + 887.6790; % [1x1] Km - Semi-major axis 
orb.e = 0; % [1x1] - Eccentricity
orb.i = deg2rad(98.9897); % [1x1] rad - Inclination
orb.n = sqrt(astroConstants(13)/(orb.a^3)); % [1x1] rad/s - Mean orbital Velocity
orb.T = 2*pi / orb.n;

ws = 10*orb.n;

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

% Try with the optimal control (marco mood)
% With this option the control on wz is really strong (goes to ws
% in like 0.2 secs but don't know how to slow it down). Also the problem is
% if I sow it down, it oscillates more (important to understand why wz
% oscillates, not much but it oscillates).
% The pointing error with actuators on oscillates, probably due to
% actuators saturation.
% s_x_max = 1;
% s_y_max = 1;
% w_x_max = 10^-2;
% w_y_max = 10^-2;
% w_z_max = 10^-2;
% Q = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2]);
% 
% u_x_max = 9e-3;
% u_y_max = 9e-3;
% u_z_max = 9e-3;
% R = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2]);


% Another option
% The two are really similar :(
% Q = diag([7e-5; 7e-5; 1.2; 1.2; 1.2]);
% R = diag([1; 1; 1]);

% Here we have a slightly less strong control on wz, but oscillates more
Q = diag([7e-5; 7e-5; 1.2; 1.2; 1.2]);
R = diag([1; 1; 1000]);

[K,S,P] = lqr(A,B,Q,R);

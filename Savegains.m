s_x_max = 1e-3;
s_y_max = 1e-3;
w_x_max = 2;
w_y_max = 2;
w_z_max = 2;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2]);


u_x_max = 9e-5;
u_y_max = 9e-5;
u_z_max = 9e-4;
Rc = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2]);

%% B
s_x_max = 5e-1;
s_y_max = 5e-1;
w_x_max = 1e-1;
w_y_max = 1e-1;
w_z_max = 1e-2;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2]);


u_x_max = 9e-3;
u_y_max = 9e-3;
u_z_max = 9e-3;
Rc = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2]);

%% Ele
s_x_max = 3;
s_y_max = 3;
w_x_max = 1e-2;
w_y_max = 1e-2;
w_z_max = 1e-2;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2]);


u_x_max = 9e-3;
u_y_max = 9e-3;
u_z_max = 9e-3;
Rc = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2]);

%% New Try
s_x_max = 2e-1;
s_y_max = 2e-1;
w_x_max = 6e-1;
w_y_max = 6e-1;
w_z_max = 5e-3;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2])*1e1;


u_x_max = 9e-3;
u_y_max = 9e-3;
u_z_max = 9e-3;
Rc = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2])*5e-2;

%% I like this

s_x_max = 2e-1;
s_y_max = 2e-1;
w_x_max = 6e-1;
w_y_max = 6e-1;
w_z_max = 5e-3;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2])*2;


u_x_max = 9e-3;
u_y_max = 9e-3;
u_z_max = 9e-3;
Rc = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2])*7e-1;

%% small difference
s_x_max = 3e-1;
s_y_max = 3e-1;
w_x_max = 6e-1;
w_y_max = 6e-1;
w_z_max = 4e-3;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2])*1.115;


u_x_max = 9e-3;
u_y_max = 9e-3;
u_z_max = 9e-3;
Rc = diag([1/u_x_max^2;1/u_y_max^2;1/u_z_max^2])*7e-1;

%% Again
% s_x_max = 3;
% s_y_max = 3;
% w_x_max = 10^-2;
% w_y_max = 10^-2;
% w_z_max = 10^-3;

s_x_max = 3e-1;
s_y_max = 3e-1;
w_x_max = 6e-1;
w_y_max = 6e-1;
w_z_max = 4e-3;
Qc = diag([1/s_x_max^2;1/s_y_max^2;1/w_x_max^2;1/w_y_max^2;1/w_z_max^2])*1.105;

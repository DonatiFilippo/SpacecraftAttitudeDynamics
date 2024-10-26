clc 
close all
clear all

%% Data
Iv = [0.07;0.055;0.025];

I = diag(Iv);
invI = inv(I);
w0 = [0.45; 0.5; 0.55];
A = [[1, 0, 0]; [0, 1, 0]; [0, 0, 1]];
%% Sim config
t0 = 0;
tf = 100;
%% Simulation Init
out = sim("Lab05.slx", "StartTime", "t0", "StopTime", "tf");

ma = zeros(1,4001);
for i = 1:4001
ma(1, i) = max(max(out.directCos(:,:,i)));
end

plot(ma)

%% Validation trash
m = ma(1, 1:3700);
[val, ind] = max(m);
out.eulerAngles(ind, :)
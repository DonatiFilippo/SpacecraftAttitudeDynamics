clc 
close all
clear all

%% Data
Iv = [0.07;0.055;0.025];

I = diag(Iv);
invI = inv(I);
w0 = [0.45; 0.52; 0.55];
A = [[1, 0, 0]; [0, sqrt(2)/2, sqrt(2)/2]; [0, -sqrt(2)/2, sqrt(2)/2]];
IC.angles = [0, 0.0000001, 0]; %The system is smart and choose Euler rappresentation based on singularities
%% Sim config
t0 = 0;
tf = 10000;
%% Simulation Init
out = sim("Lab05.slx", "StartTime", "t0", "StopTime", "tf");

%% Plot
% ma = zeros(1,140001);
% for i = 1:140001
% ma(1, i) = max(max(out.directCos(:,:,i)));
% end
val = abs(out.error);
mea = mean(val, [1,2]);
me(:) = mea(1,1,:);
figure
hold on
grid on
plot(me);
plot(out.flag1*1e-11);

figure
hold on
plot(out.eu313(1));
plot(out.eu312(1));
%% Validation trash

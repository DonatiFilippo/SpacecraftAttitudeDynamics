function [D] = DCM(angles)
%DCM Summary of this function goes here
%   we use 312
phi = angles(1);
theta = angles(2);
psi = angles(3);
D = [[cos(psi)*cos(phi) - sin(psi)*sin(phi)*sin(theta), cos(psi)*sin(phi) + sin(psi)*cos(phi)*sin(theta), -sin(psi)*cos(theta)];
     [-sin(phi)*cos(theta), cos(phi)*cos(theta), sin(theta)];
     [sin(psi)*cos(phi)+ cos(psi)*sin(phi)*sin(theta), sin(psi)*sin(phi) - cos(psi)*cos(phi)*sin(theta), cos(theta)*cos(psi)]];
end


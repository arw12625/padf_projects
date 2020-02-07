%%
clear all
close all

x = sdpvar(2,1);
y = sdpvar(1,1);

z = sdpvar(2,1);
z(1) = x(1);
z(2) = y(1);

con = [x>=0, y>=0, x(1)+y(1)<= 1, x(2) <= 2, z >= 0];

poly = Polyhedron(con);

plot(poly)

% this conversion does not appear to always work

%%
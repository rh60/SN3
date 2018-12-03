clear;clc;close;
load('../data/poisson.mat');
tri=double(msh.tri(:,1:3));

hold off
trisurf(tri,msh.x,msh.y,U);


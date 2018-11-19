clear;clc;close;
load('../data/msh.mat');
triplot(double(msh.tri),msh.x,msh.y);
axis equal

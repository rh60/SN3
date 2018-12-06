clear;clc;close;
load('../data/poisson.mat');
tri=double(msh.tri(:,1:3));
x=msh.x;y=msh.y;
hold off
trisurf(tri,x,y,U);
fu=@(x,y) exp(x+y)+1/2;
u=fu(x,y);
figure
trisurf(tri,x,y,u-U);
[m,i]=max(u-U)
[x(i),y(i)]

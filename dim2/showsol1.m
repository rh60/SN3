clear;clc;close;
load('../data/poisson.mat');
tri=double(msh.tri(:,1:3));
x=msh.x;y=msh.y;
hold off
trisurf(tri,x,y,U);
fu=@(x,y) 16*(y.*(1-y).*x.*(1-x));
u=fu(x,y);
figure
trisurf(tri,x,y,u);
[m,i]=max(u-U)
[x(i),y(i)]

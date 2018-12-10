clear;clc;close;
load('../data/ex4.mat');
tri=double(msh.tri(:,1:3));
x=msh.x;y=msh.y;
hold off
trisurf(tri,x,y,U);
fu=@(x,y) x.^2+2*y.^3;
u=fu(x,y);
figure
trisurf(tri,x,y,u);
[m,i]=max(u-U)
[x(i),y(i)]

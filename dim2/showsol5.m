clear;clc;close;
load('../data/ex5.mat');
tri=double(msh.tri(:,1:3));
x=msh.x;y=msh.y;
hold off
trisurf(tri,x,y,U);
fu=@(x,y) (x.^2+y.^2).^(1/3).*sin(2*atan2(y,x)/3+pi/3);
u=fu(x,y);
figure
trisurf(tri,x,y,u-U);
[m,i]=max(u-U)
[x(i),y(i)]

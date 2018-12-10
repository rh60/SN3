clear;clc;close;
[x, y, tri, U]=Chap8_CalculateExampleFem(5,5);
u=x.^2+2*y.^3;
trisurf(tri,x,y,U);
figure;
er=u-U;
trisurf(tri,x,y,er);
max(abs(er))
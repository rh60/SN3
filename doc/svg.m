clear;clc;close;
x=-1:0.5:1;
[X,Y]=meshgrid(x,x);
x=X(:);y=Y(:);
z=zeros(size(x));
z(13)=1;
tri=delaunay(x,y);

figure('Position',[0 0 200 100]);
trimesh(tri,x,y,z,'EdgeColor','k');
axis equal; axis off;

plot2svg
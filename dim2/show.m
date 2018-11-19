clear;clc;close;
load('../data/msh.mat');
triplot(double(msh.tri),msh.x,msh.y);
for i=1:numel(msh.x)
    text(msh.x(i),msh.y(i),int2str(i));
end
for i=1:size(msh.tri)
    x=mean(msh.x(msh.tri(i,:)));
    y=mean(msh.y(msh.tri(i,:)));
    text(x,y,int2str(i));
end

hold on
sides = {msh.yleft,msh.yright,msh.xleft,msh.xright};
colors = {'red','green','blue','black'};
for i=1:4
    s=sides{i};
    plot(msh.x(s),msh.y(s),'o','MarkerFaceColor',colors{i});
end    
axis equal
axis off
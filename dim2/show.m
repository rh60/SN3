clear;clc;close;
load('../data/msh.mat');
tri=double(msh.tri(:,1:3));
triplot(tri,msh.x,msh.y);


% for i=1:length(msh.tri)
%     x=mean(msh.x(tri(i,1:3)));
%     y=mean(msh.y(tri(i,1:3)));
%     text(x,y,int2str(i));
% end

t=msh.tri;
for it=t
    text(msh.x(it),msh.y(it),int2str(it));
end    

hold on
colors = {'red','green','blue','black'};
fn=fieldnames(msh.boundary);
for i=1:length(fn)
    s=str2double(fn{i});
    c=getfield(msh.boundary,fn{i});
    plot(msh.x(s),msh.y(s),'o','MarkerFaceColor',colors{c});
end
axis equal
axis off

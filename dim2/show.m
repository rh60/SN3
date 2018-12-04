clear;clc;close;
load('../data/poisson.mat');

showmsh(amsh)
figure
showmsh(msh)

function showmsh(msh)

tri=double(msh.tri(:,1:3));

hold off
triplot(tri,msh.x,msh.y);

for i=1:size(msh.tri,1)
    x=mean(msh.x(tri(i,1:3)));
    y=mean(msh.y(tri(i,1:3)));
    text(x,y,int2str(i),'FontSize',9);
end

for i=1:length(msh.x)    
    text(msh.x(i),msh.y(i),int2str(i),'FontSize',8);    
end    

axis equal
axis off

hold on
colors = {'red','green','blue','black'};

for i=1:length(msh.boundary)
   b=msh.boundary{i};
   ib=[b(:,1); b(end,2)];
   x=msh.x(ib);
   y=msh.y(ib);
   plot(x,y,'color',colors{i},'LineWidth',2) ;
end    

axis equal
axis off

return

for k=1:length(edg)
    ie=edg{k}(1:2);
    x=msh.x(ie);
    y=msh.y(ie);
    p=quiver(x(1),y(1),diff(x),diff(y),0,'k','LineWidth',3);             
    [edg{k}(3:5)';edg{k}(6:8)']
    pinfo(edg{k}(3:5),msh);
    pinfo(edg{k}(6:8),msh);
    delete(p);
end

function e=ij(i)
if i==3
    j=1
else
    j=i+1
end        
    e=[i,j];
end

function q=pinfo(ed,msh)
    if ed(1)==2
        e=ij(ed(3))
        ex=mean(msh.x(msh.tri(ed(2),e)));
        ey=mean(msh.y(msh.tri(ed(2),e)));
        q=text(ex,ey,int2str(ed(2:3)),'FontSize',14);
    else
        tmp=msh.boundary{ed(2)}(ed(3),:);
        ex=mean(msh.x(tmp));
        ey=mean(msh.y(tmp));
        q=text(ex,ey,int2str(ed(2:3)),'FontSize',14);
    end
    pause(2);
    delete(q);
end
end
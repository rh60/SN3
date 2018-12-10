clear;clc;close;
load('../data/ex4.mat');

a.b=zeros(size(a.b));
a.bvp.vbc{3}.value=@(x,y) 6*exp(x);
a.bvp.vbc{4}.value=@(x,y) 2*exp(x);

for i=3:4
    tmp=double(a.msh.boundary{i});
    n=size(tmp,1);
    for k=1:n
        a=CalculateContributionOnSegment(a,i,k);
    end
end

function t=affinity(P,x)
    t=P(1)*(1 - x)+P(2)*x;
end

function r=intSum(Q,values)
    r=dot(values,Q.Weights);
end

function a=CalculateContributionOnSegment(a,ib,k)
    l=length(a.Ls);
    tmp=double(a.msh.boundary{ib});
    loc2glob=tmp(k,:);
    t=loc2glob(1:2);
    x=a.msh.x(t);y=a.msh.y(t);
    h=norm([diff(x),diff(y)]);
    qx=affinity(x,a.Qs.Points);
    qy=affinity(y,a.Qs.Points);
    g=a.bvp.vbc{ib}.value(qx,qy);
    for i=1:l
        a.b(loc2glob(i)) = a.b(loc2glob(i)) + intSum(a.Qs,g.*a.Ls{i})*h;
    end
end
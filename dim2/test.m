P={[0, 0], [1, 0], [0, 1], [1/4, 0], [1/2, 0], [3/4, 0], [3/4, 1/4], [1/2, 1/2], [1/4, 3/4], [0, 3/4], [0, 1/2], [0, 1/4], [1/4, 1/4], [1/2, 1/4], [1/4, 1/2]};
n=length(P);
x=zeros(n,1);y=x;
for i=1:n
    x(i)=P{i}(1);
    y(i)=P{i}(2);
end
scatter(x,y);
for i=1:n
   text(x(i),y(i),num2str(i)); 
end
axis equal
axis off
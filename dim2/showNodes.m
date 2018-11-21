clear;clc;close;
x = [0.4, 0.6, 0.4, 0.45, 0.5, 0.55, 0.55, 0.5, 0.45, 0.4, 0.4, 0.4, 0.45, 0.5, 0.45]
y = [0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.05, 0.1, 0.15, 0.15, 0.1, 0.05, 0.05, 0.05, 0.1]
scatter(x,y);
for i=1:numel(x)
    text(x(i),y(i),int2str(i));
end
axis equal
axis off

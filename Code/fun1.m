function F = fun1(q,x,y,i);
F = [0.2*cos(q(1))+0.3*cos(q(1)+q(2))+0.1-x(i) , 0.2*sin(q(1))-0.3*sin(q(1)+q(2))-y(i)];
end
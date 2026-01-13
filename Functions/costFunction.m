function cost = costFunction(x,Ta,G,R,C,t,x0,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

y1=ModelFunction(x,Ta,G,R,C,t,x0)';
cost = rmse(y,y1);

end
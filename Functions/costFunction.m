function cost = costFunction(x, Ta, G, R, C, t, x0, y)
% costFunction

% Force vectors to be column vectors
x  = x(:);
Ta = Ta(:);
t  = t(:);
y  = y(:);

y1 = ModelFunction(x, Ta, G, R, C, t, x0);
y1 = y1(:);   % also ensure model output is a column vector

cost = rmse(y, y1);

end
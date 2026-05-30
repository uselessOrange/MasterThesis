function [y1,params]=GATrain(x, Ta, t, x0, y,lb,ub)

if nargin < 7, ub = [];
    if nargin < 6, lb = [];
    end
end

% force x, Ta, t, y to be row vectors
x = reshape(x,1,[]);
Ta = reshape(Ta,1,[]);
t = reshape(t,1,[]);
y = reshape(y,1,[]);

params = ga(@LocalcostFunc,3,[],[],[],[],lb,ub);
y1=ModelFunction(-x',Ta',params(1),params(2),params(3),t,x0);

function cost = LocalcostFunc(params)

cost = costFunction(-x',Ta',params(1),params(2),params(3),t,x0,y);

end
end
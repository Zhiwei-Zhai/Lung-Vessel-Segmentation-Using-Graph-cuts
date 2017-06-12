function y  = Sigmoid_f(x, alpha, beta)
%y  = Sigmoid_f(x, alpha, beta)
%   input: x, alpha, beta
%   output: y = 1/( 1 + exp( -alpha*(x - beta) ) )
%---------------------example
%   x = 1:0.1:10;
%   y = Sigmoid_f(x, 2, 5);
%   plot(x,y);

tmp = alpha*(x - beta);
y = 1./(1+ exp(-tmp));
end

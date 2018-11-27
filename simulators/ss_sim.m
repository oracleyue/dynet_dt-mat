function output = ss_sim(sys, data, method)

if nargin == 2
    method = 'Euler-Maruyama';
end

% unfold arguments
A = sys.A; B = sys.B; C = sys.C; D = sys.D; K = sys.K;
u = data.u; Ts = data.Ts; sigma = data.sigma; X0 = data.X0;
N = data.N; % the number of samples 

% calculate system dimensions
n = size(A,1); p = size(C,1); m = size(B,2);


% solve SDE
switch method
    case 'Euler-Maruyama'
        interval = 100;
        U = kron(u, ones(interval,1));
        dt = Ts/interval;
        
        X(1,:) = X0';
        for t = 1:1:N*interval
            dX = A*X(t,:)' + B*U(t,:)';                
            X(t+1,:) = X(t,:) + dt*dX' + ...
                sqrt(dt)*sigma*abs(X(t,:)).*randn(1,n);
            Y(t,:) = (C*X(t,:)')';
        end
        output = Y(1*interval:interval:N*interval, :);
    otherwise
        Error('Not valid method!');
end
    
% % plotting
% figure;
% plot(Ts*(1:N)', output, '-')

end 


% Local Functions
function [xout,zout] = simple_nonlinear_ring(xin,zin,input,a1d,b1d,Hidden)
d = length(xin);
Sign=ones(1,d);
%deterministic part
for i=2:d
    xout(i) = -xin(i) ...
        + Sign(i)*(1-Hidden)*xin(i-1)...
        + Sign(i)*Hidden*zin(i-1);
end
xout(1) = -xin(1)...
    + Sign(i)*(1-Hidden)*a1d*xin(d) ...
    + Sign(i)*Hidden*b1d*zin(d);
for i=1:d
    zout(i) = -zin(i) ...
        + Sign(i)*Hidden*xin(i);
end
xout(d)= xout(d)+(1-Hidden)*input;
zout(d)= zout(d)+Hidden*input;
end
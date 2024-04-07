clear; close all; clc;

% Parameters
%%%%%%%%%%%%
r=0.08; % interest rate
sigma=0.25; % volatility of underlying asset
Nt=1600; % number of time steps
Ns = 160; % number of asset price steps
Smax = 100; % maximum asset price considered
Smin = 0; % minimum asset price considered
T=1.; % Time until matrutity of the contract (1 year)
K = 35; % strike price
%%%%%%%%%%%%

dt = (T/Nt); % time steps
ds =(Smax-Smin)/Ns; %price step

% initializing the matrix of the option value
f(1:Ns+1,1:Nt+1) = 0.0;


% discretization of asset (S) and t variables
S = Smin+(0:Ns)*ds;
t=(0:Nt)*dt;

% terminal conditions (at expiration) f(S, t=0)
last_S_price = S(end);
f(1, Ns+1,1)=max(last_S_price-K,0);

% Boundary conditions 
f(1, 1:Nt+1)=0;     % f(0,t)=0 % contract worthless if stock price is 0
f(Ns+1,1:Nt+1)=Smax-K*exp(-r*t);    % f(S, t) = S-K*e^(-r*(T-t)) as S -> inf


for j=1:Nt  % loop over time
    for n=2:Ns  % loop over the stock price
        f(n, j+1)=0.5 * dt * (sigma*sigma*n*n-r*n) * f(n-1,j) +(1-dt*(sigma*sigma*n*n+r))*f(n,j)+0.5*dt*(sigma*sigma*n*n+r*n)*f(n+1,j);
    end
end


add_line = S-K;
figure(1)
plot(S, f(:,1), 'r-', S, f(:,round(Nt/2)), 'g-', S, f(:,Nt+1),'b-', S, add_line);
xlabel('S');
ylabel('Option Value');
legend('t=0', 't=t^*/2', 't=t^*');

figure(2)
mesh(t, S, f);
xlabel('t');
ylabel('Stock Price');
zlabel('Option Value');
view(45, 30);

clear; close all; clc;

function f = CallOption(K, r, T, sigma, Smax, Smin, Nt, Ns)
  % initializing the matrix of the option value
  f = zeros(Ns+1, Nt+1);


  % discretization of asset (S) and t variables
  dt=(T/Nt);
  ds =(Smax-Smin)/Ns;
  S = Smin+(0:Ns)*ds;
  t = (0:Nt)*dt;

  veti = 0:Ns;
  vetj = 0:Nt;

  % terminal conditions (at expiration) f(S, t*)
  f(:, Nt+1) = max(S-K,0);

  % Boundary conditions
  f(1,:) = 0;     % f(0,t)=0 % contract worthless if stock price is 0
  %disp(-r*dt*(Nt - vetj));
  f(Ns+1,:) = Smax - K*exp(-r*dt*(Nt - vetj));   % f(S, t) = S-K*e^(-r*(T-t)) as S -> inf

  a = 0.5*dt*(sigma^2*veti - r).*veti;
  b = 1 - dt*(sigma^2*veti.^2 + r);
  c = 0.5*dt*(sigma^2*veti+r).*veti;

  for j=Nt:-1:1  % loop over time
    for i=2:Ns  % loop over the stock price
        f(i,j)= a(i)*f(i-1,j+1) + b(i)*f(i,j+1)+ c(i)*f(i+1,j+1);
    end
  end


  t_0 = f(:,1);
  t_star = f(:, size(f, 2));
  t_half = f(:, round(Nt/2));
  figure(1)
  plot(S, t_0, 'r-', S, t_half, 'g-', S, t_star,'b-');
  xlabel('S');
  ylabel('Option Value');
  legend('t=0', 't=t^*/2', 't=t^*', 'Location', 'northwest');

  figure(2)
  mesh(t, S, f);
  xlabel('t');
  ylabel('Stock Price');
  zlabel('Option Value');
  view(45, 30);

end

function f = PutOption(K, r, T, sigma, Smax, Smin, Nt, Ns)
  % initializing the matrix of the option value
  f = zeros(Ns+1, Nt+1);


  % discretization of asset (S) and t variables
  dt=(T/Nt);
  ds =(Smax-Smin)/Ns;
  S = Smin+(0:Ns)*ds;
  t = (0:Nt)*dt;

  veti = 0:Ns;
  vetj = 0:Nt;

  % terminal conditions (at expiration) f(S, t*)
  f(:, Nt+1)=max(K-S,0);

  % Boundary conditions
  f(1,:)=K*exp(-r*dt*(Nt - vetj));    % f(S, t) = K*e^(-r*(T-t)) if S=0
  f(Ns+1,:)=0;    % f(0,t)=0 % contract worthless if S -> inf

  a = 0.5*dt*(sigma^2*veti - r).*veti;
  b = 1 - dt*(sigma^2*veti.^2 + r);
  c = 0.5*dt*(sigma^2*veti+r).*veti;

  for j=Nt:-1:1  % loop over time
    for i=2:Ns  % loop over the stock price
        f(i,j)= a(i)*f(i-1,j+1) + b(i)*f(i,j+1)+ c(i)*f(i+1,j+1);
    end
  end

  t_0 = f(:,1);
  t_star = f(:, size(f, 2));
  t_half = f(:, round(Nt/2));
  figure(3)
  plot(S, t_0, 'r-', S, t_half, 'g-', S, t_star,'b-');
  xlabel('S');
  ylabel('Option Value');
  legend('t=0', 't=t^*/2', 't=t^*', 'Location', 'northeast');

  figure(4)
  mesh(t, S, f);
  xlabel('t');
  ylabel('Stock Price');
  zlabel('Option Value');
  view(45, 30);
end



function delta = delta(f)
    h = 1 / size(f, 1);
    Nt = size(f, 2);
    Ns = size(f, 1);
    delta = zeros(Ns, Nt);
    for j = 1:Nt
      for i = 2:Ns
        delta(i, j) = (f(i, j) - f(i-1, j)) / h;
      end
    end
    delta = delta / 100;
end

%% plots the delta df at three different time stamps. Set call=1 if it's a call option.
function plotDelta(delta, call)
    Ns = size(delta, 1);
    Nt = size(delta, 2);
    figure();
    S = 2:Ns;
    t_0 = delta(2:Ns, 1);
    t_star = delta(2:Ns, Nt);
    t_half = delta(2:Ns, round(Nt/2));
    plot(S, t_0, S, t_half, S, t_star);
    if call == 1
      axis([0 size(delta, 1) 0 1])
    else
      axis([0 size(delta, 1) -1 0])
    endif

    xlabel('S');
    ylabel('Delta');
    legend('t=0', 't=t^*/2', 't=t^*');
  end


function theta = theta(f)
    h = 1 / size(f, 2);
    Nt = size(f, 2);
    Ns = size(f, 1);
    theta = zeros(Ns, Nt);
    for j = 2:Nt
      for i = 1:Ns
        theta(i, j) = (f(i, j) - f(i, j-1)) / h;
      end
    end
    theta = theta / 100;
end

function plotTheta(theta, K, call)
    Ns = size(theta, 1);
    Nt = size(theta, 2);

    figure();
    t = 1:Nt;
    if call == 1
      s_otm = theta(20, :);
      s_itm = theta(Ns - 50, :);
      s_atm = theta(K, :);
    else
      s_otm = theta(Ns - 60, :);
      s_itm = theta(20, :);
      s_atm = theta(K, :);
    endif
    plot(t, s_otm, t, s_atm, t, s_itm);
    xlabel('t');
    ylabel('Theta');
    legend('OTM', 'ATM', 'ITM');
end

% Parameters
K = 35;
r = 0.08;
T = 1;
sigma = 0.25;
Smax = 100;
Smin = 0;
Nt = 1600;
Ns = 100;

f_call = CallOption(K, r, T, sigma, Smax, Smin, Nt, Ns);
f_put = PutOption(K, r, T, sigma, Smax, Smin, Nt, Ns);

delta_call = delta(f_call);
plotDelta(delta_call, 1)
delta_put = delta(f_put);
plotDelta(delta_put, 0);

theta_call = theta(f_call);
theta_put = theta(f_put);
plotTheta(theta_call, K, 1);
plotTheta(theta_put, K, 0);

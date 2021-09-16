%% Risk of Outbreak vs Sigma
%
%  Calculate the probability of a vaccine escape outbreak occuring by 
%  conditioning on the time that the escape host emerges.
%
%%
clear; clc; close all;

% parameters for NOR function
Thresh = 50;    % threshold for a vaccine escape "outbreak"
t_run = 500;    % time to run integration over


N = 66.65e6;    % population size
tfin = 2000;    % time to run epidemic model over

mu = 0.1;               % recovery rate
beta = 2*mu/N;          % infection rate  %R number of 2%
eps = 0.2;              % ratio of infection rate of vaccinated and unvaccinated 
par = [beta, mu, eps];  % feed parameters into function

p0 = 2e-9;              % probability of vaccine escape in an unvaccinated individual
p1 = 2e-7;              % probability of vaccine escape in a vaccinated individual

m = 101;                    % no of splits for sigma
sigma = linspace(0,1,m);    % array of values for sigma from 0 to 1

% initialise array for risk plot
RiskPlotCOR = zeros(m,1);   
RiskPlotNOR = zeros(m,1);   
RiskPlotIOR = zeros(m,1);

R0counter = zeros(m,1);

for i=1:m
    sig = sigma(i);
    
    % Run Epidemic Model %
    y0 = [(1-sig)*N 100 0 sig*N 0 0];                               % inital conditions
    [t1,y] = ode45(@VIR_ode, [0, tfin], y0, [], par);               % use ode45 to perform integration
    
    S = y(:,1) + y(:,4);                                            % calculate no of susceptibles over time
    emergprob = ones(length(t1)) - (1-p0).^y(:,3).*(1-p1).^y(:,6);  % calculate cumulative probability of emergence before time t
    
    pdf = -mu.*(1-p0).^y(:,3).*(1-p1).^y(:,6).*(log(1-p0).*y(:,2) + log(1-p1).*y(:,5)); % derivative of the cdf

    
    %%% Risk Profile with IOR metric %%%
    IOR = max(1 - mu./(beta.*S), 0);            % calculate IOR metric for given sigma
    IORintegrand = IOR.*pdf;                    % multiply by pdf
    RiskPlotIOR(i) = trapz(t1, IORintegrand);   % approx integrate using trapezium rule

    
    %%% Risk Profile with NOR metric %%%
    tNOR = linspace(0, tfin, 100);              % calculate NOR metric for given sigma
    pdf2 = interp1(t1,pdf,tNOR);
    NOR = zeros( 1, length(tNOR));
    for j = 1:length(tNOR)                      
        NOR(j) = 1 - NOR_func(tNOR(j),Thresh,t_run,beta,S,t1,mu);
    end
    NORintegrand = NOR.*pdf2;                   % multiply by pdf
    RiskPlotNOR(i) = trapz(tNOR, NORintegrand); % approx integrate using trapezium rule

     
    %%% Risk Profile with COR metric %%%
    
    % R0 = beta*S/mu = 2*S/N
    % since S(0) = N, R0 is always initially greater than 1
    % we want to find the time that R0 < 1, ie S < N/2
    
    index = find(S <= N/2);             % find index of time that R0>1 first
    if isempty(index) == 1              % if R0 > 1 for all time
        IC = mu./(beta.*S(end));        % set p(tfin) = 1/R0(tfin) to match with IOR for all time after S(t) remains constant
        tspan = linspace(tfin, 0);      % solve backwards in time
        opts = odeset('RelTol',1e-15,'AbsTol',1e-16);
        [t2,p] = ode15s(@COR, tspan, IC, [], t1, S, par);
        t2 = flip(t2); p = flip(p);     % flip values to vectors are increasing, so that trapezium rule gives positive integral
        
    else                                % else if R0 < 1 for some t
        R0counter(i) = 1;               % track when R0 goes below 1 
        ind = find(S<=N*101/200, 1);    % index of time that R0 ~= 1.01
        tcrit  = t1(ind);               % match conditions at time before R0 hits 0
        
        IC = NOR_func(tcrit,Thresh,t_run,beta,S,t1,mu);             % match with NOR at some time just before R0 hits one, when NOR>0
        tspan = linspace(tcrit, 0);                                 % solve backwards from the "critial time"
        opts = odeset('RelTol',1e-13,'AbsTol',1e-16);
        [tau1,rho1] = ode15s(@COR, tspan, IC, opts, t1, S, par);
        
        IC = rho1(1);                                               % match with backwards solution
        tspan = linspace(tcrit, tfin, 2000);                        % solve forwards from the "critical time"
        [tau2,rho2] = ode15s(@COR, tspan, IC, opts, t1, S, par);
        
        tau1(1) = []; rho1(1) = [];                                 % make sure solutions don't double up at matching point
        t2 = [flip(tau1); tau2];                                    % combine time spans
        p = [flip(rho1); rho2];                                     % combine solutions for p
        
    end
    
    pdf3 = interp1(t1,pdf,t2);                  % interpolate pdf with timespan for the COR
    CORintegrand = (1-p).*pdf3;                 % combine pdf with conditional dist to form integrand
    RiskPlotCOR(i) = trapz(t2, CORintegrand);   % approx integrate using trapezium rule

end

critInd = find(R0counter == 1, 1, 'last');   % last time that R0 goes below 1


% plot probability of outbreak vs sigma 
figure(1)
plot(sigma, RiskPlotIOR, 'LineWidth', 1)
hold on
plot(sigma, RiskPlotCOR, 'LineWidth', 1)
hold on
plot(sigma, RiskPlotNOR, 'LineWidth', 1)
title('Probability of Outbreak vs Time')
xlabel('Proportion of Population Vaccinated \sigma')
ylabel('Probability')
legend('IOR', 'COR', 'NOR')
grid on

% plot(sigma(critInd), RiskPlot(critInd), 'r*')

%%
% SIR differential equations with y = (V0,I0,R0,V1,I1,R1)
% par = [beta, mu, eps]
function dydt = VIR_ode(t,y,par)
dydt = zeros(6,1);

dydt(1) = -par(1)*( y(2) + y(5))* y(1);
dydt(2) =  par(1)*(y(2)+y(5))*y(1) - par(2)*y(2);
dydt(3) =  par(2)*y(2);
dydt(4) = -par(3)*par(1)*( y(2) + y(5))* y(4);
dydt(5) =  par(3)*par(1)*( y(2) + y(5))* y(4) - par(2)*y(5);
dydt(6) =  par(2)*y(5);
end

%%
% Case Outbreak Risk COR as derived in notes
function dpdt = COR(t,p,t1,S,par)

S = interp1(t1,S,t);
% dpdt = par(1)*S*p*(1-p) + par(2)*(p-1);
dpdt = par(2)*( max( par(1)*S/par(2), 1)*p*(1-p) + (p-1) );

end 

%%
% (solve forward equations) 
% define the vector of length M=Thresh, p(t) = (p0(t), ... , pM(t)) where
% pi(t) is the probability that after t_run time, there are i infected. The
% threshold M is absorbing, and so pM(t) is the probability that at one
% time before t_run there were M infected at the same time.

function [ext,p,t] = NOR_func(t_start,Thresh,t_run,beta,S,t1,mu)

p0=zeros(Thresh+1,1);                       % vector of length Thresh+1
p0(2,1)=1;                                  % Start with one individual initially infected
tspan = linspace(t_start,t_start+t_run);    % run simulation for t_run time

[t,p] = ode45(@(t,Y) my_ode(t,Y,beta,S,t1,mu,Thresh),tspan,p0);     % Solve forward kolmogorov equations
ext = min(p(size(p,1),1),1);                                        % Determine if extinction occured, minimal non-negative solution 

function dYdt = my_ode(t,Y,beta,S,t1,mu,Thresh)

% Forward kolmogorov equations
dYdt = T_func(t,beta,mu,S,t1,Thresh)*Y;

%%%% Build generator matrix Q
function T = T_func(t,beta,mu,S,t1,N)
T=zeros(N+1,N+1);
v=linspace(0,N,N+1);

%intialise S
if t < t1(end)
   S = interp1(t1,S,t)*v;
else
   S = S(end)*v;
end

dt=mu*v;
dt(N+1) = 0;    % remove recoveries for N so it is an absorbing state

  for i=2:N                         % Define the transition matrix
      T(i,i)=-beta*S(i)-dt(i);      % diagonal entries
      T(i,i+1)=dt(i+1);             % superdiagonal entries
      T(i+1,i)=beta*S(i);           % subdiagonal entries
  end
T(1,2)=dt(2);
T(1,1) = 0;  
end
end

end
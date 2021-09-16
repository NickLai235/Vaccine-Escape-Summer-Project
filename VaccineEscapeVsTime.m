%% Probability of Vaccine Escape by time t (Deterministic)
%  
%  We can run the original deterministic model to calculate the changing
%  susceptible population over time. Then depending on when vaccine escape
%  occurs, the susceptible population at that time will determine the
%  likelihood of outbreak occuring.
%
%  We calculate the cumulative probability of vaccine escape emerging
%  before time t, and the corresponding pdf. We can then calculate the
%  probability of an outbreak occuring given the escape occurs at time t.
%  We can combine these to work out a total probability of an outbreak.
%
%%
clear; clc; close all;

% parameters for NOR function
Thresh = 50;    % threshold for a vaccine escape "outbreak"
t_run = 500;    % time to run integration over

N = 66.65e6;    % population size
tfin = 1000;    % time to run model over

mu = 0.1;               % recovery rate
beta = 2*mu/N;          % infection rate  %R number of 2%
eps = 0.2;              % ratio of infection rate of vaccinated and unvaccinated 
par = [beta, mu, eps];  % feed parameters into function

p0 = 2e-9;              % probability of vaccine escape in an unvaccinated individual
p1 = 2e-7;              % probability of vaccine escape in a vaccinated individual

%%%%%%%%%%%%%%%%%%%%%%%%
                    
sigma = [0.2 0.3 0.5];  % proportion of population initially vaccinated
                    
%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:3
    
    sig = sigma(k);
    y0 = [(1-sig)*N 100 0 sig*N 0 0];                   % inital conditions
    [t1,y] = ode45(@VIR_ode, [0, tfin], y0, [], par);   % use ode45 to perform integration
    
    S = y(:,1) + y(:,4);                                            % calculate no of susceptibles over time
    emergprob = ones(length(t1)) - (1-p0).^y(:,3).*(1-p1).^y(:,6);  % calculate cumulative probability of emergence before time t
    
    pdf = -mu.*(1-p0).^y(:,3).*(1-p1).^y(:,6).*(log(1-p0).*y(:,2) + log(1-p1).*y(:,5)); % derivative of the cdf
    
    % plot emergence probability before time t vs t %
    figure(1)
    plot(t1, emergprob, 'LineWidth', 1);
    title('Probability escape variant has emerged by time t')
    xlabel('time t')
    ylabel('Probability')
    grid on
    hold on 
    
    % plot proportion of population susceptible %
    figure(2)
    subplot(1,2,1)
    plot(t1, S./N, 'LineWidth', 1);
    title('Susceptible Population vs time')
    xlabel('time t')
    ylabel('Proportion of Total Population')
    grid on
    hold on
    
    % plot R0 vs time %
    subplot(1,2,2)
    plot(t1, (beta.*S)./mu, 'LineWidth', 1);
    title('R_0 = \beta S(t) / \mu')
    xlabel('time t')
    ylabel('R_0')
    grid on
    hold on
    
    %%
    % plot 3 different metrics to measure the probability an outbreak occurs
    % given vaccine escape appears at time t
    
    
    % parameters for NOR function %
    Thresh = 50;                    % threshold for a vaccine escape "outbreak"
    t_run = 500;                    % time to run integration over
    
    %%% run NOR simulation %%%
    tNOR = linspace(0, tfin, 100);
    NORplot = zeros(length(tNOR), 1);
    for i = 1:length(tNOR)
        NORplot(i) = 1 - NOR_func(tNOR(i),Thresh,t_run,beta,S,t1,mu);
    end
    
    %%% Run COR simulation %%%
    % R0 = beta*S/mu = 2*S/N
    % since S(0) = N, R0 is always initially greater than 1
    % we want to find the time that R0 < 1, ie S < N/2
    
    index = find(S <= N/2);             % find index of time that R0>1 first
    if isempty(index) == 1              % if R0 > 1 for all time
        IC = mu./(beta.*S(end));        % set p(tfin) = 1/R0(tfin) to match with IOR for all time after S(t) remains constant
        tspan = linspace(tfin, 0);      % solve backwards in time
        opts = odeset('RelTol',1e-15,'AbsTol',1e-16);
        [t2,p] = ode15s(@COR, tspan, IC, [], t1, S, par);
        t2 = flip(t2); p = flip(p);
        
    else                                % else if R0 < 1 for some t
        
        ind = find(S<=N*101/200, 1);    % index of time that R0 ~= 1.01
        tcrit  = t1(ind);               % match conditions at time before R0 hits 0
        
        IC = NOR_func(tcrit,Thresh,t_run,beta,S,t1,mu);             % match with NOR at some time just before R0 hits one, when NOR>0
        tspan = linspace(tcrit, 0);                                 % solve backwards from the "critial time"
        opts = odeset('RelTol',1e-13,'AbsTol',1e-16);
        [tau1,rho1] = ode15s(@COR, tspan, IC, opts, t1, S, par);
        
        IC = rho1(1);                                               % match with backwards solution
        tspan = linspace(tcrit, tfin);                              % solve forwards from the "critical time"
        [tau2,rho2] = ode15s(@COR, tspan, IC, opts, t1, S, par);
        
        tau1(1) = []; rho1(1) = [];                                 % make sure solutions don't double up at matching point
        t2 = [flip(tau1); tau2];                                    % combine time spans
        p = [flip(rho1); rho2];                                     % combine solutions for p
        
    end

    
    % plot probability of escape variant outbreak given escape appears at time t
    figure(3)
    subplot(3,1,k)
    plot( t2, 1 - p, 'LineWidth', 1)                        %plot COR
    hold on
    plot( t1, max(1 - mu./(beta.*S), 0), 'LineWidth', 1 )   %plot IOR
    hold on
    plot( tNOR, NORplot, 'LineWidth', 1)                    %plot NOR
    title( sprintf('Prob(escape variant outbreak | T_{VE} = t) vs time, \\sigma = %.2g', sig) )
    xlabel('time t')
    ylabel('Probability')
    legend('COR', 'IOR', 'NOR')
    grid on 

end 

% label values of sigma %
figure(1)
legend( sprintf('\\sigma = %.2g', sigma(1)), sprintf('\\sigma = %.2g', sigma(2)), sprintf('\\sigma = %.2g', sigma(3)), 'Location', 'northwest')
figure(2)
subplot(1,2,1)
legend( sprintf('\\sigma = %.2g', sigma(1)), sprintf('\\sigma = %.2g', sigma(2)), sprintf('\\sigma = %.2g', sigma(3)))
subplot(1,2,2)
legend( sprintf('\\sigma = %.2g', sigma(1)), sprintf('\\sigma = %.2g', sigma(2)), sprintf('\\sigma = %.2g', sigma(3)))


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

% intialise S
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



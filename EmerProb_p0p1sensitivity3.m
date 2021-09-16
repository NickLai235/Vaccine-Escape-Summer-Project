%% Emergence Probability - Sensitivity to Choices of p0 and p1
%  exploring the sensitivity of the emergence probability in changes in p0
%  and p1, the probabilities of vaccine escape in unvaccinated and
%  vaccinated hosts resp.
%  for different values of p0 and p1 we compute the value of sigma that
%  attains the maximum emergence probability, and plot as a contour map.
%
%%
clear; clc; close all;
% let y be the vector y = (V0,I0,R0,V1,I1,R1) splitting the population into
% unvaccinated and susceptible, unvaccinated and infected, unvaccinated and
% removed, and analagously for vaccinated.

N = 66.65e6;        % population size
I0 = 100;           % one initial infection

mu = 0.1;           % recovery rate
beta = 2*mu/N;      % infection rate  %R number of 2%
eps = 0.2;          % ratio of infection rate of vaccinated and unvaccinated

n = 17;                     % no of splits for p0 and p1
xx = linspace(-10, -6, n);   % x and y axis, log10(p0) and log10(p1)
P = 10.^xx;                 % array of values for p0 and p1

m = 101;                    % no of splits for sigma
sig = linspace(0,1,m);      % array of values for sigma from 0 to 1, excluding boundary values 
emergprob = zeros(m,1);     % initialise array of emergence probabilities

 
argmaxEmerProb = zeros(n,n);% intialise array of argmax{P(emergence)} 
maxEmerProb = zeros(n,n);   % intialise array of max{P(emergence)}

for j = 1:n             % for n different values of p0
    for k = 1:n         % for n different values of p1
        
        for l = 1:m
            Rnext = R_inf(sig(l), mu, beta, N, I0, eps);                % calculate final size for each value of sigma (this is a vector)
            
            emergprob(l) = 1 - (1-P(k))^Rnext(2)*(1-P(j))^Rnext(1);     % emergence probability for given sigma, as per notes, and given p0 and p1
        end
        
        [M, I] = max(emergprob);
        argmaxEmerProb(k,j) = sig(I); % compute argmax{Prob(emergence)}
        maxEmerProb(k,j) = M;         % compute max{Prob(emergence)}
    end 
end

% plot contour map of argmax{Prob(emergence)} with color bar
figure(1)
contourf( xx, xx, argmaxEmerProb);                
hold on
plot(log10(2e-9), log10(2e-7), 'r*')
xlabel('log_{10}(p0)');
ylabel('log_{10}(p1)');
title('Sensitivity of Model to p0 and p1');
cb = contourcbar;
cb.XLabel.String = 'Proportion of Population Vaccinated \sigma';

% plot contour map of max{Prob(emergence)} with color bar
figure(2)
contourf( xx, xx, maxEmerProb);                   
hold on
plot(log10(2e-9), log10(2e-7), 'r*')
xlabel('log_{10}(p0)');
ylabel('log_{10}(p1)');
title('Sensitivity of Model to p0 and p1');
cb = contourcbar;
cb.XLabel.String = 'Probability of Vaccine Escape';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function f = R_inf(sig, mu, beta, N, I0, eps) 
% given these parameters, output is a vector of the final size of the
% infected unvaccinated and vaccinated populations, (R0(inf),R1(inf)).

f = zeros(2,1);

if sig<1 && sig>=0  % if not all the population is vaccinated
    syms x         
    finalsizeeqn0 = x == (mu/beta)*log(x/((1-sig)*N)) - sig*N/((1-sig)*N)^eps * x^eps +I0 + N;   % final size equation for x=V0(inf) as derived in notes
    V0inf = vpasolve(finalsizeeqn0, x, 1);      % solve for V0(inf)
    V1inf = sig*N/((1-sig)*N)^eps * V0inf^eps;  % calculate V1(inf)

    R0inf = (1-sig)*N - V0inf + I0;             % since V0+I0+R0=(1-sig)N +1 and I0(inf)=0, taking limit t->inf, we can find R0(inf)
    R1inf = (sig)*N - V1inf;                    % since V1+I1+R1=(sig)N and I1(inf)=0, taking limit t->inf, we can find R1(inf)
    
elseif sig == 1     % if all the population is vaccinated
    syms x
    finalsizeeqn1 = x ==(mu/(beta*eps))*log(x/N) + N + I0;   % final size equation for x=V1(inf) as derived in notes
    V0inf = 0;                              % sigma = 0 -> V0 == 0 for all t 
    V1inf = vpasolve(finalsizeeqn1, x, 1);  % solve for V1(inf)
    
    R0inf = (1-sig)*N - V0inf + I0;         % since V0=0, I0+R0=1 and I0(inf)=0, -> R0(inf) = 1
    R1inf = (sig)*N - V1inf;                % since V1+I1+R1=N and I1(inf)=0, we can find R1(inf)
end

f(1) = R0inf;
f(2) = R1inf;
end




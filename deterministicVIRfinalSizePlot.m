%% Final Size for Vaccinated Hosts.
%  Solving a set of ODEs to describe a compartmentalised model of
%  vaccinated and unvaccinated populations for the final size of the
%  epidemic, and plotting vs sigma, the proportion of intially vaccinated
%  individuals

%%
clear; clc; close all; 

% let y be the vector y = (V0,I0,R0,V1,I1,R1) splitting the population into
% unvaccinated and susceptible, unvaccinated and infected, unvaccinated and
% removed, and analagously for vaccinated.

N = 66.65e6;         % population size
I0 = 100;            % 100 initial infections

mu = 0.1;            % recovery rate
beta = 2*mu/N;       % infection rate, setting R number to 2
eps = [0.1 0.2 0.3]; % ratio of infection rate of vaccinated and unvaccinated 
p0 = 2e-9;           % probability of vaccine escape in an unvaccinated individual
p1 = 2e-7;           % probability of vaccine escape in a vaccinated individual

m = 101;                    % no of splits in sigma
sig = linspace(0,1,m);      % array of values for sigma from 0 to 1


finalsizeplot = zeros(m,1); % initialise array of final sizes
propcases = zeros(m,1);     % initialise array of proportion of cases in vaccinated individuals
emergprob = zeros(m,1);     % initialise array of emergence probabilities

for i = 1:3
    Rnovax = R_inf(sig(1), mu, beta, N, I0, eps(i));             % calculate final size with no vaccine, sig(1) = 0
    finalsizenovax = Rnovax(1) + Rnovax(2);                      
    
    for j = 1:m
        Rnext = R_inf(sig(j), mu, beta, N, I0, eps(i));          % calculate final size for each value of sigma (this is a vector)
         
        finalsizeplot(j) = (Rnext(1) + Rnext(2))/finalsizenovax; % final size of infected population for given sigma, vaccinated and unvaccinated, relative to without vaccine
        propcases(j) = Rnext(2)/(Rnext(1)+Rnext(2));             % proportion of cases amongst vaccinated individuals
        emergprob(j) = 1 - (1-p1)^Rnext(2)*(1-p0)^Rnext(1);      % probability of the emergence of vaccine escape for given sigma, as per notes
    end

    % plot final size 
    figure(1)
    plot(sig(1:m), finalsizeplot)
    title('Final Size vs Sigma');
    xlabel('Proportion of Population Vaccinated \sigma');
    ylabel('Final size relative to in abscence of vaccine');
    grid on
    hold on   
    % plot proportion of cases in vaccinated hosts
    figure(2)
    plot(sig(1:m), propcases)
    title('Proportion of Cases in Vaccinated Hosts vs Sigma');
    xlabel('Proportion of Population Vaccinated \sigma');
    ylabel('Proportion of Cases in Vaccinated Individuals');
    grid on
    hold on
    % plot emergence probability
    figure(3)
    plot(sig(1:m), emergprob)
    title(sprintf('Probability of Vaccine Escape vs Sigma (p0 = %.3g, p1 = %.3g)', p0, p1));
    xlabel('Proportion of Population Vaccinated \sigma');
    ylabel('Probability of Emergence');
    grid on
    hold on
    
end

% label values of epilson to 3sf
figure(1)
legend( sprintf('\\epsilon = %.3g', eps(1)), sprintf('\\epsilon = %.3g', eps(2)), sprintf('\\epsilon = %.3g', eps(3)), 'Location', 'northeast')
figure(2)
legend( sprintf('\\epsilon = %.3g', eps(1)), sprintf('\\epsilon = %.3g', eps(2)), sprintf('\\epsilon = %.3g', eps(3)), 'Location', 'northwest')
figure(3)
legend( sprintf('\\epsilon = %.3g', eps(1)), sprintf('\\epsilon = %.3g', eps(2)), sprintf('\\epsilon = %.3g', eps(3)), 'Location', 'northwest')
hold off



%%
function f = R_inf(sig, mu, beta, N, I0, eps) 
% given these parameters output is a vector of the final size of the infected unvaccinated and vaccinated populations, (R0(inf),R1(inf))

f = zeros(2,1);

if sig<1 && sig>=0 % if not all the population is vaccinated
    syms x         % x = V0(inf)
    finalsizeeqn0 = x == (mu/beta)*log(x/((1-sig)*N)) - sig*N/((1-sig)*N)^eps * x^eps +I0 + N;   % final size equation for V0(inf) as derived in notes
    V0inf = vpasolve(finalsizeeqn0, x, 1);      % solve for V0(inf)
    V1inf = sig*N/((1-sig)*N)^eps * V0inf^eps;  % calculate V1(inf)

    f(1) = (1-sig)*N - V0inf + I0;              % since V0+I0+R0=(1-sig)N +1 and I0(inf)=0, we can find R0(inf)
    f(2) = (sig)*N - V1inf;                     % since V1+I1+R1=(sig)N and I1(inf)=0, we can find R1(inf)
    
elseif sig == 1 % if all the population is vaccinated
    syms x      % x = V1(inf)
    finalsizeeqn1 = x ==(mu/(beta*eps))*log(x/N) + N + I0;   % final size equation for V1(inf) as derived in notes
    V0inf = 0;                              % sigma = 0 -> V0 == 0 for all t 
    V1inf = vpasolve(finalsizeeqn1, x, 1);  % solve for V1(inf)
    
    f(1) = (1-sig)*N - V0inf + I0;          % since V0=0, I0+R0=1 and I0(inf)=0, -> R0(inf) = 1
    f(2) = (sig)*N - V1inf;                 % since V1+I1+R1=N and I1(inf)=0, we can find R1(inf)
end
end





%% Final Size for Vaccinated Hosts 2 Dose Model
%  Solving a set of ODEs to describe a compartmentalised model of
%  individuals vaccinated with 2, 1, and 0 doses for the final size of the
%  epidemic, and plotting vs theta, the proportion of vaccination
%  administered as two doses
%
%  fix M = N doses so there are enough for everyone to have one single dose
%
%%
clear; clc; close all; 

% let y be the vector y = (V0,I0,R0,V1,I1,R1,V2,I2,R2) splitting the population into
% susceptible, infected, and removed with 0,1,2 doses respectively 

N = 66.65e6;            % population size
M = N;                  % total number of doses available
I0 = 100;               % 100 initial infections

mu = 0.1;               % recovery rate
beta = 2*mu/N;          % infection rate, setting R number to 2
eps1 = [0.1, 0.2, 0.3]; % reduction in infection rate for 1 dose  
eps2 = 0.1;             % reduction in infection rate for 2 dose  

p0 = 2e-9;              % probability of vaccine escape in an unvaccinated individual
p1 = 2e-7;              % probability of vaccine escape in a single dose vaccinated individual
p2 = p0;                % probability of vaccine escape in a double dose vaccinated individual

m = 101;                                    % no of splits in sigma
theta = linspace( 0, min( 2*N/M-1, 1), m);  % array of values for theta, the proportion of vaccine given as a single dose


finalsizeplot = zeros(m,1);     % initialise array of final sizes
onedosepropcases = zeros(m,1);  % initialise array of proportion of cases in 1-dose individuals
twodosepropcases = zeros(m,1);  % initialise array of proportion of cases in 2-dose individuals

onedoseemergprob = zeros(m,1);  % initialise array of vaccine escape prob in 1-dose host
twodoseemergprob = zeros(m,1);  % initialise array of vaccine escape prob in 2-dose host
emergprob = zeros(m,1);         % initialise array of vaccine escape probabilities

for i = 1:3
%    Rnovax = R_inf(theta(1), mu, beta, N, I0, M, eps1, eps2);             % calculate final size with no vaccine, sig(1) = 0
%    finalsizenovax = Rnovax(1) + Rnovax(2);                      
    
    for j = 1:m
        Rnext = R_inf(theta(j), mu, beta, N, I0, M, eps1(i), eps2);             % calculate final size for each value of theta (this is a vector)
         
        finalsizeplot(j) = (Rnext(1) + Rnext(2)) + Rnext(3);                    % final size of infected population for given sigma, vaccinated and unvaccinated, relative to without vaccine
        onedosepropcases(j) = Rnext(2)/(Rnext(1)+Rnext(2)+Rnext(3));            % proportion of cases amongst single dose vaccinated individuals
        twodosepropcases(j) = Rnext(3)/(Rnext(1)+Rnext(2)+Rnext(3));            % proportion of cases amongst single dose vaccinated individuals
        
        onedoseemergprob(j) = 1 - (1-p1)^Rnext(2);                              % probability of vaccine escape in 1-dose host
        twodoseemergprob(j) = 1 - (1-p2)^Rnext(3);                              % probability of vaccine escape in 2-dose host
        emergprob(j) = 1 - (1-p2)^Rnext(3)*(1-p1)^Rnext(2)*(1-p0)^Rnext(1);     % probability of the emergence of vaccine escape for given theta, as per notes
    end

    % plot final size 
    figure(1)
    plot(theta(1:m), finalsizeplot, 'LineWidth', 1)
    title('Final Size vs \theta');
    xlabel('Proportion of vaccines given as single dose, \theta');
    ylabel('Final size relative to in abscence of vaccine');
    grid on
    hold on
    
    % plot proportion of cases in vaccinated hosts
    figure(2)
    subplot(1,2,1);
    plot(theta(1:m), onedosepropcases, 'LineWidth', 1)
    title('Proportion of Cases in 1-Dose Hosts vs \theta');
    xlabel('Proportion of vaccine given as single dose, \theta');
    ylabel('Proportion of cases in 1-dose hosts');
    grid on
    hold on
    subplot(1,2,2);
    plot(theta(1:m), twodosepropcases, 'LineWidth', 1)
    title('Proportion of Cases in 2-Dose Hosts vs \theta');
    xlabel('Proportion of vaccine given as single dose, \theta');
    ylabel('Proportion of cases in 2-dose hosts');
    grid on
    hold on
    
    % plot probability of vaccine escape in 1-dose hosts and 2-dose hosts separately
    figure(3)
    subplot(1,2,1);
    plot(theta(1:m), onedoseemergprob, 'LineWidth', 1)
    title(sprintf('P(Vaccine Escape in 1-dose host) vs \\theta (p_1 = %.3g)', p1));
    xlabel('Proportion of vaccine given as single dose, \theta');
    ylabel('Probability');
    grid on
    hold on
    subplot(1,2,2);
    plot(theta(1:m), twodoseemergprob, 'LineWidth', 1)
    title(sprintf('P(Vaccine Escape in 2-dose host) vs \\theta (p_2 = %.3g)', p2));
    xlabel('Proportion of vaccine given as single dose, \theta');
    ylabel('Probability');
    grid on
    hold on
    
    % plot total probability of vaccine escape
    figure(4)
    plot(theta(1:m), emergprob, 'LineWidth', 1)
    title(sprintf('P(Vaccine Escape) vs \\theta (p_0 = p_2 = %.3g, p_1 = %.3g)', p0, p1));
    xlabel('Proportion of vaccine given as single dose, \theta');
    ylabel('Probability');
    grid on
    hold on
    
end

% label values of epilson to 3sf
figure(1)
legend( sprintf('\\epsilon_1 = %.3g', eps1(1)), sprintf('\\epsilon_1 = %.3g', eps1(2)), sprintf('\\epsilon_1 = %.3g', eps1(3)), 'Location', 'northeast')
figure(2)
subplot(1,2,1)
legend( sprintf('\\epsilon_1 = %.3g', eps1(1)), sprintf('\\epsilon_1 = %.3g', eps1(2)), sprintf('\\epsilon_1 = %.3g', eps1(3)), 'Location', 'northwest')
figure(3)
subplot(1,2,1)
legend( sprintf('\\epsilon_1 = %.3g', eps1(1)), sprintf('\\epsilon_1 = %.3g', eps1(2)), sprintf('\\epsilon_1 = %.3g', eps1(3)), 'Location', 'northeast')
subplot(1,2,2)
legend( sprintf('\\epsilon_1 = %.3g', eps1(1)), sprintf('\\epsilon_1 = %.3g', eps1(2)), sprintf('\\epsilon_1 = %.3g', eps1(3)), 'Location', 'northeast')
figure(4)
legend( sprintf('\\epsilon_1 = %.3g', eps1(1)), sprintf('\\epsilon_1 = %.3g', eps1(2)), sprintf('\\epsilon_1 = %.3g', eps1(3)), 'Location', 'northeast')
hold off


%%
function f = R_inf(theta, mu, beta, N, I0, M, eps1, eps2) 
% given these parameters output is a vector of the final size of infected individuals vaccinated with 0,1,2 doses, (R0(inf),R1(inf),R2(inf))

f = zeros(3,1);
IC = [ N - (1+theta)*M/2 , theta*M , (1-theta)*M/2 ];   % set initial conditions

if IC(1)>0  % if there is a positive proportion of the population unvaccinated
    
    syms x  % x = V0(inf)
    
    finalsizeeqn0 = 0 == - x - IC(2)/(IC(1))^eps1 * x^eps1 - IC(3)/(IC(1))^eps2 * x^eps2 + (mu/beta)*log(x/IC(1)) + I0 + N ;   % final size equation for V0(inf) as derived in notes
    
    V0inf = vpasolve(finalsizeeqn0, x, 1);      % solve for V0(inf)
    V1inf = IC(2)/(IC(1))^eps1 * V0inf^eps1;    % calculate V1(inf)
    V2inf = IC(3)/(IC(1))^eps2 * V0inf^eps2;    % calculate V2(inf)
    
elseif IC(1) == 0   % if there are no unvaccinated individuals
    
    syms x          % x = V1(inf)
    
    finalsizeeqn1 = 0 == - x - IC(3)/(IC(2))^(eps2/eps1) * x^(eps2/eps1) + (mu/(beta*eps1))*log(x/IC(2)) + N + I0;      % final size equation for V1(inf) as derived in notes
    
    V0inf = 0;                                              % V0 == 0 for all t 
    V1inf = vpasolve(finalsizeeqn1, x, 1);                  % solve for V1(inf)
    V2inf = IC(3)/(IC(2))^(eps2/eps1) * V1inf^(eps2/eps1);  % calculate V2(inf)
 
end

f(1) = IC(1) - V0inf + I0;                  % since V0+I0+R0=(1-sig)N +1 and I0(inf)=0, we can find R0(inf)
f(2) = IC(2) - V1inf;                       % since V1+I1+R1=(sig)N and I1(inf)=0, we can find R1(inf)
f(3) = IC(3) - V2inf;                       % since V1+I1+R1=(sig)N and I1(inf)=0, we can find R1(inf)

end





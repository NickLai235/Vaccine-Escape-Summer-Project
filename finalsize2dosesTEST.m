%% Final Size Equation for Infected Vaccinated Population 1 Dose
%  solving the system of differential equations numerically using ode45 in
%  Matlab, and plotting the final size against sigma
%%
clear; clc; close all; 

% let y be the vector y = (V0,I0,R0,V1,I1,R1,V2,I2,R2) splitting the population into
% susceptible, infected, and removed with 0,1,2 doses respectively 

N = 66.65e6;            % population size
M = .5*N;                  % total number of doses available
I0 = 100;               % 100 initial infections
tfin = 1000;

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
        par = [beta, mu, eps1(i), eps2];      %feed parameters into function      
        y0 = [ N - (1+theta(j))*M/2 I0 0 theta(j)*M 0 0 (1-theta(j))*M/2 0 0 ];    %set initial conditions   
        
        [t,y] = ode45(@SIR_ode, [0, tfin], y0, [], par); %use ode45 to perform integration
        Rnext(1) = y(end,3);                          %total no. of infected unvaccinated hosts
        Rnext(2) = y(end,6);                          %total no. of infected vaccinated hosts
        Rnext(3) = y(end,9);
         
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
%SIR differential equations with y = (V0,I0,R0,V1,I1,R1)
%par = [beta, mu, eps1(i), eps2]
function dydt = SIR_ode(t,y,par)
dydt = zeros(9,1);

dydt(1) = -par(1)*( y(2) + y(5) + y(8) )* y(1);
dydt(2) =  par(1)*( y(2) + y(5) + y(8) )*y(1) - par(2)*y(2);
dydt(3) =  par(2)*y(2);
dydt(4) = -par(3)*par(1)*( y(2) + y(5) + y(8) )* y(4);
dydt(5) =  par(3)*par(1)*( y(2) + y(5) + y(8) )* y(4) - par(2)*y(5);
dydt(6) =  par(2)*y(5);
dydt(7) = -par(4)*par(1)*( y(2) + y(5) + y(8) )* y(7);
dydt(8) =  par(4)*par(1)*( y(2) + y(5) + y(8) )* y(7) - par(2)*y(8);
dydt(9) =  par(2)*y(8);
end

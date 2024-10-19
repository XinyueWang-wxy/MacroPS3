clear

%% 2.1 Calibration
% Defining the known parameters
n = 0.01*ones(1,5);
gamma = 0.015*ones(1,5);
rho = [zeros(1,3),-0.25,0.167];
sigma = [2,1,4,2,2];
paras = zeros(4,5);
% Calibration with the four steady state equations
for i = 1:5
paras(:,i) = parasolve(n(i),gamma(i),rho(i),sigma(i));
end
tb1 = table(["Benchmark","IES = 1.00","IES = 0.25", "ES = 0.8","ES = 1.2"]',n',gamma',rho',sigma',paras(1,:)',paras(2,:)',paras(3,:)',paras(4,:)');
tb1.Properties.VariableNames = {'var','n','gamma','rho','sigma','beta','delta', 'theta','B'};


%% Calculating starting k0_hat
k0_hat = zeros(5,1);
for i =1:5
k0_hat(i) = k0_solve(paras(3,i),paras(4,i),rho(i));
end
% Reporting
tb2 = table(["Benchmark","IES = 1.00","IES = 0.25", "ES = 0.8","ES = 1.2"]',n', ...
gamma',rho',sigma',k0_hat);
tb2.Properties.VariableNames = {'var','n','gamma','rho','sigma','k0_hat'};

%% Plot
% (1)
T = 100; %Setting T*
TOL = 10^-16; %Setting tolerance level
k = zeros(T+1,5);
y = zeros(T+1,5);
y_der = zeros(T+1,5);
c_bar = zeros(1,5);
for i =1:5
k(1,i) = k0_hat(i);
y(1,i) = fk(k(1,i),paras(3,i),rho(i),paras(4,i));
intval2 = (0.7*fk(k(1,i),paras(3,i),rho(i),paras(4,i))+(1-paras(2,i))* ...
k(1,i))/(1+n(i))*(1+gamma(i));
intval = [0,intval2];
c_bar(i) = 0.3*y(1,i);
while abs(intval(1)-intval(2)) > TOL
k1_tilde = sum(intval)/2;
k(2,i) = k1_tilde;
for t = 1:T-1
y(t+1,i) = fk(k(t+1,i),paras(3,i),rho(i),paras(4,i));
y_der(t+1,i) = paras(3,i)*paras(4,i)^rho(i)*(k(t+1,i)/ ...
y(t+1,i))^(rho(i)-1);
% Solving k_t+2 for the Euler Equation WRT t (7)
term1=real((y(t,i)+(1-paras(2,i))*k(t,i)-(1+n(i))*(1+gamma(i))* k(t+1,i)-c_bar(i))/...
((1+gamma(i))*(paras(1,i)*(y_der(t+1,i)+ (1-paras(2,i))))^(-1/sigma(i))));
term2=real(y(t+1,i)+(1-paras(2,i))*k(t+1,i)-c_bar(i)-term1);
k_next = real(term2/((1+n(i))*(1+gamma(i))));
if k_next<0
k(T+1,i) = 0;
break
%We set up a break system to prevent any trajectory that collapase to kt=0 before T
else
k(t+2,i) = k_next;
end
end

%(2)
y(T+1,i) = fk(k(T+1,i),paras(3,i),rho(i),paras(4,i));
if k(T+1,i)<3
intval(1) = k1_tilde;
else
intval(2) = k1_tilde;
end
end
end

capital_output = zeros(T+1,5);
for i =1:5
for t = 1:T+1
capital_output(t,i) = k(t,i)/y(t,i);
end
end
investment_rate = zeros(T,5);
for i =1:5
for t = 1:T
investment_rate(t,i) = ((1+gamma(i))*(1+n(i))*k(t+1,i)- ...
(1-paras(2,i))*k(t,i))/y(t,i);
end
end

mpk = zeros(T+1,5);
for i =1:5
for t = 1:T+1
mpk(t,i) = paras(3,i)*(paras(4,i))^(rho(i))*(k(t,i)/y(t,i))^(rho(i)-1);
end
end

gdp_pc = zeros(T+1,5);
for i =1:5
for t = 1:T+1
% We normalize A0 as 1 for simplicity, but in the end, this
% normalization doesn't matter since we are only using the log GDP
% per capita and the A0 term will end up as a constant.
gdp_pc(t,i) = y(t,i)*(1+gamma(i))^(t-1);
end
end

gdp_pcg = zeros(T,5);
for i =1:5
for t = 1:T
gdp_pcg(t,i) = (gdp_pc(t+1,i)-gdp_pc(t,i))/gdp_pc(t,i);
end
end


%% Plotting for each set of parameters
title_name = ["Benchmark","IES=1.00","IES=0.25","ES=0.80","E=1.20"];
for i = 1:5
figure
yyaxis left
plot(1:T+1, capital_output(:,i),'LineWidth',2.0);
xlabel("Time");
yyaxis right;
plot(1:T, investment_rate(:,i),'LineWidth',2.0);
hold on
plot(1:T+1, mpk(:,i),'LineWidth',2.0);
hold on
plot(1:T, gdp_pcg(:,i),'LineWidth',2.0);
hold on
title (title_name(i));
legend("Capital-output ratio","Investment rate","Marginal production of capital", ...
"Growth rate of GDP per capita",'Location',"best");
end

%% (c)

figure
for i = 1:5
plot(1:T+1, capital_output(:,i),'LineWidth',1.5);
hold on
end
ylabel("Capital Output ratio");
xlabel("Time");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure 1')


figure
for i = 1:5
plot(1:T, investment_rate(:,i),'LineWidth',1.5);
hold on
end
ylabel("Investment rate");
xlabel("Time");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure 2')


figure
for i = 1:5
plot(1:T, mpk(1:100,i),'LineWidth',1.5);
hold on
end
ylabel("Marginal production of capital");
xlabel("Time");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure 3')


figure
for i = 1:5
plot(1:T, gdp_pcg(:,i),'LineWidth',1.5);
hold on
end
ylabel("Growth rate of GDP per capita");
xlabel("Time");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure 4')

lggdp_pc = log(gdp_pc);
figure
for i = 1:5
plot(lggdp_pc(:,i), capital_output(:,i),'LineWidth',1.5);
hold on
end
ylabel("Capital Output Ratio");
xlabel("Log GDP per capita");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure 5')


figure
for i = 1:5
plot(lggdp_pc(1:100,i), investment_rate(:,i),'LineWidth',1.5)
hold on
end
ylabel("Investment Rate");
xlabel("Log GDP per capita");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure 6')


figure
for i = 1:5
plot(lggdp_pc(:,i), mpk(:,i),'LineWidth',1.5)
hold on
end
ylabel("Marginal Production of Capital");
xlabel("Log GDP per capita");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure7')

figure
for i = 1:5
plot(lggdp_pc(1:100,i), gdp_pcg(:,i),'LineWidth',1.5)
hold on
end
ylabel("GDP per capita Growth Rate");
xlabel("Log GDP per capita");
legend("Benchmark","IES=1.00","IES=0.25","ES=0.80","ES=1.20",'Location',"best");
title('Figure 8')



%% (d) Report T, yT/y0
% The growth of output per capita y_t/y_0
T90 = 44;
for i =1:5
grow(i,:) = (1+gamma(i))^T90*y(T90,i)/y(1,i);
end
% The average growth rate of output per capita
for i = 1:5
av_grow(i,:) = (grow(i,:)^(1/T90)-1);
end


% Decomposition of output growth
total_grow = log(grow);
for i = 1:5
tech_grow(i,:) = T90*log(1+gamma(i))/total_grow(i,:);
% tech_grow(i,:) = ((1+gamma(i))^T90)/grow(i,:);
end
tran_grow = 1-tech_grow;
tb3 = table(["Benchmark","IES = 1.00","IES = 0.25", "ES = 0.8","ES = 1.2"]', ...
grow,av_grow,tech_grow,tran_grow);
tb3.Properties.VariableNames = {' ','yt/y0','Average Growth Rate',...
'Growth from Tech','Growth from Transitions'}



%% Function
function S = parasolve(n, gamma, rho, sigma)
syms betta delta theta B
eqn1 = betta*(1+gamma)^(-sigma)*[theta*B^rho*(3)^(rho-1)+1-delta] == 1;
eqn2 = [(1+n)*(1+gamma)-(1-delta)]*3 == 0.15;
eqn3 = theta*(B^rho)*(3^rho) == 0.33;
if rho == 0
eqn4 = 3^(1-theta)/B == 3;
else
eqn4 = 3/(B*(theta*(3^rho)+(1-theta))^(1/rho)) == 3;
end
paras = vpasolve(eqn1,eqn2,eqn3,eqn4,[betta delta theta B]);
S = ones(4,1);
S(1) = paras.betta;
S(2) = paras.delta;
S(3) = paras.theta;
S(4) = paras.B;
end
function k0 = k0_solve(theta,B,rho)
syms k
if rho == 0
eqn = k^(1-theta)/B == 0.36;
else
eqn = k/(B*(theta*k^rho+(1-theta))^(1/rho)) == 0.36;
end
k0 = vpasolve(eqn, k);
end
function fk1 = fk(k,theta,rho,B)
if rho == 0
fk1 = B*k^theta;
else
fk1 = B*(theta*k^rho+(1-theta))^(1/rho);
end
end

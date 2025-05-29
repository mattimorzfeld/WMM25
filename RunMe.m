clearvars
close all
clc
Colors = brewermap(8,'Dark2');



%% Data
%% ------------------------------------------------------------------------
%d =    [4.1; 2.5;   2537];  % data: vp, vs, rho_m
%s =    [0.2; 0.3;    167];  % standard deviation of the noise
d =    [4.1; 2.5;   2562];  % data: vp, vs, rho_m
s =    [0.2; 0.3;    196];  % standard deviation of the noise

H = [1;1;1]; % select which data to use
nData = sum(H);
AssignData
%% ------------------------------------------------------------------------

%% Model
%% ------------------------------------------------------------------------
% Parameter bounds
%   asp  phi   Water content k mu rho_min
%lb = [0.03  0.05 0  75.6 25.6 2689]';     % lower bound
%ub = [0.99  0.50 1  80   40   2900]';     % upper bound
lb = [0.01  0.01 0  35 20 2600]';     % lower bound
ub = [0.99  0.40 1  80   50   3100]';     % upper bound
n = length(ub);                           % number of unknown parameters
%% ------------------------------------------------------------------------


%% Set up the MCMC and inversion
%% ------------------------------------------------------------------------
% Make a function for the log posterior
logpi=@(x)myLogPi(x,lb,ub,d,s,H);

%% Run the emcee hammer
Xo = zeros(n,1);
numModels = 0;
% initialize with parameters that satisfy all constraints
go = 1;
counter = 0;
while go == 1
    Xo = lb+(ub-lb).*rand(n,1);
    if isfinite(myLogPi(Xo,lb,ub,d,s,H))
        break
    end
end
Xo(3) = randi(2,1)-1;

% Cold start
Nsteps = 1e4;
[X, ~, LogPi,~]=myRWM(Nsteps,Xo,logpi,H);
RMSE = sqrt(-2*LogPi/length(d));
X = X(:,RMSE<3);
Xo = X(:,randi(length(X),1));

% Warm start
Nsteps = 5e5;
[X, D, LogPi,~]=myRWM(Nsteps,Xo,logpi,H);

% MCMC aftermath
BurnIn = 1e4;
X = X(:,BurnIn:end);
D = D(:,BurnIn:end);
LogPi = LogPi(BurnIn:end);
RMSE = sqrt(-2*LogPi/length(d));
%% ------------------------------------------------------------------------


%% Plots
%% ------------------------------------------------------------------------
% PlotScript
%% ------------------------------------------------------------------------

%% Display results
%% --------------------------------------------
m = mean(X,2);
sp = std(X,[],2);

disp(' '), disp(' ')
fprintf('Asp. ratio: %g +/- %g \n',m(1),sp(1))
fprintf('Porosity: %g +/- %g \n',m(2),sp(2))
fprintf('K: %g +/- %g \n',m(4),sp(4))
fprintf('mu: %g +/- %g \n',m(5),sp(5))
fprintf('rho_min: %g +/- %g \n',m(6),sp(6))
%% --------------------------------------------


%% Automatable figure
%% --------------------------------------------
%% --------------------------------------------
probNoH2O = sum(X(3,:)==0)/length(X(3,:));
probH2O = sum(X(3,:)==1)/length(X(3,:));

dryInds = find(X(3,:)==0);
wetInds = find(X(3,:)==1);


% Data fits
figure
subplot(321),hold on
histogram(D(1,:),50,'Normalization','pdf','FaceColor',Colors(5,:))
histogram(D(2,:),50,'Normalization','pdf','FaceColor',Colors(1,:))
errorbar(d(1),0,2*s(1),'horizontal','.','MarkerSize',30,...
    'CapSize',18,'LineWidth',4,'Color',Colors(8,:))
errorbar(d(2),0,2*s(2),'horizontal','.','MarkerSize',30,...
        'CapSize',18,'LineWidth',4,'Color',Colors(7,:))
set(gca,'FontSize',16)
legend('v_p','v_s')
box on
grid off
xlabel('V_p and V_s[km/s]'),ylabel('pdf')
% xlim([3 5])
set(gcf,'Color','w')
fontname("Arial")
 
subplot(322),hold on
histogram(D(3,:),50,'Normalization','pdf','FaceColor',Colors(8,:))
errorbar(d(3),0,2*s(3),'horizontal','k.','MarkerSize',30,...
    'CapSize',18,'LineWidth',4)
set(gca,'FontSize',16)
box on
grid off
xlabel('\rho_b [kg/m^3]'),ylabel('pdf')
set(gcf,'Color','w')
fontname("Arial")

subplot(326),hold on
bar([0 1],[probNoH2O probH2O],'FaceColor','b')
labels = {'Dry', 'Wet'};
set(gca, 'XTick', 1:2, 'XTickLabel', labels)

xlabel('Water saturation')
ylabel('Probability')
box on
grid off
ylim([0 1])
xlim([-.5 1.5])
xticks([0 1])
set(gca,'FontSize',16)
set(gcf,'Color','w')
text(0, probNoH2O, strcat(num2str(round(probNoH2O,2)*100),'%'),'FontSize',16,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(1, probH2O, strcat(num2str(round(probH2O,2)*100),'%'),'FontSize',16,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
fontname("Arial")

subplot(323),hold on
histogram(X(1,dryInds),20,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',Colors(3,:),'LineWidth',2)
histogram(X(1,wetInds),20,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',Colors(2,:),'LineWidth',2)
histogram(X(1,:),20,'Normalization','pdf','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
legend('Dry','Wet','Dry & wet')
box on
grid off
xlabel('\alpha'),ylabel('pdf')
set(gca,'FontSize',16)
set(gcf,'Color','w')
fontname("Arial")

subplot(324),hold on
histogram(X(2,dryInds),20,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',Colors(3,:),'LineWidth',2)
histogram(X(2,wetInds),20,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',Colors(2,:),'LineWidth',2)
histogram(X(2,:),20,'Normalization','pdf','DisplayStyle','stairs','EdgeColor','k','LineWidth',2)
legend('Dry','Wet','Dry & wet','Location','NorthEast')
box on
grid off
xlabel('\phi'),ylabel('pdf')
set(gca,'FontSize',16)
set(gcf,'Color','w')
fontname("Arial")

subplot(325),hold on
ndhist(X(1,:),X(2,:),'filter','bins',0.8);
set(gca,'FontSize',16)
%colormap(crameri('davos'));
colormap sky;
colorbar
box on
xlabel('\alpha')
ylabel('\phi')
set(gca,'FontSize',16)
set(gcf,'Color','w')
box on
fontname("Arial")

%% --------------------------------------------
%% --------------------------------------------

set(gcf,'Position',[608 1 1313 976])

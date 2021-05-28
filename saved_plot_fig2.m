function [] = saved_plot_fig2(save_name)

% Creates a plot as the one reported in Figure 2
% save_name = name of the mat file containing the results from main;

% Load saved data
load(save_name,'Gen_U','Gen_S','mm','NN')


figure()
% Plot the excess risk related with the coarsest scale
la = loglog(mm,Gen_U(1,:),'b-s','LineWidth',1.5); hold on
lb = loglog(mm,Gen_S(1,:),'r-s','LineWidth',1.5);

% Plot the excess risk related with the finest scale
lc = loglog(mm,Gen_U(length(NN),:),'b-p','LineWidth',1.5);
ld = loglog(mm,Gen_S(length(NN),:),'r-p','LineWidth',1.5);

% Expected decay
lf = loglog(mm,0.05*mm.^(-1/2),'k--','LineWidth',1.5);

% Legend
legend([la,lc,lb,ld,lf],['Unsupervised, N = ',num2str(NN(1))],...
    ['Unsupervised, N = ',num2str(NN(end))],['Supervised, N = ',num2str(NN(1))],...
    ['Supervised, N = ',num2str(NN(end))],'m^{-1/2}')

% Further settings
xlim([mm(1),mm(end),])
grid on
xlabel('Sample size, m')
ylabel('Excess risk')
set(gca,'Fontsize',12)

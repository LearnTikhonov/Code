function [] = saved_plot_fig1(save_name,Nk)

% Creates a plot as the one reported in Figure 1
% save_name = name of the mat file containing the results from main;
% Nk =  pick one of the signal resolutions reported in the saved data

% Load saved data
load(save_name,'Gen_U','Gen_S','Gen_U_par','Gen_S_par','mm','NN')

figure()

% Plot the error risk in both scenarios
la = loglog(mm,Gen_U(Nk,:),'b-s','LineWidth',1.5); hold on
lb = loglog(mm,Gen_S(Nk,:),'r-s','LineWidth',1.5);

% Expected decay
lf = loglog(mm,0.05*mm.^(-1/2),'k--','LineWidth',1.5);

% Plot error bars based on standard deviation
svdU = zeros(size(mm));
svdS = zeros(size(mm));
for i = 1:length(mm)
    parNT = Gen_U_par(Nk,i,:);
    svdU(i) = sqrt(var(parNT));
    [~] = loglog([mm(i),mm(i)],[Gen_U(Nk,i)-svdU(i),Gen_U(Nk,i)+svdU(i)],'b-s');
    
    parT = Gen_S_par(Nk,i,:);
    svdS(i) = sqrt(var(parT));
    [~] = loglog([mm(i),mm(i)],[Gen_S(Nk,i)-svdS(i),Gen_S(Nk,i)+svdS(i)],'r-s');
end

% Superimpose colored regions
upnT = (Gen_U(Nk,:)-svdU);
donT = (Gen_U(Nk,:)+svdU);
lkk = fill([(mm'); flipud((mm'))],[donT';flipud(upnT')],[0 0 1], 'linestyle', 'none');
set(lkk,'facealpha',.2)

upT = (Gen_S(Nk,:)-svdS);
doT = (Gen_S(Nk,:)+svdS);
lkk = fill([(mm'); flipud((mm'))],[doT';flipud(upT')],[1 0 0], 'linestyle', 'none');
set(lkk,'facealpha',.2)

% Legend
legend([la,lb,lf],['Unsupervised, N = ',num2str(NN())],...
    ['Supervised, N = ',num2str(NN(Nk))],'m^{-1/2}')

% Further settings
xlim([mm(1),mm(end),])
grid on
xlabel('Sample size, m')
ylabel('Excess risk')
set(gca,'Fontsize',12)

end



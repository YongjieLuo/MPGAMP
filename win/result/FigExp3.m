clear all; close all;
clc;

load DataExp3_N1000.mat;

MP_mean = mean(err_relative{1},2);
MP_std = std(err_relative{1},0,2);
OMP_mean = mean(err_relative{2},2);
OMP_std = std(err_relative{2},0,2);
BP_mean = mean(err_relative{3},2);
BP_std = std(err_relative{3},0,2);
AMP_mean = mean(err_relative{4},2);
AMP_std = std(err_relative{4},0,2);
GAMP_mean = mean(err_relative{5},2);
GAMP_std = std(err_relative{5},0,2);
MPGAMP_mean = mean(err_relative{6},2);
MPGAMP_std = std(err_relative{6},0,2);

xaxis = baseM+(0:addM);

fig_mse = figure;
h_MP = plot(xaxis, MP_mean, 'm-v', 'LineWidth', 1.0);
hold on;
h_OMP = plot(xaxis, OMP_mean, 'k-+', 'LineWidth', 1.0);
hold on;
h_BP = plot(xaxis, BP_mean, 'c-s', 'LineWidth', 1.0);
hold on;
h_AMP = plot(xaxis, AMP_mean, 'g-x', 'LineWidth', 1.0);
hold on;
h_GAMP = plot(xaxis, GAMP_mean, 'y-o', 'LineWidth', 1.0);
hold on;
h_MPGAMP = plot(xaxis, MPGAMP_mean, 'r-^', 'LineWidth', 1.0);


xlabel('Number of Measurements','FontSize',12,'FontName','Arial'); 
ylabel('Average Relative Error','FontSize',12,'FontName','Arial');
box on; grid on;
axis([95,250,-0.1,1.1]);

h_lgd = legend([h_MP, h_OMP, h_BP, h_AMP, h_GAMP, h_MPGAMP],...
                'MP','OMP','BP','AMP','GAMP','MPGAMP','Location','NorthEast');
set(h_lgd, 'FontSize', 12);
set(h_lgd, 'FontName','Arial');
%set(h_lgd, 'Units', 'centimeters');
%pos = get(h_lgd, 'position');
%set(h_lgd, 'position', pos+[0,0,5,5])

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 640 400]);

fig_name = ['FigExp3'];
print(fig_mse, '-depsc', fig_name)

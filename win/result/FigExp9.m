%
clear all; close all;
clc;

load DataExp9_phase_transition

mesh_delta = length(delta);
mesh_rho = length(rho);
M = zeros(mesh_delta,1);
K = zeros(mesh_delta, mesh_rho);

sucesses_GAMP = zeros(mesh_rho, mesh_delta);

threshold = 0.1;

for i = 1:mesh_delta
    M(i) = floor(N*delta(i));
    
    for j = 1:mesh_rho
        K(i,j) = ceil(M(i)*rho(j));
        
        for r = 1:trialNum
            if errL2_relative{2}(i,j,r) < threshold
                sucesses_GAMP(i,j) = sucesses_GAMP(i,j) + 1;
            end
        end
    end
end

sucesses_rate_GAMP = sucesses_GAMP/trialNum;

fig_pt_GAMP = figure
pcolor(delta,rho,sucesses_rate_GAMP); 
colormap (1-gray); colorbar; shading interp; 
xlabel('\delta'); ylabel('\rho');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 500 250]);

fig_name = ['FigExp9'];
print(fig_pt_GAMP, '-depsc', fig_name)




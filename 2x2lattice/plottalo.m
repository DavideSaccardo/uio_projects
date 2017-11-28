function [] = plottalo(name, name_size, xax, yax, ax_size, tick_size)
%
%filename = 'hist_2.400000.txt';
title(name,'Interpreter','teX','Fontsize',name_size);
xlabel(xax,'interpreter','latex','fontsize',ax_size);
ylabel(yax,'interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','tex');
set(gca,'FontSize',tick_size);
set(gca, 'FontName', 'Times');
grid on;
%ll=legend(gca,'show','with Coulomb potential $\omega_r = 5$','without Coulomb potential','location','northeast');
%set(ll,'Interpreter','Latex');
end


%% Make plots
n_plot = floor(20/dt);
subplot(2,3,1);
hold on;
plot(vTime(1:n_plot),monetary_shock(1:n_plot),'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Monetary Policy Shock','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 vTime(n_plot)]);
grid on;
hold off;

subplot(2,3,2);
hold on;
plot(vTime(1:n_plot),inflation(1:n_plot),'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Inflation','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 vTime(n_plot)]);
grid on;
hold off;

subplot(2,3,3);
hold on;
plot(vTime(1:n_plot),consumption(1:n_plot),'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Consumption','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 vTime(n_plot)]);
grid on;
hold off;

subplot(2,3,4);
hold on;
plot(vTime(1:n_plot),Y(1:n_plot),'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('GDP','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 vTime(n_plot)]);
grid on;
hold off;

subplot(2,3,5);
hold on;
plot(vTime(1:n_plot),lab_sup(1:n_plot),'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Labor Supply','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 vTime(n_plot)]);
grid on;
hold off;

subplot(2,3,6);
hold on;
plot(vTime(1:n_plot),wage(1:n_plot),'linewidth',1.5,'linestyle',line_style,'color',color);
set(gcf,'color','w');
title('Wage','interpreter','latex','fontsize',14);
ylabel('$\%$ deviation','interpreter','latex');
xlim([1 vTime(n_plot)]);
grid on;
hold off;

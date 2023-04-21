

%%
figure
maker_idx = 1:2:length(2:2:250);
nn = 2:2:250;
hold on
load TimeScaleNSymmetric2
plot(2:2:250, tau_com,'r-*','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% H/L = 2
hold on
load TimeScaleNSymmetric1
plot(2:2:250, tau_com,'b-d','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% H/L=1
hold on
load TimeScaleNSymmetric3
plot(2:2:250, tau_com,'m-o','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% H/L = 1/1000
hold on
plot(2:2:250,(2+3/4*2)*nn-1-0.91*2,'k-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% H/L = 2
hold on
plot(2:2:250,(2+3/4*1)*nn-1-0.91*1,'k-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% H/L=1
hold on
plot(2:2:250,(2+3/4*1/1000)*nn-1-0.91*1/1000,'k-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% H/L = 1/1000

xlabel('\fontsize{25}\it n','FontName', 'Times New Roman','FontWeight','bold')
ylabel('\fontsize{20} Generalized RC Time','FontName', 'Times New Roman','FontWeight','bold')
legend({'$ H/\mathcal{L} = 2$ ','$ H/\mathcal{L} = 1$','$H/\mathcal{L} =1/1000$',...
    'Ref.29: $(2+\frac{3H}{4\mathcal{L}})n-1-0.91\frac{H}{\mathcal{L}}$'},'Location',...
    'northwest','fontsize',18,'Interpreter','latex')

legend('boxoff')
set(gca,'FontName','Times New Roman','FontSize',18,...
    'GridColor','k','FontWeight','bold','LineWidth',1.5)
box on
set(gcf,'unit','normalized','position',[0.1,0.0,0.450,0.7]);
filename=['fig_tau_symm_com1','.eps'];
print(gcf,'-depsc',filename)

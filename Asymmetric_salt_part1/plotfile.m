%% ...plotting ..................
figure
load AsydataPNP_asymm2000
% load AsydataPNP_asymmError2000
hold on
N_half=4000;
maker_idx = [1:160:length(x),length(x)];
% maker_idx = [1:20:length(x),length(x)];
plot(x(1,1:N_half),pphi(2,1:N_half),'MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=0.5
plot(x(1,1:N_half),pphi(4,1:N_half),'r:','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=2
plot(x(1,1:N_half),pphi(5,1:N_half),'b-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=5
plot(x(1,1:N_half),pphi(6,1:N_half),'k--','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)%t=10
% plot(x(1,1:N_half),pphi(9,1:N_half),'--','LineWidth',2.5)% t=40
% load FDMdataPNPasymError2000
load FDMdataPNPasym2000
plot(x(1,1:N_half),pphi(2,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.)% t=0.5
plot(x(1,1:N_half),pphi(4,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.)% t=2
plot(x(1,1:N_half),pphi(5,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.)% t=5
plot(x(1,1:N_half),pphi(6,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.)%t=10
% plot(x(1,1:N_half),pphi(9,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.5)% t=40
%         title(['\Deltax ',num2str(1/4000)])
xlabel('\fontsize{30}\it x','FontName', 'Times New Roman','FontWeight','bold')
ylabel('\fontsize{25}Potential Distribution','FontName', 'Times New Roman','FontWeight','bold')
legend({'\it t = 0.5 ','\it t = 2','\it t = 5','\it t = 10','Num.'},'Location','southeast','FontSize',25,'FontName', 'Times New Roman','FontWeight','bold')
legend('boxoff')
set(gca,'FontName','Times New Roman','FontSize',22,...
    'GridColor','k','FontWeight','bold','LineWidth',1.5)
box on
%         set(gcf,'unit','normalized','position',[0.1,0.0,0.42,0.6]);
set(gcf,'unit','normalized','position',[0.1,0.0,0.45,0.7]);
ylim([-0.25 0.25])
set(gca,'YTick',[-0.25  -0.1 0 0.1 0.25 ])
str = {'\fontsize{35} (d)'};
dim = [0.17 0 0.0 0.9];% x0,y0,width,height
h = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor',[1 1 1]);
h.LineStyle = 'none';
filename=['fig_potential_asymm_com','.eps'];
print(gcf,'-depsc',filename)


%% ...plotting ..................
figure
load AsydataPNP_asymm2000
hold on
N_half=4000;
maker_idx = [1:80:length(x),length(x)];
% plot(x(1,1:N_half),pp(2,1:N_half),'MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% t=0.5
% plot(x(1,1:N_half),pp(4,1:N_half),'r:','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=2
plot(x(1,1:N_half),pp(5,1:N_half),'MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% t=5
% plot(x(1,1:N_half),pp(6,1:N_half),'b-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)%t=10
% plot(x(1,1:N_half),pphi(9,1:N_half),'--','LineWidth',2.5)% t=40
load FDMdataPNPasym2000
% plot(x(1,1:N_half),pp(2,1:N_half),'-.','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.)% t=0.5
% plot(x(1,1:N_half),pp(4,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.)% t=2
plot(x(1,1:N_half),pp(5,1:N_half),'-.','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',2)% t=5
% plot(x(1,1:N_half),pp(6,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.)%t=10
% plot(x(1,1:N_half),pphi(9,1:N_half),'^','MarkerSize',8,'MarkerEdge',[0,0, 0],'MarkerIndices',maker_idx,'LineWidth',1.5)% t=40
%         title(['\Deltax ',num2str(1/4000)])
xlabel('\fontsize{30}\it x','FontName', 'Times New Roman','FontWeight','bold')
ylabel('\fontsize{20}Cation Concentration','FontName', 'Times New Roman','FontWeight','bold')
legend({'\it t= 5 ','Num.'},'Location','northeast','FontSize',18,'FontName', 'Times New Roman','FontWeight','bold')
legend('boxoff')
set(gca,'FontName','Times New Roman','FontSize',18,...
    'GridColor','k','FontWeight','bold','LineWidth',1.5)
box on
%         set(gcf,'unit','normalized','position',[0.1,0.0,0.42,0.6]);
set(gcf,'unit','normalized','position',[0.1,0.0,0.450,0.7]);
filename=['fig_cation_asymm_com','.eps'];
print(gcf,'-depsc',filename)

%% ...plotting ..................
figure
load AsydataPNP_asymm2000
hold on
N=4000;
maker_idx = [1:90:length(x),length(x)];
plot(x(1,1:N),pphi(1,1:N),'k--','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=0
hold on
plot(x(1,1:N),pphi(2,1:N),'r-d','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=0.5
% plot(x,pphi(3,:),'-s','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2)% t=1
plot(x(1,1:N),pphi(4,1:N),'g-*','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=2
plot(x(1,1:N),pphi(5,1:N),'m-o','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=5
plot(x(1,1:N),pphi(6,1:N),'-+','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)%t=10
plot(x(1,1:N),pphi(7,1:N),'r-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)%t=15
% plot(x(1,1:N),pphi(8,1:N),'-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)%t=20
% plot(x(1,1:N),pphi(9,1:N),'LineWidth',2)% t=40
% plot(x(1,1:N),pphi(10,1:N),'--','LineWidth',2.5)% t=60
% plot(x(1,1:N_half),zeros(1,length(x(1,1:N_half))),'--','LineWidth',2)
%         title(['Error estimator, nx = ',num2str(model.nx)])
xlabel('\fontsize{35}\it x','FontName', 'Times New Roman','FontWeight','bold')
ylabel('\fontsize{28}Potential Distribution','FontName', 'Times New Roman','FontWeight','bold')
legend({'\it t = 0','\it t = 0.5 ','\it t = 2','\it t = 5','\it t = 10','\it t = 15'},'Location','northwest','FontSize',24,'FontName', 'Times New Roman','FontWeight','bold')
legend('boxoff')
set(gca,'FontName','Times New Roman','FontSize',25,...
    'GridColor','k','FontWeight','bold','LineWidth',1.5)
box on
%         set(gcf,'unit','normalized','position',[0.1,0.0,0.42,0.6]);
set(gcf,'unit','normalized','position',[0.1,0.0,0.50,0.7]);
ylim([-0.25 0.25])
set(gca,'YTick',[-0.25  -0.1 0 0.1 0.25 ])
str = {'\fontsize{40} (a)'};
dim = [0.65 0 0.0 0.32];% x0,y0,width,height
h = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor',[1 1 1]);
h.LineStyle = 'none';
filename=['fig_potential_asymm','.eps'];
print(gcf,'-depsc',filename)



%% ...plotting ..................
figure
load AsydataPNP_asymm2000
maker_idx = [1:1:5,5:5:20,60:4:90,90:10:length(tt_real),length(tt_real)];
semilogx(tt_real(1,:),Q_total(1:length(tt_real),1)','d-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_total(1:length(tt_real),2)','rs-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_total(1:length(tt_real),3)','b*-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_total(1:length(tt_real),4)','m^-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_total(1:length(tt_real),5)','o-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
% for i_n=1:floor(n_ed/2)%length(t_store)+1
% %     plot(tt_real(1,:),Q_total(1:length(tt_real),i_n)','LineWidth',2)
% semilogx(tt_real(1,:),Q_total(1:length(tt_real),i_n)','LineWidth',2)
% hold on
% end
%         title(['Asymmetric salt ','2,', '-1'])
xlabel('\fontsize{32}Time','FontName', 'Times New Roman','FontWeight','bold')
ylabel('\fontsize{28}Total Diffuse Charge','FontName', 'Times New Roman','FontWeight','bold')
legend({'\it i = 1','\it i = 2', '\it i = 3','\it i = 4','\it i = 5'},'Location','northwest','FontSize',25,'FontName', 'Times New Roman','FontWeight','bold')
legend('boxoff')
% xlim([0 20])
set(gca,'FontName','Times New Roman','FontSize',25,...
    'GridColor','k','FontWeight','bold','LineWidth',1.5)
box on
set(gcf,'unit','normalized','position',[0.1,0.0,0.50,0.7]);
xlim([0.01 20])
xx =[0.01 0.1  1 10 20];
set(gca,'XTick',[0.01 0.1  1 10 20],'XTickLabel',xx)
ylim([-0.001 0.006])
% yy = [-0.001 0 0.002 0.004 0.006 ];
set(gca,'YTick',[-0.001 0 0.002 0.004 0.006 ])
str = {'\fontsize{40} (b)'};
dim = [0.65 0 0.0 0.32];% x0,y0,width,height
h = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor',[1 1 1]);
h.LineStyle = 'none';
filename=['fig_totalcharge_asymmL','.eps'];
print(gcf,'-depsc',filename)


%% ...plotting ..................
figure
load AsydataPNP_asymm2000
maker_idx = [1:1:5,5:10:60,60:4:90,90:20:length(tt_real),length(tt_real)];
semilogx(tt_real(1,:),Q_totalR(1:length(tt_real),1)','d-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_totalR(1:length(tt_real),2)','rs-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_totalR(1:length(tt_real),3)','b*-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_totalR(1:length(tt_real),4)','m^-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
semilogx(tt_real(1,:),Q_totalR(1:length(tt_real),5)','o-','MarkerSize',10,'MarkerIndices',maker_idx,'LineWidth',2.5)
hold on
% for i_n=1:floor(n_ed/2)%length(t_store)+1
% %     plot(tt_real(1,:),Q_total(1:length(tt_real),i_n)','LineWidth',2)
% semilogx(tt_real(1,:),Q_total(1:length(tt_real),i_n)','LineWidth',2)
% hold on
% end
%         title(['Asymmetric salt ','2,', '-1'])
xlabel('\fontsize{32}Time','FontName', 'Times New Roman','FontWeight','bold')
ylabel('\fontsize{28}Total Diffuse Charge','FontName', 'Times New Roman','FontWeight','bold')
legend({'\it i = 1','\it i = 2', '\it i = 3','\it i = 4','\it i = 5'},'Location','southwest','FontSize',25,'FontName', 'Times New Roman','FontWeight','bold')
legend('boxoff')
set(gca,'FontName','Times New Roman','FontSize',25,...
    'GridColor','k','FontWeight','bold','LineWidth',1.5)
box on
set(gcf,'unit','normalized','position',[0.1,0.0,0.5,0.7]);
xlim([0.01 20])
xx =[0.01 0.1  1 10 20];
set(gca,'XTick',[0.01 0.1  1 10 20],'XTickLabel',xx)
ylim([-0.006 0.001])
% yy = [-0.006 -0.004 -0.002 0  0.001 ];
set(gca,'YTick',[-0.006 -0.004 -0.002 0  0.001 ])
str = {'\fontsize{40} (c)'};
dim = [0.65 0 0.0 0.32];% x0,y0,width,height
h = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor',[1 1 1]);
h.LineStyle = 'none';
filename=['fig_totalcharge_asymmR','.eps'];
print(gcf,'-depsc',filename)



%% calculate errors
load FDMdataPNPasymError2000
t_store =tt_real2;
N_half =4000;
pp_all(1:length(t_store),:) = pp(2:length(t_store)+1,:);
nn_all(1:length(t_store),:) = nn(2:length(t_store)+1,:);
pphi_all(1:length(t_store),:) = pphi(2:length(t_store)+1,:);
load AsydataPNP_asymmError2000
pp_all(length(t_store)+1:2*(length(t_store)),:) = pp(2:length(t_store)+1,:);
nn_all(length(t_store)+1:2*(length(t_store)),:) = nn(2:length(t_store)+1,:);
pphi_all(length(t_store)+1:2*(length(t_store)),:) = pphi(2:length(t_store)+1,:);
Err_pp = max(abs(pp_all(1:length(t_store),1:N_half) - pp_all(length(t_store)+1:2*(length(t_store)),1:N_half)),[],2);%./max(abs(pp_all(1:length(t_store),1:N_half)),[],2);
Err_nn = max(abs(nn_all(1:length(t_store),1:N_half) - nn_all(length(t_store)+1:2*(length(t_store)),1:N_half)),[],2);%./max(abs(nn_all(1:length(t_store),1:N_half)),[],2);
Err_pphi = max(abs(pphi_all(1:length(t_store),1:N_half) - pphi_all(length(t_store)+1:2*(length(t_store)),1:N_half)),[],2);%./max(abs(pphi_all(1:length(t_store),1:N_half)),[],2);
figure
% t_store = [0.5,1,2,5,10,15,20,40,60,70,80,90,100];
maker_idx = [1:10:length(t_store),length(t_store)];
semilogy(t_store(1,1:end),Err_pp(1:end,1),'r--','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=0.5
hold on
semilogy(t_store(1,1:end),Err_nn(1:end,1),'b-.','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=0.5
hold on
semilogy(t_store(1,1:end),Err_pphi(1:end,1),'k','MarkerSize',8,'MarkerIndices',maker_idx,'LineWidth',2.5)% t=0.5
%         title(['Error estimator, nx = ',num2str(model.nx)])
%        title('Error of Potential distribution')
xlabel('\fontsize{28} Time','FontName', 'Times New Roman','FontWeight','bold')
ylabel('\fontsize{25} Error','FontName', 'Times New Roman','FontWeight','bold')
legend({'\it c_+', '\it c_-','\it \Phi'},'Location','southeast','FontSize',28,'FontName', 'Times New Roman','FontWeight','bold')
legend('boxoff')
set(gca,'FontName','Times New Roman','FontSize',22,...
    'GridColor','k','FontWeight','bold','LineWidth',1.5)
set(gcf,'unit','normalized','position',[0.1,0.0,0.45,0.7]);
ylim([0.0001 0.1])
% set(gca,'YTick',[0.0001 0.0005  0.005 0.05 ])
str = {'\fontsize{35} (c)'};
dim = [0.16 0 0.0 0.9];% x0,y0,width,height
h = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor',[1 1 1]);
h.LineStyle = 'none';
filename=['fig_asymm_abserr','.eps'];
print(gcf,'-depsc',filename)
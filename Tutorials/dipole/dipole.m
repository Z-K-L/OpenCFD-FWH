clc;
clear
%% Stationary dipole in a moving medium with AoA
format long
data=importdata('./FWH-dipole/matlab.control');
c0=data(1)
M0=data(2)
rou0=data(3)
AoA=data(4) 
AoA=AoA*pi/180;
M1=M0*cos(AoA);
M2=M0*sin(AoA);
A=data(5) %振幅 m^2/s
f=data(6) %频率
omiga=2*pi*f;
n_t=round(data(7))
t(:)=linspace(2/n_t,2,n_t);

r_FWH=round(data(8))
r_Obs=round(data(9))
theta=pi/2;
phi=linspace(0,2*pi,201);

x1 = zeros(length(theta), length(phi));
x2 = zeros(length(theta), length(phi));
x3 = zeros(length(theta), length(phi));
for ii=1:length(theta)
    for j=1:length(phi)
x1(ii,j)=r_Obs*sin(theta(ii))*cos(phi(j));
x2(ii,j)=r_Obs*sin(theta(ii))*sin(phi(j));
x3(ii,j)=r_Obs*cos(theta(ii));
    end
end

dx1 = 1e-3;
dx2 = 1e-3;

R_s_mx2 = zeros(length(theta), length(phi));
R_mx2 = zeros(length(theta), length(phi));
R_s_px2 = zeros(length(theta), length(phi));
R_px2 = zeros(length(theta), length(phi));
R_spx1_px2= zeros(length(theta), length(phi));
Rpx1_px2= zeros(length(theta), length(phi));
R_spx1_mx2= zeros(length(theta), length(phi));
Rpx1_mx2= zeros(length(theta), length(phi));
R_smx1_px2= zeros(length(theta), length(phi));
Rmx1_px2= zeros(length(theta), length(phi));
R_smx1_mx2= zeros(length(theta), length(phi));
Rmx1_mx2= zeros(length(theta), length(phi));
R_spx2_px2= zeros(length(theta), length(phi));
Rpx2_px2= zeros(length(theta), length(phi));
R_spx2_mx2= zeros(length(theta), length(phi));
Rpx2_mx2= zeros(length(theta), length(phi));
R_smx2_px2= zeros(length(theta), length(phi));
Rmx2_px2= zeros(length(theta), length(phi));
R_smx2_mx2= zeros(length(theta), length(phi));
Rmx2_mx2= zeros(length(theta), length(phi));
for ii=1:length(theta)
    for j=1:length(phi)
R_s_px2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)+dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)+dx2).^2+x3(ii,j).^2));
R_px2(ii,j)=(-M1*x1(ii,j)-M2*(x2(ii,j)+dx2)+R_s_px2(ii,j))/(1-M0^2);
R_s_mx2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)-dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)-dx2).^2+x3(ii,j).^2));
R_mx2(ii,j)=(-M1*x1(ii,j)-M2*(x2(ii,j)-dx2)+R_s_mx2(ii,j))/(1-M0^2);


R_spx1_px2(ii,j)= sqrt((M1*(x1(ii,j)+dx1)+M2*(x2(ii,j)+dx2)).^2+(1-M0^2)*((x1(ii,j)+dx1).^2+(x2(ii,j)+dx2).^2+x3(ii,j).^2));
Rpx1_px2(ii,j)= (-M1*(x1(ii,j)+dx1)-M2*(x2(ii,j)+dx2)+R_spx1_px2(ii,j))/(1-M0^2);
R_spx1_mx2(ii,j)= sqrt((M1*(x1(ii,j)+dx1)+M2*(x2(ii,j)-dx2)).^2+(1-M0^2)*((x1(ii,j)+dx1).^2+(x2(ii,j)-dx2).^2+x3(ii,j).^2));
Rpx1_mx2(ii,j)= (-M1*(x1(ii,j)+dx1)-M2*(x2(ii,j)-dx2)+R_spx1_mx2(ii,j))/(1-M0^2);

R_smx1_px2(ii,j)= sqrt((M1*(x1(ii,j)-dx1)+M2*(x2(ii,j)+dx2)).^2+(1-M0^2)*((x1(ii,j)-dx1).^2+(x2(ii,j)+dx2).^2+x3(ii,j).^2));
Rmx1_px2(ii,j)= (-M1*(x1(ii,j)-dx1)-M2*(x2(ii,j)+dx2)+R_smx1_px2(ii,j))/(1-M0^2);
R_smx1_mx2(ii,j)= sqrt((M1*(x1(ii,j)-dx1)+M2*(x2(ii,j)-dx2)).^2+(1-M0^2)*((x1(ii,j)-dx1).^2+(x2(ii,j)-dx2).^2+x3(ii,j).^2));
Rmx1_mx2(ii,j)= (-M1*(x1(ii,j)-dx1)-M2*(x2(ii,j)-dx2)+R_smx1_mx2(ii,j))/(1-M0^2);

R_spx2_px2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)+dx2+dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)+dx2+dx2).^2+x3(ii,j).^2));
Rpx2_px2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)+dx2+dx2)+R_spx2_px2(ii,j))/(1-M0^2);
R_spx2_mx2(ii,j)= sqrt((M1*x1(ii,j)+M2*x2(ii,j)).^2+(1-M0^2)*(x1(ii,j).^2+x2(ii,j).^2+x3(ii,j).^2));
Rpx2_mx2(ii,j)= (-M1*x1(ii,j)-M2*x2(ii,j)+R_spx2_mx2(ii,j))/(1-M0^2);

R_smx2_px2(ii,j)= sqrt((M1*x1(ii,j)+M2*x2(ii,j)).^2+(1-M0^2)*(x1(ii,j).^2+x2(ii,j).^2+x3(ii,j).^2));
Rmx2_px2(ii,j)= (-M1*x1(ii,j)-M2*x2(ii,j)+R_smx2_px2(ii,j))/(1-M0^2);
R_smx2_mx2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)-dx2-dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)-dx2-dx2).^2+x3(ii,j).^2));
Rmx2_mx2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)-dx2-dx2)+R_smx2_mx2(ii,j))/(1-M0^2);
    end
end

fai0_mx2 = zeros(length(theta), length(phi),length(t));
fai0_px2 = zeros(length(theta), length(phi),length(t));
fai0px1_px2 = zeros(length(theta), length(phi),length(t));
fai0px1_mx2 = zeros(length(theta), length(phi),length(t));
fai0mx1_px2 = zeros(length(theta), length(phi),length(t));
fai0mx1_mx2 = zeros(length(theta), length(phi),length(t));
fai0px2_px2 = zeros(length(theta), length(phi),length(t));
fai0px2_mx2 = zeros(length(theta), length(phi),length(t));
fai0mx2_px2 = zeros(length(theta), length(phi),length(t));
fai0mx2_mx2 = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
    for j=1:length(phi)
        for k=1:length(t)
  fai0_px2(ii,j,k)=A/(4*pi*R_s_px2(ii,j))*exp(1i*omiga*(t(k)-R_px2(ii,j)/c0));
  fai0_mx2(ii,j,k)=A/(4*pi*R_s_mx2(ii,j))*exp(1i*omiga*(t(k)-R_mx2(ii,j)/c0));

  fai0px1_px2(ii,j,k)=A/(4*pi*R_spx1_px2(ii,j))*exp(1i*omiga*(t(k)-Rpx1_px2(ii,j)/c0));
  fai0px1_mx2(ii,j,k)=A/(4*pi*R_spx1_mx2(ii,j))*exp(1i*omiga*(t(k)-Rpx1_mx2(ii,j)/c0));

  fai0mx1_px2(ii,j,k)=A/(4*pi*R_smx1_px2(ii,j))*exp(1i*omiga*(t(k)-Rmx1_px2(ii,j)/c0));
  fai0mx1_mx2(ii,j,k)=A/(4*pi*R_smx1_mx2(ii,j))*exp(1i*omiga*(t(k)-Rmx1_mx2(ii,j)/c0));

  fai0px2_px2(ii,j,k)=A/(4*pi*R_spx2_px2(ii,j))*exp(1i*omiga*(t(k)-Rpx2_px2(ii,j)/c0));
  fai0px2_mx2(ii,j,k)=A/(4*pi*R_spx2_mx2(ii,j))*exp(1i*omiga*(t(k)-Rpx2_mx2(ii,j)/c0));

  fai0mx2_px2(ii,j,k)=A/(4*pi*R_smx2_px2(ii,j))*exp(1i*omiga*(t(k)-Rmx2_px2(ii,j)/c0));
  fai0mx2_mx2(ii,j,k)=A/(4*pi*R_smx2_mx2(ii,j))*exp(1i*omiga*(t(k)-Rmx2_mx2(ii,j)/c0));
        end
    end
end
clear R_s_mx2;clear R_mx2;clear R_s_px2;clear R_px2;
clear R_spx1_px2;clear Rpx1_px2;clear R_spx1_mx2;clear Rpx1_mx2;
clear R_smx1_px2;clear Rmx1_px2;clear R_smx1_mx2;clear Rmx1_mx2;
clear R_spx2_px2;clear Rpx2_px2;clear R_spx2_mx2;clear Rpx2_mx2;
clear R_smx2_px2;clear Rmx2_px2;clear R_smx2_mx2;clear Rmx2_mx2;
%% compute potential
fai = zeros(length(theta), length(phi),length(t));
faipx1 = zeros(length(theta), length(phi),length(t));
faimx1 = zeros(length(theta), length(phi),length(t));
faipx2 = zeros(length(theta), length(phi),length(t));
faimx2 = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    for k=1:length(t)
    fai(ii,j,k)= (fai0_px2(ii,j,k)-fai0_mx2(ii,j,k))/(2*dx2);
    faipx1(ii,j,k)= (fai0px1_px2(ii,j,k)-fai0px1_mx2(ii,j,k))/(2*dx2);
    faimx1(ii,j,k)= (fai0mx1_px2(ii,j,k)-fai0mx1_mx2(ii,j,k))/(2*dx2);
    faipx2(ii,j,k)= (fai0px2_px2(ii,j,k)-fai0px2_mx2(ii,j,k))/(2*dx2);
    faimx2(ii,j,k)= (fai0mx2_px2(ii,j,k)-fai0mx2_mx2(ii,j,k))/(2*dx2);
    end
  end
end
clear fai0_px2;clear fai0_mx2;
clear fai0px1_px2;clear fai0px1_mx2;
clear fai0mx1_px2;clear fai0mx1_mx2;
clear fai0px2_px2;clear fai0px2_mx2;
clear fai0mx2_px2;clear fai0mx2_mx2;
%% compute Vel
% u = zeros(length(theta), length(phi),length(t));
% v = zeros(length(theta), length(phi),length(t));
% w = zeros(length(theta), length(phi),length(t));
% for ii=1:length(theta)
%     for k=1:length(t)
%     u(ii,:,k)= real(gradient(fai(ii,:,k), x1(ii,:)));
%     v(ii,:,k)= real(gradient(fai(ii,:,k), x2(ii,:)));
%     w(ii,:,k)= real(gradient(fai(ii,:,k), x3(ii,:)));
%     end
% end
%% compute p 
p1 = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    for k=1:length(t)
    p1(ii,j,k)= -rou0*(M1*c0*(faipx1(ii,j,k)-faimx1(ii,j,k))/(2*dx1)+ M2*c0*(faipx2(ii,j,k)-faimx2(ii,j,k))/(2*dx2));
    end
  end
end
clear faipx1;clear faimx1;clear faipx2;clear faimx2;

pt = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    tmp(:)=fai(ii,j,:);
    pt(ii,j,:)= -rou0*gradient(tmp, t);
  end
end
clear tmp;
clear fai;

p=real(p1+pt);
clear p1;clear pt;
% rou=p/c0^2;
%% compute RMSP

RMSP = zeros(length(theta), length(phi));
for ii=1:length(theta)
    for j=1:length(phi)
    tmp(:)=p(ii,j,:);
    RMSP(ii,j) = rms(tmp);
    end
end
clear tmp
figure(1)
set(groot,'defaultLineLineWidth',3.5);
p1=polarplot(phi,RMSP(1,:));
hold on
%读取FWH结果并绘图
n_Obs=round(data(10))
U0=M0*c0;
for ii=1:n_Obs
    file_number = sprintf('%03d', ii); 
    file_name = strcat('./FWH-dipole/FWH_result-mpi/p_observer-',file_number, '.log');
    FWH(ii,:,:)=importdata(file_name);
    FWH_t(ii,:)=FWH(ii,:,1)/U0;
    FWH_p(ii,:)=FWH(ii,:,2)*rou0*U0^2;
    RMSP_FWH(ii) = rms(FWH_p(ii,:));
end
p2=polarplot(linspace(2*pi/n_Obs,2*pi,n_Obs),RMSP_FWH,'o','LineWidth',4,'MarkerSize',8);
rlim([0 0.027])
set(gca,'FontName','Times New Roman');
set(gca,'linewidth',3,'fontsize',25);
legend([p1 p2],'analytic solution','OpenCFD-FWH_WT solution', ...
    'Interpreter','none','linewidth',2,'fontsize',25)
% close(1);

figure(2)
tmp_p=p(1,(length(phi)-1)/4*3+1,:);
p1=plot(t,tmp_p(:),'-*');
hold on

p2=plot(FWH_t(round(n_Obs*0.75),:),FWH_p(round(n_Obs*0.75),:),'o','LineWidth',4,'MarkerSize',8);
set(groot,'defaultLineLineWidth',5);
set(gca,'FontName','Times New Roman');
set(gca,'linewidth',3,'fontsize',35);
legend([p1 p2],'analytic result','OpenCFD-FWH_WT result', 'Interpreter', ...
    'none','linewidth',2,'fontsize',25,'location','bestoutside')
box off
xlim([0.5 1.5])
xlabel('\itt/s');
ylabel("\itp'/Pa");
% close(2);




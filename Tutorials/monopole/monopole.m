clc;
clear
%% Stationary monopole in a moving medium with AoA
format long
data=importdata('./FWH-monopole/matlab.control');
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
t(:)=linspace(2.5/n_t,2.5,n_t);

r_FWH=round(data(8))
r_Obs=round(data(9))
theta=pi/2;
phi=linspace(0,2*pi,145);

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

R_s = zeros(length(theta), length(phi));
R = zeros(length(theta), length(phi));
dx1 = 1e-10;
dx2 = 1e-10;

R_spx1= zeros(length(theta), length(phi));
Rpx1= zeros(length(theta), length(phi));
R_smx1= zeros(length(theta), length(phi));
Rmx1= zeros(length(theta), length(phi));
R_spx2= zeros(length(theta), length(phi));
Rpx2= zeros(length(theta), length(phi));
R_smx2= zeros(length(theta), length(phi));
Rmx2= zeros(length(theta), length(phi));
for ii=1:length(theta)
    for j=1:length(phi)
R_s(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j))).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)).^2+x3(ii,j).^2));
R(ii,j)=(-M1*x1(ii,j)-M2*(x2(ii,j))+R_s(ii,j))/(1-M0^2);


R_spx1(ii,j)= sqrt((M1*(x1(ii,j)+dx1)+M2*(x2(ii,j))).^2+(1-M0^2)*((x1(ii,j)+dx1).^2+(x2(ii,j)).^2+x3(ii,j).^2));
Rpx1(ii,j)= (-M1*(x1(ii,j)+dx1)-M2*(x2(ii,j))+R_spx1(ii,j))/(1-M0^2);

R_smx1(ii,j)= sqrt((M1*(x1(ii,j)-dx1)+M2*(x2(ii,j))).^2+(1-M0^2)*((x1(ii,j)-dx1).^2+(x2(ii,j)).^2+x3(ii,j).^2));
Rmx1(ii,j)= (-M1*(x1(ii,j)-dx1)-M2*(x2(ii,j))+R_smx1(ii,j))/(1-M0^2);

R_spx2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)+dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)+dx2).^2+x3(ii,j).^2));
Rpx2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)+dx2)+R_spx2(ii,j))/(1-M0^2);

R_smx2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)-dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)-dx2).^2+x3(ii,j).^2));
Rmx2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)-dx2)+R_smx2(ii,j))/(1-M0^2);
    end
end

%% compute potential
fai = zeros(length(theta), length(phi),length(t));
faipx1 = zeros(length(theta), length(phi),length(t));
faimx1 = zeros(length(theta), length(phi),length(t));
faipx2 = zeros(length(theta), length(phi),length(t));
faimx2 = zeros(length(theta), length(phi),length(t));
faipx3 = zeros(length(theta), length(phi),length(t));
faimx3 = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    for k=1:length(t)
    fai(ii,j,k)= A/(4*pi*R_s(ii,j))*exp(1i*omiga*(t(k)-R(ii,j)/c0));
    faipx1(ii,j,k)= A/(4*pi*R_spx1(ii,j))*exp(1i*omiga*(t(k)-Rpx1(ii,j)/c0));
    faimx1(ii,j,k)= A/(4*pi*R_smx1(ii,j))*exp(1i*omiga*(t(k)-Rmx1(ii,j)/c0));
    faipx2(ii,j,k)= A/(4*pi*R_spx2(ii,j))*exp(1i*omiga*(t(k)-Rpx2(ii,j)/c0));
    faimx2(ii,j,k)= A/(4*pi*R_smx2(ii,j))*exp(1i*omiga*(t(k)-Rmx2(ii,j)/c0));
    end
  end
end
clear R_s;
clear R_spx1;clear Rpx1
clear R_smx1;clear Rmx1;
clear R_spx2;clear Rpx2;
clear R_smx2;clear Rmx2;
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
    pt(ii,j,:)= -rou0*1i*omiga*fai(ii,j,:);
  end
end
clear tmp;
clear fai;

p=real(p1+pt);
clear p1;clear pt;

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
set(groot,'defaultLineLineWidth',5);
p1=polarplot(phi,RMSP(1,:));
hold on
%读取FWH结果并绘图
n_Obs=round(data(10))
U0=M0*c0;
for ii=1:n_Obs
    file_number = sprintf('%03d', ii); 
    file_name = strcat('./FWH-monopole/FWH_result-mpi/p_observer-',file_number, '.log');
    FWH(ii,:,:)=importdata(file_name);
    FWH_t(ii,:)=FWH(ii,:,1)/U0; 
    FWH_p(ii,:)=FWH(ii,:,2)*rou0*U0^2;
    RMSP_FWH(ii) = rms(FWH_p(ii,:));
end
p2=polarplot(linspace(2*pi/n_Obs,2*pi,n_Obs),RMSP_FWH,'o','LineWidth',4,'MarkerSize',8);
rlim([0 0.12])
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
xlim([0.75 2.25])
xlabel('\itt/s');
ylabel("\itp'/Pa");
clc;
clear;
tic;
%% Stationary dipole in a moving medium with AoA
format long
c0=340;
M0=0.5;
rou0=1;
AoA=10; AoA=AoA*pi/180;
M1=M0*cos(AoA);
M2=M0*sin(AoA);
A=1.5; %振幅 m^2/s
f=7.5; %频率
omiga=2*pi*f;
n_t=1e3;
t(:)=linspace(2/n_t,2,n_t);

r_FWH=2;
r_Obs=50;
n_Obs=20;
n_theta=18;
theta=linspace(pi/2/n_theta,pi-pi/2/n_theta,n_theta); % cell center
% theta*180/pi
n_phi=36;
phi=linspace(pi/n_phi,2*pi-pi/n_phi,n_phi); % cell center
% phi*180/pi

% theta=pi/2;
% phi=linspace(0,2*pi,201);

x1 = zeros(length(theta), length(phi));
x2 = zeros(length(theta), length(phi));
x3 = zeros(length(theta), length(phi));
for ii=1:length(theta)
    for j=1:length(phi)
x1(ii,j)=r_FWH*sin(theta(ii))*cos(phi(j));
x2(ii,j)=r_FWH*sin(theta(ii))*sin(phi(j));
x3(ii,j)=r_FWH*cos(theta(ii));
    end
end
% set(groot,'defaultLineLineWidth',1);
% plot3(x1,x2,x3,'-o')


R_s = zeros(length(theta), length(phi));
R = zeros(length(theta), length(phi));
dx1 = 1e-3;
dx2 = 1e-3;
dx3 = 1e-3;

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
R_spx3_px2= zeros(length(theta), length(phi));
Rpx3_px2= zeros(length(theta), length(phi));
R_spx3_mx2= zeros(length(theta), length(phi));
Rpx3_mx2= zeros(length(theta), length(phi));
R_smx3_px2= zeros(length(theta), length(phi));
Rmx3_px2= zeros(length(theta), length(phi));
R_smx3_mx2= zeros(length(theta), length(phi));
Rmx3_mx2= zeros(length(theta), length(phi));
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

R_spx3_px2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)+dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)+dx2).^2+(x3(ii,j)+dx3).^2));
Rpx3_px2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)+dx2)+R_spx3_px2(ii,j))/(1-M0^2);
R_spx3_mx2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)-dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)-dx2).^2+(x3(ii,j)+dx3).^2));
Rpx3_mx2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)-dx2)+R_spx3_mx2(ii,j))/(1-M0^2);

R_smx3_px2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)+dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)+dx2).^2+(x3(ii,j)-dx3).^2));
Rmx3_px2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)+dx2)+R_smx3_px2(ii,j))/(1-M0^2);
R_smx3_mx2(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j)-dx2)).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)-dx2).^2+(x3(ii,j)-dx3).^2));
Rmx3_mx2(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j)-dx2)+R_smx3_mx2(ii,j))/(1-M0^2);

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
fai0px3_px2 = zeros(length(theta), length(phi),length(t));
fai0px3_mx2 = zeros(length(theta), length(phi),length(t));
fai0mx3_px2 = zeros(length(theta), length(phi),length(t));
fai0mx3_mx2 = zeros(length(theta), length(phi),length(t));
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
  
  fai0px3_px2(ii,j,k)=A/(4*pi*R_spx3_px2(ii,j))*exp(1i*omiga*(t(k)-Rpx3_px2(ii,j)/c0));
  fai0px3_mx2(ii,j,k)=A/(4*pi*R_spx3_mx2(ii,j))*exp(1i*omiga*(t(k)-Rpx3_mx2(ii,j)/c0));

  fai0mx3_px2(ii,j,k)=A/(4*pi*R_smx3_px2(ii,j))*exp(1i*omiga*(t(k)-Rmx3_px2(ii,j)/c0));
  fai0mx3_mx2(ii,j,k)=A/(4*pi*R_smx3_mx2(ii,j))*exp(1i*omiga*(t(k)-Rmx3_mx2(ii,j)/c0));
        end
    end
end
a=1
clear R_s_mx2;clear R_mx2;clear R_s_px2;clear R_px2;
clear R_spx1_px2;clear Rpx1_px2;clear R_spx1_mx2;clear Rpx1_mx2;
clear R_smx1_px2;clear Rmx1_px2;clear R_smx1_mx2;clear Rmx1_mx2;
clear R_spx2_px2;clear Rpx2_px2;clear R_spx2_mx2;clear Rpx2_mx2;
clear R_smx2_px2;clear Rmx2_px2;clear R_smx2_mx2;clear Rmx2_mx2;
clear R_spx3_px2;clear Rpx3_px2;clear R_spx3_mx2;clear Rpx3_mx2;
clear R_smx3_px2;clear Rmx3_px2;clear R_smx3_mx2;clear Rmx3_mx2;
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
    fai(ii,j,k)= (fai0_px2(ii,j,k)-fai0_mx2(ii,j,k))/(2*dx2);
    faipx1(ii,j,k)= (fai0px1_px2(ii,j,k)-fai0px1_mx2(ii,j,k))/(2*dx2);
    faimx1(ii,j,k)= (fai0mx1_px2(ii,j,k)-fai0mx1_mx2(ii,j,k))/(2*dx2);
    faipx2(ii,j,k)= (fai0px2_px2(ii,j,k)-fai0px2_mx2(ii,j,k))/(2*dx2);
    faimx2(ii,j,k)= (fai0mx2_px2(ii,j,k)-fai0mx2_mx2(ii,j,k))/(2*dx2);
    faipx3(ii,j,k)= (fai0px3_px2(ii,j,k)-fai0px3_mx2(ii,j,k))/(2*dx2);
    faimx3(ii,j,k)= (fai0mx3_px2(ii,j,k)-fai0mx3_mx2(ii,j,k))/(2*dx2);
    end
  end
end

a=2
clear fai0_px2;clear fai0_mx2;
clear fai0px1_px2;clear fai0px1_mx2;
clear fai0mx1_px2;clear fai0mx1_mx2;
clear fai0px2_px2;clear fai0px2_mx2;
clear fai0mx2_px2;clear fai0mx2_mx2;
clear fai0px3_px2;clear fai0px3_mx2;
clear fai0mx3_px2;clear fai0mx3_mx2;
%% compute Vel
u = zeros(length(theta), length(phi),length(t));
v = zeros(length(theta), length(phi),length(t));
w = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    for k=1:length(t)
    u(ii,j,k)= real((faipx1(ii,j,k)-faimx1(ii,j,k))/(2*dx1));
    v(ii,j,k)= real((faipx2(ii,j,k)-faimx2(ii,j,k))/(2*dx2));
    w(ii,j,k)= real((faipx3(ii,j,k)-faimx3(ii,j,k))/(2*dx3));
    end
  end
end
a=3
clear faipx1;clear faimx1;clear faipx2;clear faimx2;clear faipx3;clear faimx3;

%% compute p & rou
p1 = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    for k=1:length(t)
    p1(ii,j,k)= -rou0*(M1*c0*u(ii,j,k)+ M2*c0*v(ii,j,k));
    end
  end
end

a=4
pt = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    tmp(:)=fai(ii,j,:);
    pt(ii,j,:)= real(-rou0*gradient(tmp, t));
  end
end
clear tmp;
clear fai;
a=5
p=p1+pt;
clear p1;clear pt;
rou=p/c0^2;
%% 写入matlab.control
fileID = fopen('./FWH-dipole/matlab.control', 'w'); % 打开文件，以写入模式打开  
fprintf(fileID, '%f\n', c0);
fprintf(fileID, '%f\n', M0);
fprintf(fileID, '%f\n', rou0);
fprintf(fileID, '%f\n', AoA*180/pi);
fprintf(fileID, '%f\n', A);
fprintf(fileID, '%f\n', f);
fprintf(fileID, '%f\n', n_t);
fprintf(fileID, '%f\n', r_FWH);
fprintf(fileID, '%f\n', r_Obs);
fprintf(fileID, '%f\n', n_Obs);
fclose(fileID);
%% Observers.dat.dat
fileID = fopen('./FWH-dipole/Observers.dat', 'w'); % 打开文件，以写入模式打开  
theta_Obs=pi/2;
phi_Obs=linspace(2*pi/n_Obs,2*pi,n_Obs);
for j=1:length(phi_Obs)
    fprintf(fileID, '%25.16E', r_Obs*sin(theta_Obs)*cos(phi_Obs(j)));  
    fprintf(fileID, '%25.16E', r_Obs*sin(theta_Obs)*sin(phi_Obs(j)));  
%     fprintf(fileID, '%25.16E\n', r_Obs*cos(theta_Obs));  
    fprintf(fileID, '  %.1fd0\n', 0);
%     if j>1
%         plot(r_Obs*sin(theta_Obs)*cos(phi_Obs(j)),r_Obs*sin(theta_Obs)*sin(phi_Obs(j)),"o")
%         hold on
%     else
%         plot(r_Obs*sin(theta_Obs)*cos(phi_Obs(j)),r_Obs*sin(theta_Obs)*sin(phi_Obs(j)),"*")
%         hold on
%     end
end
fclose(fileID);
% axis equal
% legend(Location="bestoutside")
%% 写入control.fwh
U0=M0*c0;
fileID = fopen('./FWH-dipole/control.fwh', 'w'); % 打开文件，以写入模式打开  
fprintf(fileID, '%s\n', '$control_FWH');  
fprintf(fileID, '  %s', 'Ma=');fprintf(fileID, '%f', M0);fprintf(fileID, '%s\n', 'd0');
fprintf(fileID, '  %s', 'AoA=');fprintf(fileID, '%f', AoA*180/pi);fprintf(fileID, '%s\n', 'd0');
fprintf(fileID, '  %s', 'delta_t=');fprintf(fileID, '%f', (t(2)-t(1))/(1/U0));fprintf(fileID, '%s\n', 'd0');
fprintf(fileID, '  %s', 'Num_Obs=');fprintf(fileID, '%d\n', n_Obs);
fprintf(fileID, '  %s', 'Kstep_start=');fprintf(fileID, '%d\n', 1);
fprintf(fileID, '  %s', 'Kstep_end=');fprintf(fileID, '%d\n', length(t));
fprintf(fileID, '  %s', 'delta_step=');fprintf(fileID, '%d\n', 1);
fprintf(fileID, '  %s', 'NUM_THREADS=');fprintf(fileID, '%d\n', 12);
fprintf(fileID, '  %s', 'FWH_data_Format=');fprintf(fileID, '%d\n', 1);
fprintf(fileID, '%s\n', '$end'); 
fclose(fileID);
%% 写入FWH_Surface_Geo.dat
fileID = fopen('./FWH-dipole/FWH_Surface_Geo.dat', 'w'); % 打开文件，以写入模式打开  

fprintf(fileID, '%s\n', ' variables=x,y,z,n1,n2,n3,dS');  
fprintf(fileID, '   %d\n', 1); 
fprintf(fileID, '   %d %d %d\n',length(theta),length(phi),1);
d_theta=theta(2)-theta(1);
d_phi=phi(2)-phi(1);
for ii=1:length(theta)
    for j=1:length(phi)
    dS(ii,j)=r_FWH^2*sin(theta(ii))*d_phi*d_theta;
    fprintf(fileID, '%25.16E', x1(ii,j));  
    fprintf(fileID, '%25.16E', x2(ii,j));  
    fprintf(fileID, '%25.16E', x3(ii,j));  
    fprintf(fileID, '%25.16E', x1(ii,j)/r_FWH); %n1
    fprintf(fileID, '%25.16E', x2(ii,j)/r_FWH); %n2
    fprintf(fileID, '%25.16E', x3(ii,j)/r_FWH); %n3
    fprintf(fileID, '%25.16E\n', dS(ii,j));  
    end
end
test=sum(dS(:,:));
delta_S=sum(test)-4*pi*r_FWH^2
fclose(fileID);

%% 写入FWH-xxxxxxx.dat
for k=1:length(t)
    file_number = sprintf('%08d', k); 
    file_name = strcat('./FWH-dipole/FWH-',file_number, '.dat');
    fileID = fopen(file_name, 'w'); 
     for ii=1:length(theta)
        for j=1:length(phi)
        fprintf(fileID, '%25.16E',1+rou(ii,j,k)/rou0); %无量纲化
        fprintf(fileID, '%25.16E',(u(ii,j,k)+U0*cos(AoA))/U0); %叠加来流速度
        fprintf(fileID, '%25.16E',(v(ii,j,k)+U0*sin(AoA))/U0); %叠加来流速度
        fprintf(fileID, '%25.16E',w(ii,j,k)/U0);
        fprintf(fileID, '%25.16E\n',p(ii,j,k)/(rou0*U0^2));
        end
    end 
    fclose(fileID);
end
toc
% RMSP = zeros(length(theta), length(phi));
% for ii=1:length(theta)
%     for j=1:length(phi)
%     tmp(:)=p(ii,j,:);
%     RMSP(ii,j) = rms(tmp);
%     end
% end
% clear tmp
% set(groot,'defaultLineLineWidth',3.5);
% polarplot(phi,RMSP(1,:));

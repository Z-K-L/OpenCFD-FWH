clc;
clear;
tic;
%% Stationary monopole in a moving medium with AoA
format long
c0=340;
M0=0.6;
rou0=1;
AoA=45; AoA=AoA*pi/180;
M1=M0*cos(AoA);
M2=M0*sin(AoA);
A=1; %振幅 m^2/s
f=5; %频率
omiga=2*pi*f;
n_t=1e3;
t(:)=linspace(2.5/n_t,2.5,n_t);

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
set(groot,'defaultLineLineWidth',1);
% plot3(x1,x2,x3,'o')
% axis equal

dx1 = 1e-8;
dx2 = 1e-8;
dx3 = 1e-8;


R_s = zeros(length(theta), length(phi));
R = zeros(length(theta), length(phi));
R_spx1= zeros(length(theta), length(phi));
Rpx1= zeros(length(theta), length(phi));
R_smx1= zeros(length(theta), length(phi));
Rmx1= zeros(length(theta), length(phi));
R_spx2= zeros(length(theta), length(phi));
Rpx2= zeros(length(theta), length(phi));
R_smx2= zeros(length(theta), length(phi));
Rmx2= zeros(length(theta), length(phi));
R_spx3= zeros(length(theta), length(phi));
Rpx3= zeros(length(theta), length(phi));
R_smx3= zeros(length(theta), length(phi));
Rmx3= zeros(length(theta), length(phi));
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

R_spx3(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j))).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)).^2+(x3(ii,j)+dx3).^2));
Rpx3(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j))+R_spx3(ii,j))/(1-M0^2);

R_smx3(ii,j)= sqrt((M1*x1(ii,j)+M2*(x2(ii,j))).^2+(1-M0^2)*(x1(ii,j).^2+(x2(ii,j)).^2+(x3(ii,j)-dx3).^2));
Rmx3(ii,j)= (-M1*x1(ii,j)-M2*(x2(ii,j))+R_smx3(ii,j))/(1-M0^2);


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
    faipx3(ii,j,k)= A/(4*pi*R_spx3(ii,j))*exp(1i*omiga*(t(k)-Rpx3(ii,j)/c0));
    faimx3(ii,j,k)= A/(4*pi*R_smx3(ii,j))*exp(1i*omiga*(t(k)-Rmx3(ii,j)/c0));
    end
  end
end
clear R_s;clear R;
clear R_spx1;clear Rpx1
clear R_smx1;clear Rmx1;
clear R_spx2;clear Rpx2;
clear R_smx2;clear Rmx2;
clear R_spx3;clear Rpx3;
clear R_smx3;clear Rmx3;
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

pt = zeros(length(theta), length(phi),length(t));
for ii=1:length(theta)
  for j=1:length(phi)
    pt(ii,j,:)= -real(rou0*1i*omiga*fai(ii,j,:));
  end
end
clear tmp;
clear fai;

p=p1+pt;
clear p1;clear pt;
rou=p/c0^2;
%% 写入matlab.control
fileID = fopen('./FWH-monopole/matlab.control', 'w'); % 打开文件，以写入模式打开  
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
fileID = fopen('./FWH-monopole/Observers.dat', 'w'); % 打开文件，以写入模式打开  
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
fileID = fopen('./FWH-monopole/control.fwh', 'w'); % 打开文件，以写入模式打开  
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
fileID = fopen('./FWH-monopole/FWH_Surface_Geo.dat', 'w'); % 打开文件，以写入模式打开  

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
    file_name = strcat('./FWH-monopole/FWH-',file_number, '.dat');
    fileID = fopen(file_name, 'w'); 
     for ii=1:length(theta)
        for j=1:length(phi)
        fprintf(fileID, '%25.16E',1+rou(ii,j,k)/rou0);    %无量纲化
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

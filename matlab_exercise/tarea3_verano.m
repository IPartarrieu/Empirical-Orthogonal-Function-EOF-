clear all
load('data_tarea3_2022.mat')
%% EOF compuesta entre altura de Isoterma 0 y temperatura superficial
%% calculamos anomalías estandarizadas
for i=1:12
    iso_anom(:,:,i:12:size(ISO,3))=(ISO(:,:,i:12:end)-mean(ISO(:,:,i:12:end),3))./std(ISO(:,:,i:12:end),0,3);
    t2_anom(:,:,i:12:size(T2,3))=(T2(:,:,i:12:end)-mean(T2(:,:,i:12:end),3))./std(T2(:,:,i:12:end),0,3);
end
clear ISO T2
%% sacamos promedios para los meses de verano DJF
a=1;
for i=12:12:756-12
    iso_anom_ver(:,:,a)=mean(iso_anom(:,:,i:i+2),3);
    t2_anom_ver(:,:,a)=mean(t2_anom(:,:,i:i+2),3);
    a=a+1;
end
% tenemos que trabajar con los datos sin tendencia
iso_anom_ver_dt=detrend3(iso_anom_ver);
t2_anom_ver_dt=detrend3(t2_anom_ver);
%
% determinamos X Y y T para trabajar
[X,Y,T]=size(iso_anom_ver_dt);
%
%
%
%
%% EOF combinada, caso hacia abajo
S = reshape(permute(iso_anom_ver_dt,[3 1 2]),T,X*Y);
P = reshape(permute(t2_anom_ver_dt,[3 1 2]),T,X*Y);
F=[S';P'];
[L1,A1,E1,error]=EOF(F);

var=round((L1(1)/sum(L1))*100,2);
% significancia
n=62;
% % figure()
% % scatter(1:10,L1(1:10),'filled')
% % hold on
% % plot(L1(1:10)+L1(1:10)*sqrt(2/n),'+r','linewidth',2)
% % plot(L1(1:10)-L1(1:10)*sqrt(2/n),'+r','linewidth',2)
% % grid minor

%
% Graficamos A: Los componentes principales

figure()

subplot(311)
plot(1959:2020,-A1(:,1)/std(A1(:,1)),'linewidth',1.5,'Color','b')
grid minor
xlim([1959 2020])
ylim([-3 3])
title('1ra Componente Principal (CP)','fontsize',14)
ylabel('Z-Score','FontSize',14)
xlabel('Tiempo','FontSize',14)

corr_iso=corr3(iso_anom_ver_dt,-A1(:,1)/std(A1(:,1)));
corr_t2=corr3(t2_anom_ver_dt,-A1(:,1)/std(A1(:,1)));

subplot(312)
contourf(lon,lat,corr_iso')
hold on
borders('countries','center',180,'linewidth',1.5,'color','k')
colormap(redbluecmap)
colorbar
caxis([-1 1])
title('Correlación PTP entre ISO 0 y CP','fontsize',14)
ylabel('Latitud','FontSize',14)
xlabel('Longitud','fontsize',14)

subplot(313)
contourf(lon,lat,corr_t2')
hold on
borders('countries','center',180,'linewidth',1.5,'color','k')
colormap(redbluecmap)
colorbar
caxis([-1 1])
title('Correlación PTP entre T2 y CP','fontsize',14)
ylabel('Latitud','FontSize',14)
xlabel('Longitud','fontsize',14)
sgtitle(['Verano - P.V: ',num2str(var),'%'],'fontsize',20)
% se debería multiplicar por -1 para que se pueda apreciar mejor la figura

%% Caso "hacia al lado"
% clear F
% F=[S' P'];
% 
% 
% % En este caso, las variables no vienen acopladas, asi que las varianzas y
% % componentes principales se sumaran
% 
% L_iso=L2(1:62,1);
% L_t2=L2(63:end,1);
% A_iso=A2(1:62,1:62);
% A_t2=A2(63:end,63:end);
% 
% figure()
% %for i=1:5
% %    subplot(5,1,i)
%     plot(1959:2020,-A_iso(:,1)/std(A_iso(:,1)),'linewidth',1.5)
%     hold on
%     plot(1959:2020,-A_t2(:,1)/std(A_t2(:,1)),'LineWidth',1.5)
%     grid minor
%     legend('ISO','T2')
%     ylim([-5 5])
%     xlim([1959 2020])
% %end

%% sacar la correlación PTP
for i=1:61
    for j=1:49
        PTP_corr1(i,j)=corr(squeeze(iso_anom_ver_dt(i,j,:)),squeeze(t2_anom_ver_dt(i,j,:)));
    end
end
coef_det1=PTP_corr1.^2;

figure()
subplot(2,1,1)
contourf(lon,lat,PTP_corr1')
hold on
borders('countries','center',180,'linewidth',1.5,'color','k')
colormap("turbo")
colorbar
caxis([-1 1])
title('Series Originales','fontsize',14)
ylabel('Latitud','FontSize',14)
xlabel('Longitud','fontsize',14)

subplot(2,1,2)
contourf(lon,lat,PTP_corr2')
hold on
borders('countries','center',180,'linewidth',1.5,'color','k')
colormap("turbo")
colorbar
caxis([-1 1])
title('Series Reconstruidas','fontsize',14)
ylabel('Latitud','FontSize',14)
xlabel('Longitud','fontsize',14)
sgtitle('PTP Correlation','fontsize',20)

% varianza compartida entre ambas variables

%% Reconstrucción
% F=A*E'
F_r=A1(:,1:3)*E1(:,1:3)';
F_r=F_r';
S_r=F_r(1:length(F_r(:,1))/2,:);
P_r=F_r((length(F_r(:,1))/2)+1:end,:);

for i=1:62
    iso_r(:,:,i)=reshape(S_r(:,i),[61 49]);
    t2_r(:,:,i)=reshape(P_r(:,i),[61 49]);
end

for i=1:61
    for j=1:49
        PTP_corr2(i,j)=corr(squeeze(iso_r(i,j,:)),squeeze(t2_r(i,j,:)));
    end
end

coef_det2=PTP_corr2.^2;

figure()

subplot(2,1,1)
contourf(lon,lat,coef_det1')
hold on
borders('countries','center',180,'linewidth',1.5,'color','k')
colormap("turbo")
colorbar
caxis([0 1])
title('Series Originales','fontsize',14)
ylabel('Latitud','FontSize',14)
xlabel('Longitud','fontsize',14)

subplot(2,1,2)
contourf(lon,lat,coef_det2')
hold on
borders('countries','center',180,'linewidth',1.5,'color','k')
colormap("turbo")
colorbar
caxis([0 1])
title('Serie Reconstruida','fontsize',14)
ylabel('Latitud','FontSize',14)
xlabel('Longitud','fontsize',14)
sgtitle('Coeficiente de Determinación','fontsize',20)
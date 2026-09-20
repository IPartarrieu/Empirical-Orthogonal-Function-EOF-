clear all
load('data_tarea3_2022.mat')
%% EOF compuesta entre altura de Isoterma 0 y temperatura superficial
%% calculamos anomalías estandarizadas
for i=1:12
    iso_anom(:,:,i:12:size(ISO,3))=(ISO(:,:,i:12:end)-mean(ISO(:,:,i:12:end),3))./std(ISO(:,:,i:12:end),0,3);
    t2_anom(:,:,i:12:size(T2,3))=(T2(:,:,i:12:end)-mean(T2(:,:,i:12:end),3))./std(T2(:,:,i:12:end),0,3);
end
clear ISO T2
%% sacamos promedios para los meses de primavera SON
a=1;
for i=9:12:756
    iso_anom_pri(:,:,a)=mean(iso_anom(:,:,i:i+2),3);
    t2_anom_pri(:,:,a)=mean(t2_anom(:,:,i:i+2),3);
    a=a+1;
end
% tenemos que trabajar con los datos sin tendencia
iso_anom_pri_dt=detrend3(iso_anom_pri);
t2_anom_pri_dt=detrend3(t2_anom_pri);
%
% determinamos X Y y T para trabajar
[X,Y,T]=size(iso_anom_pri_dt);
%
%
%
%
%% EOF combinada, caso hacia abajo
S = reshape(permute(iso_anom_pri_dt,[3 1 2]),T,X*Y);
P = reshape(permute(t2_anom_pri_dt,[3 1 2]),T,X*Y);
F=[S';P'];
[L1,A1,E1,error]=EOF(F);
A1=-A1;
var=round((L1(1)/sum(L1))*100,2);

% significancia
n=63;
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
plot(1959:2021,A1(:,1)/std(A1(:,1)),'linewidth',1.5,'Color','b')
grid minor
xlim([1959 2021])
ylim([-3 3])
title('1ra Componente Principal (CP)','fontsize',14)
ylabel('Z-Score','FontSize',14)
xlabel('Tiempo','FontSize',14)

corr_iso=corr3(iso_anom_pri_dt,A1(:,1)/std(A1(:,1)));
corr_t2=corr3(t2_anom_pri_dt,A1(:,1)/std(A1(:,1)));

subplot(312)
contourf(lon,lat,corr_iso')
hold on
borders('countries','center',180,'linewidth',1.5,'color','k')
colormap(redbluecmap)
colorbar
caxis([-1 1])
title('Correlación PTP: ISO 0 y CP','fontsize',14)
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
sgtitle(['Primavera - P.V: ',num2str(var),'%'],'fontsize',20)

% se debería multiplicar por -1 para que se pueda apreciar mejor la figura

%% Caso "hacia al lado"
% clear F
% F=[S' P'];
% [L2,A2,E2,error2]=EOF(F);
% 
% % En este caso, las variables no vienen acopladas, asi que las varianzas y
% % componentes principales se sumaran
% 
% L_iso=L2(1:63,1);
% L_t2=L2(64:end,1);
% A_iso=A2(1:63,1:63);
% A_t2=A2(64:end,64:end);
% 
% figure()
% %for i=1:5
% %    subplot(5,1,i)
%     plot(1959:2021,A_iso(:,1)/std(A_iso(:,1)),'linewidth',1.5)
%     hold on
%     plot(1959:2021,A_t2(:,1)/std(A_t2(:,1)),'LineWidth',1.5)
%     grid minor
%     legend('ISO','T2')
%     ylim([-5 5])
%     xlim([1959 2021])
% %end
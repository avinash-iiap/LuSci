clear;
filename='compiled_seeing.csv';
fid=fopen(filename,'r');
filearray=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %s','Delimiter',',');
fclose(fid);
dn=datenum(filearray{1});
for i=2:6
    gsee(:,i-1)=filearray{i};
end
for i=7:10
    per_see(:,i-6)=filearray{i};
end
for i=11:14
    corr_coeff(:,i-10)=filearray{i};
end
filename=char(filearray{end});

[dn_unique,n]=unique(dn);
for i=1:length(dn_unique)
    dn_sort(i)=dn(n(i));
    gsee_sort(i,1:4)=gsee(n(i),1:4);
    per_see_sort(i,:)=per_see(n(i),:);
    corr_coeff_sort(i,:)=corr_coeff(n(i),:);
    filename_sort(i,:)=filename(n(i),:);
end
layers={'3 m','12 m','48 m','196 m','Total Seeing'};

xq=dn_sort(1):dn_sort(end);
vq1 = interp1(dn_sort,gsee_sort,xq,'pchip');
vq1(vq1<0)=0;
figure(1);
plot(xq,vq1,dn_sort,gsee_sort,'o');
box on; grid on;
datetick('x',20);
axis tight;
legend('3 m','12 m','48 m','196 m','Total Seeing');
xlabel('Date of observation');
ylabel('Median Seeing in arcseconds');
title('Atmospheric seeing at four pivot point altitudes');
print(gcf,'compiled_seeing','-dpsc2');

h = figure(2);
set(h, 'DefaultTextFontSize', 8);
gsee_3m=[];gsee_12m=[];gsee_48m=[];gsee_192m=[];gsee_total=[];grp=[];
for i=1:length(dn_sort)
    filename_see=strcat(filename_sort(i,:),'-seeing_data.txt');
    fid=fopen(filename_see);
    filearray=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f','Delimiter',' ');
    gsee_3m=vertcat(gsee_3m,filearray{end-4}(2:end));
    gsee_12m=vertcat(gsee_12m,filearray{end-3}(2:end));
    gsee_48m=vertcat(gsee_48m,filearray{end-2}(2:end));
    gsee_192m=vertcat(gsee_192m,filearray{end-1}(2:end));
    gsee_total=vertcat(gsee_total,filearray{end}(2:end));
    date_str=datestr(dn_sort(i));
    str_array=repmat(date_str,length(filearray{end}(2:end)),1);
    grp=vertcat(grp,str_array);
end
subplot(2,1,1);
boxplot(gsee_3m,grp,'labelorientation','inline','extrememode','compress','DataLim',[0,3]);
title(sprintf('Boxplot of atmospheric seeing variation at an altitude of %s ', layers{1}));
ylabel('Seeing in arcseconds');
subplot(2,1,2);
boxplot(gsee_12m,grp,'labelorientation','inline','extrememode','compress','DataLim',[0,3]);
title(sprintf('Boxplot of atmospheric seeing variation at an altitude of %s', layers{2}));
ylabel('Seeing in arcseconds');
print(gcf,'compiled_seeing','-dpsc2','-append');

h = figure(3);
set(h, 'DefaultTextFontSize', 8);
subplot(2,1,1);
boxplot(gsee_48m,grp,'labelorientation','inline','extrememode','compress','DataLim',[0,3]);
title(sprintf('Boxplot of atmospheric seeing variation at an altitude of %s', layers{3}));
ylabel('Seeing in arcseconds');
subplot(2,1,2);
boxplot(gsee_192m,grp,'labelorientation','inline','extrememode','compress','DataLim',[0,3]);
title(sprintf('Boxplot of atmospheric seeing variation at an altitude of %s', layers{4}));
ylabel('Seeing in arcseconds');
print(gcf,'compiled_seeing','-dpsc2','-append');

% h = figure(4);
% set(h, 'DefaultTextFontSize', 8);
% boxplot(gsee_total,grp,'labelorientation','inline','extrememode','compress','DataLim',[0,3]);
% title(strcat('Boxplot of total seeing over time'));
% ylabel('Seeing in arcseconds');
% print(gcf,'compiled_seeing','-dpsc2','-append');

figure(4);
for i=1:4
%     if i==5
%         subplot(3,2,[i i+1]);
%     else
        subplot(2,2,i);
%     end
    histfit(gsee(:,i));
%     if i==5
%         title(strcat('Integrated Seeing distribution, Median=',num2str(median(gsee(:,i))),' arcseconds'));
%     else
        title(sprintf('Seeing histogram at %s altitude, Median = %s%s', layers{i}, num2str(median(gsee(:,i))), char(34)));
%     end
    xlabel('Seeing in arcseconds');
end
print(gcf,'compiled_seeing','-dpsc2','-append');

% figure(6);
% for i=1:4 
%     corr_corr=corrcoef(corr_coeff(:,i),per_see(:,i));
%     corr=corr_corr(2);
%     subplot(2,2,i);
%     scatter(corr_coeff(:,i),per_see(:,i));
%     lsline
%     box on; grid on;
%     axis tight;
%     title(strcat('Correlation coefficient, and Seeing at ',layers(i),' as a percentage of the total seeing, R=',num2str(corr)),'FontSize',6);
%     xlabel('Correlation coefficient of layer seeing with total seeing (R)','FontSize',6);
%     ylabel('Percentage of layer seeing with respect to total seeing (%)','FontSize',6);
% end
% print(gcf,'compiled_seeing','-dpsc2','-append');
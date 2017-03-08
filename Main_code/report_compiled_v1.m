%---------------LuSci Compiled Report File---------------------------------
% After computing the covariances from the matlab file
% (voltage_to_covariance_500samp.m) and the turbulence profile from the IDL
% program (profrest.pro), this file computes the seeing values and makes a
% ps report consisting of the following:
%
% Page 1:
% 1. RAW Voltage data obtained on the specified date 
% 2. Turbulence profile at 'pivot points' (Default value of '3 m','12
% m','48 m','196 m') in the atmosphere (specified in the IDL program)
% 3. Integrated seeing values % till a specified altitude
%
% Page 2:
% Integrated Seeing values (with more clarity)
%
% Page 3:
% Mean Covariance Scatter Plot with respect to baselines of the photodiodes
% 
% Page 4:
% A 5 minute sample of the voltage time-series which shows the allowed and
% rejected voltages of channel 1. Histogram of the voltage distribution is
% also given.
%
% Page 5:
% The complete variance time-series which shows the allowed and
% rejected voltages of channel 1. Histogram of the variance distribution is
% also given.
%
% Page 6:
% The complete covariance time-series which shows the allowed and
% rejected voltages of the first baseline. Histogram of the covariance 
% distribution is also given.
%
% Page 7 and 8:
% The complete power spectrum of all channels with an x-axis range of half 
% the sampling frequency

clear;
nchan=6;
offset=0;range=5;   % For zoomed in analysis, offset is the start time in hours and range is the datapoints to be analyzed after the offset in minutes
offset_v=0; % Hours of covariance and variance data to be analyzed
chan=1;
cov_anal=1;
fpath=mfilename('fullpath');
filepath=fpath(1:length(mfilename('fullpath'))-length(mfilename));
addpath(filepath);
fclose('all');

if(exist('var_report.mat'))==2
    load('var_report.mat');
end
% delete('var_report.mat');
range_v=floor(length(cov_raw)/60);

c=size(tot);
time=zeros(c(1)*c(3),1);
for i=1:c(3)
    for j=1:c(1)
        time((c(1)*(i-1))+j)=tot(j,1,i);
    end
end
total=zeros(c(1)*c(3),1:nchan);
for i=1:c(3)
    for j=1:c(1)
        for k=1:nchan
            total((c(1)*(i-1))+j,k)=tot(j,k+1,i);
        end
    end
end
total_raw=zeros(length(total),nchan+1);
total_raw(:,1)=time;
total_raw(:,2:nchan+1)=total(:,1:nchan);
a=1;
index=[];
for i=1:length(time)
    if time(i)==0
        index(a)=i;
        a=a+1;
    end
end
time(index)=[];
total(index,:)=[];
time=(time./(24*60*60))+dn;
MJD_epoch='Nov 16, 1858,17:30';

d=tp_mjd_importfile(strcat(output_path,output_filename,'.tp'));
tptime=d+datenum(MJD_epoch);
[yy mm dd HH MM SS]=datevec(tptime(1,1));
if HH<7
    DateString = datestr(tptime(1,1) - 1,'dd mmmm yyyy');
else
    DateString = datestr(tptime(1,1),'dd mmmm yyyy');
end

% RAW Voltage Values
time=time(time<=max(tptime));
total=total(1:length(time),:);
time_start=tptime(1)-time(1);
time(1:ceil(time_start*43200000))=[];
total(1:ceil(time_start*43200000),:)=[];
figure(1);
subplot(3,1,1); 
plot(time,-(total));
datetick('x',15);
box on

axis tight;
[PATHSTR,NAME,EXT] = fileparts(input_file);
if(isempty(input_dark)~=1 && isempty(input_sky)~=1)
    title(sprintf('%s Dark and Sky Subtracted RAW voltage; Source File: %s',DateString,NAME),'Interpreter','none');
elseif(isempty(input_dark)~=1)
    title(sprintf('%s Dark subtracted RAW voltage; Source File: %s',DateString,NAME),'Interpreter','none');
else
    title(sprintf('%s RAW voltage; Source File: %s',DateString,NAME),'Interpreter','none');
end
ylabel('Voltage in Volts');

% Pivot point Cn^2 display
lusci_pivot=pivot_importfile(strcat(output_path,output_filename,'.tp'));
for i=1:5
    pivot(:,i)=10.^lusci_pivot(:,i);
end
subplot(3,1,2);

dl = diff([tptime pivot]); 
tol = min(dl(:,1))*2; 
jumpind = dl(:,1)>=tol;                  
blocks = cumsum(jumpind); 
p=size(pivot);
color_str=['r','g','b','c','m','y','k','w'];
for i=0:blocks(end)
    for j=1:p(2)
        pivot_temp=pivot(:,j);
        semilogy(tptime(blocks==i),pivot_temp(blocks==i),'.-','LineWidth',0.5,'Color',color_str(j));
        hold on;
    end
end
box on
datetick('x',15);
axis tight;
title('Pivot point C_n^2');
ylabel('C_n^2 values at pivot points');
legend('3 m','12 m','48 m','196 m','784 m');
hold off

cn2_int=gsee_importfile(strcat(output_path,output_filename,'.tp'));
for i=1:4
    gsee(:,i) = (cn2_int(:,i)./6.8e-13).^0.6;
end

% Integrated Seeing values
subplot(3,1,3);
g=size(gsee);
color_str=['r','g','b','c','k','y','k','w'];
for i=0:blocks(end)
    for j=1:g(2)
        gsee_temp=gsee(:,j);
        plot(tptime(blocks==i),gsee_temp(blocks==i),'.-','LineWidth',0.5,'Color',color_str(j));
        hold on;
    end
end
box on
datetick('x',15);
axis tight;
legend('3 m','12 m','48 m','196 m');
xlabel('Time of observation');
ylabel('Seeing in arcseconds');
title('Integrated Seeing values');
hold off;
print(gcf,strcat(output_path,output_filename),'-dpsc2');

figure(2);
for i=0:blocks(end)
    for j=1:g(2)
        gsee_temp=gsee(:,j);
        plot(tptime(blocks==i),gsee_temp(blocks==i),'.-','LineWidth',0.5,'Color',color_str(j));
        hold on;
    end
end
box on
datetick('x',15);
axis tight;
legend('3 m','12 m','48 m','196 m');
xlabel('Time of observation');
ylabel('Seeing in arcseconds');
title('Integrated Seeing values');
hold off;
print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');

% Covariance Scatter Plot
figure(3);
cov_import=cov_importfile(strcat(output_path,output_filename,'.cov'));
var_scatter=mean(mean(cov_import(:,2:nchan+1)));
detpos=[0,0.19,0.23,0.25,0.28,0.40];
a=0;
for i=2:nchan+1
    for j=i+1:nchan+1
        a=a+1;
        bsl(a)=detpos(j-1)-detpos(i-1);
        cov_scatter(a)=median(cov_import(:,a+nchan+1));
    end
end
scatter([0 bsl],[var_scatter cov_scatter],'fill');
box on
title('Covariance Scatter plot');
xlabel('Baseline in metres');
ylabel('Covariance');
print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');
saveas(gcf,strcat(output_path,output_filename,'-scatter'),'epsc');

% Sample voltage time series and histogram
figure(4);
start=offset*30000*60+1;
stop=offset*30000*60+range*30000;
total_raw=total_raw(start:stop,:);
a=0;b=0;
for i=1:length(total_raw)
    if total_raw(i,1)~=0
        a=a+1;
        total_valid(a,:)=total_raw(i,:);
    else
        b=b+1;
        total_rej(b,:)=total_raw(i,:);
        total_rej(b,1)=total_raw(2,1)+(i-2)*0.002;
    end
end
total_valid_time=total_valid(:,1);
total_rej_time=total_rej(:,1);
color_str=['r','g','b','c','k','y','k','w'];
subplot(2,2,[1,2]);
plot(total_valid_time,total_valid(:,chan+1),'.','LineWidth',0.5,'Color',color_str(1));
hold on;
plot(total_rej_time,total_rej(:,chan+1),'.','LineWidth',0.5,'Color',color_str(2));
axis tight;
legend('Valid Data','Rejected Data');
title(strcat(num2str(range),' minute chunk of channel',{' '},num2str(chan),' for',{' '},num2str(std_factor),' sigma outlier rejection'));
xlabel('Time in seconds');
ylabel('Voltage');
hold off;
subplot(2,2,3);
histfit(total_raw(:,chan+1));
title('RAW Voltage Distribution before rejection');
xlabel('Voltage');
ylabel('Histogram');
subplot(2,2,4);
histfit(total_valid(:,chan+1));
title(strcat('Voltage distribution after',{' '},num2str(std_factor),' sigma rejection'));
xlabel('Voltage');
ylabel('Histogram');
print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');

% Sample variance time series and histogram
a=0;b=0;
start=offset_v*60+1;
stop=(offset_v+range_v)*60;
var_raw=var_raw(start:stop,:);
cov_raw=cov_raw(start:stop,:);
s=size(var_raw);
for i=1:s(1)
    if var_raw(i,1)~=0
        a=a+1;
        var_valid(a,:)=var_raw(i,:);
    else
        b=b+1;
        var_rej(b,:)=var_raw(i,:);
        var_rej(b,1)=var_raw(2,1)+(i-2);
    end
end
a=0;b=0;
s=size(cov_raw);
for i=1:s(1)
    if cov_raw(i,1)~=0
        a=a+1;
        cov_valid(a,:)=cov_raw(i,:);
    else
        b=b+1;
        cov_rej(b,:)=cov_raw(i,:);
        cov_rej(b,1)=cov_raw(2,1)+(i-2);
    end
end
figure(5);
subplot(2,2,[1,2]);
plot(var_valid(:,1),var_valid(:,chan+1),'.','LineWidth',0.5,'Color',color_str(1));
hold on;
plot(var_rej(:,1),var_rej(:,chan+1),'.','LineWidth',0.5,'Color',color_str(2));
axis tight;
legend('Valid Data','Rejected Data');
title(strcat(num2str(range_v),' hour chunk of channel',{' '},num2str(chan),' variance for',{' '},num2str(std_factor),' sigma outlier rejection'));
xlabel('Time in minutes');
ylabel('Variance');
hold off;
subplot(2,2,3);
histfit(var_raw(:,chan+1));
title('RAW Variance Distribution before rejection');
xlabel('Variance');
ylabel('Histogram');
subplot(2,2,4);
histfit(var_valid(:,chan+1));
title(strcat('Variance distribution after',{' '},num2str(std_factor),' sigma rejection'));
xlabel('Variance');
ylabel('Histogram');
print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');

% Sample covariance time series and histogram
figure(6);
subplot(2,2,[1,2]);
plot(cov_valid(:,1),cov_valid(:,cov_anal+1),'.','LineWidth',0.5,'Color',color_str(1));
hold on;
plot(cov_rej(:,1),cov_rej(:,cov_anal+1),'.','LineWidth',0.5,'Color',color_str(2));
axis tight;
legend('Valid Data','Rejected Data');
title(strcat(num2str(range_v),' hour chunk of covariance number',{' '},num2str(cov_anal),' for',{' '},num2str(std_factor),' sigma outlier rejection'));
xlabel('Time in minutes');
ylabel('Covariance');
hold off;
subplot(2,2,3);
histfit(cov_raw(:,chan+1));
title('RAW Covariance Distribution before rejection');
xlabel('Covariance');
ylabel('Histogram');
subplot(2,2,4);
histfit(cov_valid(:,chan+1));
title(strcat('Covariance distribution after',{' '},num2str(std_factor),' sigma rejection'));
xlabel('Covariance');
ylabel('Histogram');
print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');

% Variance seeing values comparison (commented out because variance seeing
% is horribly inaccurate)
comp=intersect(d,cov_import(:,1));        
a=1;
for i=1:length(cov_import)
    if a<=length(comp)
        if comp(a)==cov_import(i,1)
            J(a)=median(cov_import(i,2:nchan+1))./2e5;
            a=a+1;
        end
    end
end
gsee_int = (J./6.8e-13).^0.6;
% figure(7); 
% color_str=['r','g','b','c','k','y','k','w'];
% for i=0:blocks(end)
%     for j=1:g(2)
%         gsee_temp=gsee(:,j);
%         plot(tptime(blocks==i),gsee_temp(blocks==i),'.-','LineWidth',0.5,'Color',color_str(j));
%         hold on;
%     end
% end
% % tptime_int=comp+datenum(MJD_epoch);
% % dl = diff([tptime_int (gsee_int)']); 
% % tol = min(dl(:,1))*2; 
% % jumpind = dl(:,1)>=tol;                  
% % blocks = cumsum(jumpind);
% % for i=0:blocks(end)
% %     plot(tptime_int(blocks==i),gsee_int(blocks==i),'.-','LineWidth',0.5,'Color',color_str(g(2)+1));
% %     legend('Total Seeing');
% %     hold on;
% % end
% box on
% legend('3 m','12 m','48 m','196 m');
% datetick('x',15);
% axis tight;
% xlabel('Time of observation');
% ylabel('Seeing in arcseconds');
% title('Integrated Seeing values (Black line represents total seeing computed from variance)');
% hold off;
% print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');

%Histogram of seeing at different altitudes
s=size(gsee);
fig=figure(7);
orient(fig,'portrait');
set(fig,'units','normalized','outerposition',[0 0 1 1])
layers={'3 m','12 m','48 m','196 m'};

for i=1:s(2)
    per_see(i)=(median(gsee(:,i))./median(gsee_int))*100;
    std_see(i)=std(gsee(:,i));
end
gsee_corr=corrcoef([gsee gsee_int']);
subplot(2,2,1);
histfit(gsee(:,1));
disp_str={sprintf('Median seeing = %0.3f arcseconds',median(gsee(:,1))),sprintf('SD in seeing = %0.3f arcseconds',std_see(1))};
annotation('textbox',[0.34 0.8 0.1 0.1],'String',disp_str,'FontSize',4);
title(strcat('Seeing distribution at',' ',layers(1)));
ylabel('Seeing in arcseconds');
subplot(2,2,2);
histfit(gsee(:,2));
disp_str={sprintf('Median seeing = %0.3f arcseconds',median(gsee(:,2))),sprintf('SD in seeing = %0.3f arcseconds',std_see(2))};
annotation('textbox',[0.785 0.8 0.1 0.1],'String',disp_str,'FontSize',4);
title(strcat('Seeing distribution at',' ',layers(2)));
ylabel('Seeing in arcseconds');
subplot(2,2,3);
histfit(gsee(:,3));
disp_str={sprintf('Median seeing = %0.3f arcseconds',median(gsee(:,3))),sprintf('SD in seeing = %0.3f arcseconds',std_see(3))};
annotation('textbox',[0.34 0.325 0.1 0.1],'String',disp_str,'FontSize',4);
title(strcat('Seeing distribution at',' ',layers(3)));
ylabel('Seeing in arcseconds');
subplot(2,2,4);
histfit(gsee(:,4));
disp_str={sprintf('Median seeing = %0.3f arcseconds',median(gsee(:,4))),sprintf('SD in seeing = %0.3f arcseconds',std_see(4))};
annotation('textbox',[0.785 0.325 0.1 0.1],'String',disp_str,'FontSize',4);
title(strcat('Seeing distribution at',' ',layers(4)));
ylabel('Seeing in arcseconds');
% subplot(3,2,[5 6]);
% histfit(gsee_int);
% title(strcat('Integrated Seeing distribution'));
% ylabel('Seeing in arcseconds');
% disp_str={sprintf('Median seeing = %0.3f arcsecond',median(gsee_int)),sprintf('SD in seeing = %0.3f arcsecond',std(gsee_int))};
% annotation('textbox',[0.775 0.2 0.1 0.1],'String',disp_str,'FontSize',4);
% print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');

% Power spectrum of 6 channels
sample=samplelvm_import(input_file);
nchan=6;
figure(8);
for i=1:nchan/2
    x=sample(:,i);
    subplot(3,1,i);
    [psdestx_tot,Fxx]=periodogram(x,rectwin(length(x)),length(x),500);
    plot(Fxx,10*log10(psdestx_tot)); grid on;
    axis tight;
    xlabel('Hz'); ylabel('Power/Frequency (dB/Hz)');
    title(strcat('Power Spectral Density estimate for channel no.',num2str(i)));
end
print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');
figure(10);
for i=(nchan/2)+1:nchan
    x=sample(:,i);
    subplot(3,1,i-nchan/2);
    [psdestx_tot,Fxx]=periodogram(x,rectwin(length(x)),length(x),500);
    plot(Fxx,10*log10(psdestx_tot)); grid on;
    axis tight;
    xlabel('Hz'); ylabel('Power/Frequency (dB/Hz)');
    title(strcat('Power Spectral Density estimate for channel no.',num2str(i)));
end
print(gcf,strcat(output_path,output_filename),'-dpsc2','-append');

% Appending layer seeing values and percentages in compilation file - for
% use in the campaign mode pipeline
if exist('compiled_seeing.txt','file')
    if isempty(strfind(fileread('compiled_seeing.txt'),output_filename))
        fid=fopen('compiled_seeing.txt','a');
        str=fprintf(fid,['\n',datestr(tptime(1),'dd-mmm-yyyy'),'\t',num2str(median(gsee)),'\t',num2str(median(gsee_int)),'\t',num2str(per_see),'\t',num2str(gsee_corr(end,1:end-1)),'\t',output_filename]);
        fclose('all');
    end
else
    fid=fopen('compiled_seeing.txt','w');
    if HH<7
        fprintf(fid,[datestr(tptime(1) - 1,'dd-mmm-yyyy'),'\t',num2str(median(gsee)),'\t',num2str(median(gsee_int)),'\t',num2str(per_see),'\t',num2str(gsee_corr(end,1:end-1)),'\t',output_filename]);
    else
        fprintf(fid,[datestr(tptime(1),'dd-mmm-yyyy'),'\t',num2str(median(gsee)),'\t',num2str(median(gsee_int)),'\t',num2str(per_see),'\t',num2str(gsee_corr(end,1:end-1)),'\t',output_filename]);
    end
    fclose('all');
end

% Heading datestring with Timestamp, AM, Cn2, seeing - for use in campaign
% mode
sd_filename=strcat(output_path,output_filename,'-seeing_data.txt');
fid=fopen(sd_filename,'w');
if HH<7
    fprintf(fid,[datestr(tptime(1,1)-1,'dd-mmm-yyyy'),'\n']);
else
    fprintf(fid,[datestr(tptime(1),'dd-mmm-yyyy'),'\n']);
end
for i=1:length(tptime)
    k=datestr(tptime(i));
    fprintf(fid,[k(13:length(k)),' ']);
    for j=1:5
        fprintf(fid,[num2str(pivot(i,j)),' ']);
    end
    for j=1:4
        fprintf(fid,[num2str(gsee(i,j)),' ']);
    end
    fprintf(fid,[num2str(gsee_int(i)),' ']);
    fprintf(fid,'\n');
end
fclose('all');
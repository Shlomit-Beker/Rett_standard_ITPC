%% plot TFR data. run after timeFreq_Rett.m
% Shlomit Beker shlomitbeker@gmail.com

figure;
if Group == 1
    suptitle(' NT - instantaneous phase')
else
    suptitle(' RETT instantaneous phase')
end
set(gca,'fontsize', 14);
hold on;

%%
for k = 1:6
% data is a matrix of frequency (rows) x time (columns)

TITLES = {'ISI 450 std','ISI 450 dev','ISI 900 std','ISI 900 dev','ISI 1800 std','ISI 1800 dev'};

clear Timeoi x y x2 y2 Data 

N = 500;
timeReduction= 183:410;
%timeReduction = [length(timeoi)./5*2:length(timeoi)./4*3];
Timeoi = timeoi(round(timeReduction));

Data = data{k}(:,round(timeReduction));  %the original data (and time epoch) had too many Nans in it for the interpolation function. 
[n, m] = size(Data);
[x,y] = meshgrid(Timeoi,freqoi); % low-res grid
[x2,y2] = meshgrid(Timeoi(1):1/N/5:Timeoi(end),freqoi(1):.01:freqoi(end));  %high-res grid
dataInterp = interp2(x,y,Data, x2,y2, 'linear'); %interpolate up

% if plotP == 1
%     dataInterp = interp2(x,y,pVal, x2,y2, 'linear'); %interpolate up
% end

presOrder = [1,4,2,5,3,6];

subplot(2,3,presOrder(k))

figure;
f = surf(x2,y2,dataInterp);
%f=surf(x,y,data); %no interpolation
f.EdgeColor = 'none';
f.FaceColor = 'interp';
f.FaceLighting = 'gouraud';
set(gca,'ydir','normal')
%set(gca,'YScale','log');

% ylabel('Frequency, Hz (log)')
% xlabel('Time, Sec.')
% title(TITLES{k});
box on
colormap jet; 
ax = gca;
caxis([0 0.3]); 
% 
% caxis(ax.CLim)
view(0,90)
yticks([2 5 10 15 20 25])
axis tight
hold on;
% stimuli lines

z = get(f,'ZData');
%set(f,'ZData',z-10);
z_max = max(max(get(f,'Zdata')));
%hold on;
% for k = 2:length(lines)
%     line([lines(k), lines(k)],[y2(1,1),y2(end,1)],[z_max,z_max]...
%         ,'Color','w','LineWidth',2,'LineStyle','--');
% end
end
%%
plotP = input(prompt);
promptK = 'choose K (1-6)';
if plotP ==1
figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,4,[1 4 2 8])
     h =  imagesc(Timeoi(1):1/N/5:Timeoi(end),freqoi(1):.01:freqoi(end),dataInterp)
    hold on; 
    dataInterp(isnan(dataInterp)==1)=1;
    pSig = dataInterp <0.05;
    [row,col,v] = find(pSig);   
    contour(Timeoi(1):1/N/5:Timeoi(end),freqoi(1):.01:freqoi(end),pSig,'k');
    colormap jet;
    colorbar;
    caxis([-0.1 0.1]);
    set(gca,'YDir','normal')
%     for k = 2:length(lines)
%         line([lines(k), lines(k)],[y2(1,1),y2(end,1)],[z_max,z_max]...
%         ,'Color','w','LineWidth',2,'LineStyle','--');
%     end
    ylabel('Frequency (Hz)');
    xlabel('Time (Sec.)');
    title('TD-ASD ITPC difference')
    set(gca,'fontsize', 14);
else
hold on; 
subplot(3,3,[1 4])


%plot(mean(dataInterp,2),freqoi)
freqVec = freqoi(1):.01:freqoi(end); 
%freqVec = freqoi(1):1/500/5:freqoi(end);
plot(mean(dataInterp,2),freqVec,'LineWidth',2);

axis tight
ylabel('Frequency, Hz')
xlabel('ITPC value')
set(gca,'fontsize', 12);
%%
hold on; 
dataInterp(dataInterp==1) = NaN;
timeVec = Timeoi(1):1/N/5:Timeoi(end);
ax1 = subplot(4,4,[9, 12]);
area(ax1,timeVec, nanmean(dataInterp,1),'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
box off

xlabel('Time, sec.')
ylabel('ITPC value')
set(gca,'fontsize', 12);
axis tight

ylim([0 0.02]);
set(gca,'fontsize', 14);

end
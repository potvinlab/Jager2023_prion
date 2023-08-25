%%Figure 2

%%
%Experiments used 
data_dir = '';

load([data_dir '/Ch_expt1.mat'])
load([data_dir '/Ch_expt2.mat'])
load([data_dir '/Ch_expt3.mat'])

%% pool experiments
Ch_expt2 = ensure_divTime_in_absolute_time(Ch_expt2);
Ch_expt3 = ensure_divTime_in_absolute_time(Ch_expt3);
ch_pool_start_1 = pool_experiments({  Ch_expt2 Ch_expt3 Ch_expt1},[1 1 1]);

%%
%Figure 2b - YFP concentration relative to loss event
data = ch_pool_start_1;
startTime = 1 ;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig] =prionLossDistribution(data,0,1, startTime);



fieldToAlign = 'YFPAvg';
lossTime = lostTimeNormal;
lossIndex = lostIndexNormal;
minLossTime = 40;
[alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

figure;
hold on
plot((-1000:1000)*8,nanmedian(alignedTraces(:,:),1),'o-', 'color',[0.7490 ,   0.6863  ,  0.8000]) %YFPavg

xlim([-10 1]*8)
ylim ([0 40])
xlabel('Time relative to loss event (min)')
ylabel(['Median YFP concentration (FU)'])
ylimits = ylim;
line([0 0],[ylimits(1) ylimits(2)], 'linewidth',2, 'color', [0.80,0.80,0.80],'linestyle','--')
box off
% Shaded std bars------------------------------------
% define boundaries of shaded error bars for current sample 
hold on
std_traces = nanstd(alignedTraces(:,:),1,1)./sqrt(sum(alignedTraces(:,:)>0,1));
mean_traces = nanmedian(alignedTraces(:,:),1);
std_traces(isnan(std_traces))=0;
mean_traces(isnan(mean_traces))=0;
x = (-1000:1000)*8;
lo = mean_traces - 2*std_traces;
hi = mean_traces + 2*std_traces;
% Draw patch of shaded error bars, set color, transparency
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)
disp(['Number of cells: ' num2str(max(nansum(alignedTraces(:,:)>0)))])
%%
%Figure 2c - Histogram prion loss in cell cycle
startTime = 1 ;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig] =prionLossDistribution(data,0,1, startTime);

neverIndex = neverIndexNormal;
lossIndex =  [lostIndexNormal lostIndexBig ];
lossTime =  [lostTimeNormal lostTimeBig ];
autofluo = 0; 
dataQ = prionPartErrorSignTot_keepOutliers(data,lossIndex,lossTime, neverIndex, autofluo);

fieldToAlign = 'cellPhase';
lossTime = lostTimeNormal;
lossIndex = lostIndexNormal;
minLossTime = 30;

% histogram of cell cycle loss events 
[alignedTraces] = alignFieldAtLossBeforeAfter(dataQ,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

indexCellsNotAtDivision = alignedTraces(:,1001) >0;
alignedCellCycle = alignedTraces;

figure; histogram(alignedTraces(:,1001),'BinWidth',.1,'facecolor', [0.7490 ,   0.6863  ,  0.8000])

xlabel({'Cell cycle position', 'at time of loss'})
ylabel('Number of cells')
disp(['Number of cells: ' num2str(sum(lossTime>30))])


%%
%Figure 2e - Absolute partitioning errors

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_at_div_ste] = align_trace_part_div_cen(data);


figure;plot(x_part, abs_part_error_in_div,'-o')

hold on;
mean_traces = abs_part_error_in_div;
x = x_part;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*abs_part_error_at_div_ste;
hi = mean_traces + 2*abs_part_error_at_div_ste;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)


xlim([-8 8])
ylim([0 0.2])
xlabel('Divisions relative to the loss')
ylabel('Absolute partitioning error')

ylimits = ylim;
xlimits = xlim;

h=line([0 0],[ylimits(1) ylimits(2)]);
set(h,'linestyle','--','color','k')

%%
%Figure 2f - Partitionng errors at cell division

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_at_div_ste] = align_trace_part_div_cen(data);

figure;plot(x_part_at_div, part_error_in_div_at_div,'-o')
hold on;
mean_traces = part_error_in_div_at_div;
x = x_part_at_div;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*part_error_in_div_at_div_ste;
hi = mean_traces + 2*part_error_in_div_at_div_ste;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)


plot(x_part_not_at_div, part_error_in_div_not_at_div,'-o')

mean_traces = part_error_in_div_not_at_div;
x = x_part_not_at_div;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*part_error_in_div_not_at_div_ste;
hi = mean_traces + 2*part_error_in_div_not_at_div_ste;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)


xlim([-5 5])
ylim([-0.2 0.2])
xlabel('Divisions relative to the loss')
ylabel('Partitioning error')

ylimits = ylim;
xlimits = xlim;

h=line([0 0],[ylimits(1) ylimits(2)]);
set(h,'linestyle','--','color','k')

h=line([xlimits(1) xlimits(2)],[0 0]);
set(h,'linestyle','--','color','k')
legend('Lost at cell division', 'Lost at diff. time of the cell cycle')
%%
%Figure 2g - Mean y positions of aggregates 
data=dataQ;
spotCen = 2;
lossTime = lostTimeNormal;
lossIndex = lostIndexNormal;
minLossTime = 30;


[alignedTraces] = alignSpotCenAtLossBeforeAfter(data,spotCen,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

%at least 10 points

time_min = (-1000:1000)*8 ;
validTimePoints = sum(~isnan(alignedTraces),1)>10 & time_min<0;
averageTrace = nanmean(alignedTraces(~indexCellsNotAtDivision,:),1);


figure;
% Using 0.1031746um/px, 
plot(time_min(validTimePoints), -averageTrace(validTimePoints)*0.1031746)

hold on;

mean_traces = -averageTrace(validTimePoints)*0.1031746;
x = time_min(validTimePoints);

ste = nanstd(alignedTraces(~indexCellsNotAtDivision,:),[],1)./sqrt(sum(~isnan(alignedTraces(~indexCellsNotAtDivision,:)),1));
ste = ste(validTimePoints)*0.1031746;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*ste;
hi = mean_traces + 2*ste;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)


validTimePoints = sum(~isnan(alignedTraces),1)>10 & time_min<0;
averageTrace = nanmean(alignedTraces(indexCellsNotAtDivision,:),1);
time_min = (-1000:1000)*8 ;
plot(time_min(validTimePoints), -averageTrace(validTimePoints)*0.1031746)


mean_traces = -averageTrace(validTimePoints)*0.1031746;
x = time_min(validTimePoints);

ste = nanstd(alignedTraces(indexCellsNotAtDivision,:),[],1)./sqrt(sum(~isnan(alignedTraces(indexCellsNotAtDivision,:)),1));
ste = ste(validTimePoints)*0.1031746;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*ste;
hi = mean_traces + 2*ste;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

 xlim([-400 50])
xlabel('Time relative to loss event (min)')
ylabel(['Mean y position of aggregates (\mum)'])
ylimits = ([-6 -2]);
line([0 0],[ylimits(1) ylimits(2)], 'linewidth',0.5, 'color', 'k','linestyle','--')
legend('Lost at cell division','Not lost at division','location','best')


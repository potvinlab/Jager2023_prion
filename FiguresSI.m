%%Supplemental figures

data_dir = '';

load([data_dir '/Ch_expt1.mat'])
load([data_dir '/Ch_expt1_CV.mat'])

load([data_dir '/Ch_expt2.mat'])
load([data_dir '/Ch_expt3.mat'])

load([data_dir '/Ml_expt3.mat'])
load([data_dir '/Ml_expt2.mat'])
load([data_dir '/Ml_expt1.mat'])

load([data_dir '/Lh_expt1.mat'])
Lh_expt1 = ensure_divTime_in_absolute_time(Lh_expt1);
load([data_dir '/Lh_expt2.mat'])
Lh_expt2 = ensure_divTime_in_absolute_time(Lh_expt2);
load([data_dir '/Lh_expt3.mat'])
Lh_expt3 = ensure_divTime_in_absolute_time(Lh_expt3);

load([data_dir '/Ch_ultra_R2.mat'])
load([data_dir '/Ch_ultra_R3.mat'])
load([data_dir '/Ch_ultra_R4.mat'])

load([data_dir '/Ch_expt1_no_curation.mat'])
load([data_dir '/Ch_expt2_no_curation.mat'])
load([data_dir '/Ch_expt3_no_curation.mat'])

load([data_dir '/Lh_expt4_no_curation.mat'])
load([data_dir '/Lh_expt5_no_curation.mat'])
load([data_dir '/Lh_expt6_no_curation.mat'])

%% pool experiments
Ch_expt2 = ensure_divTime_in_absolute_time(Ch_expt2);
Ch_expt3 = ensure_divTime_in_absolute_time(Ch_expt3);
ch_pool = pool_experiments({  Ch_expt2 Ch_expt3 Ch_expt1},[40 40 40]);
ch_pool_start_1 = pool_experiments({  Ch_expt2 Ch_expt3 Ch_expt1},[1 1 1]);

a6_pool = pool_experiments([Lh_expt1 Lh_expt2 Lh_expt3],[35 40 36]);
a6_pool_start1 = pool_experiments([Lh_expt1 Lh_expt2 Lh_expt3],[1 1 1]);

Ml_expt1 = ensure_divTime_in_absolute_time(Ml_expt1);
Ml_expt3 = ensure_divTime_in_absolute_time(Ml_expt3);
Ml_expt2 = ensure_divTime_in_absolute_time(Ml_expt2);

a9_pool = pool_experiments({Ml_expt1 Ml_expt3 Ml_expt2} ,[40 40 35]);
a9_pool_start_1 = pool_experiments({Ml_expt1 Ml_expt3 Ml_expt2} ,[1 1 1]);
%%
%Figure S2
%a - mean YFP (FU)
data = Ch_expt1;

[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,1);
fieldToAlign = 'YFPAvg';

lossIndex = [lostIndexNormal keptIndexNormal];

minLossTime = 1;

lossTime = ones(1,length(lossIndex))*50;

[alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));
figure;
plot(((-1000:1000)+50)*8/60,nanmedian(alignedTraces(:,:),1))

yfpTrace = nanmedian(alignedTraces(:,:),1)./nanmedian(nanmedian(alignedTraces(:,:),1));

xlim([0 15])
xlabel('Time (hours)')
ylabel(['Mean [YFP] (FU)'])

%b - Fraction of cells with prions, start time without equilibrium
data = Ch_expt1; 
startTime = 0;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1, startTime);

figure;plot(frame_index_small*8/60, loss_curve_small)
box off
xlim([0 1000/60])
ylim([0 1])
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
box off

%c - Fractions of cells with prions, start time with equilibrium 
data = Ch_expt1; 
startTime = 40;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1, startTime);

figure;plot(frame_index_small*8/60, loss_curve_small)
box off
xlim([0 1000/60])
ylim([0 1])
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
box off

%% Figure S3
%
%a - cell traces ordered by prion loss
data = Ch_expt1;
 startTime = 1;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig] =prionLossDistribution(data,0,1, startTime);

lostTime = [lostTimeNormal lostTimeBig keptTimeNormal keptTimeBig];
lostIndex = [lostIndexNormal lostIndexBig keptIndexNormal keptIndexBig];
plotTracesOrderedLoss(data,lostTime,lostIndex)

%b - individual cell trace
data = Ch_expt1;
lostTime = [lostTimeNormal lostTimeBig keptTimeNormal keptTimeBig];
lostIndex = [lostIndexNormal lostIndexBig keptIndexNormal keptIndexBig];
plotTracesOrderedLoss(data,lostTime,lostIndex)

ylim([900 900.5])
ylabel('Cell index')
set(gca,'yticklabel',{})
set(gca,'box','off')
set(gca,'box','off')
set(gcf,'PaperPositionMode','auto')
%c - Images of cell trace over time

%d - Number of spots added or lost between time points
%Ch SSB PrD lineage
data = Ch_expt1;
startTime = 40;
[dist] = spot_change_dist(data,startTime);
figure;histogram(dist)
xlim([-3 3])
xlabel({'Number spots added (or lost)', 'between time points'})
ylabel('Number of time points')

%Ch SSB PrD stable lineage
data = Ch_ultra_R2;
startTime = 40;
[dist] = spot_change_dist(data,startTime);
figure;histogram(dist)
xlim([-3 3])
xlabel({'Number spots added (or lost)', 'between time points'})
ylabel('Number of time points')

%% Figure S4
% 
%a - methods for classification of aggregates
dataCell= Ch_expt1_CV;
min_mov_hist = [];
peak_over_ave_hist = [];
cv_hist = [];
avg_dist_old = [];
avg_dist_small = [];
stdallnonzeros = @(x) std(x(x>0),0,'all');
meanallnonzeros = @(x) mean(x(x>0),'all');
for i=1:length(dataCell.spotNo)
    ind = find(dataCell.spotNo{i});
    %only keep traces that start at beginning for now
    if dataCell.frameNo{i}(1)<20   && median(dataCell.YFPAvg{i}(ind)) > 10
        
        
        % 3 out of the first 15 frames have spot, start in prion
        
        if sum(dataCell.spotNo{i}(5:20)>=1)>=3 
            
            % if tracking old pole aggregates, select cells that have spots
            % that move less than 3 pixels on average for 18 frames
            
            min_move = calc_min_spot_move(dataCell.spotCen{i});
            min_move = cumsum(min_move);
            min_move = min_move(19:end)-min_move(1:end-18);
            min_mov_hist = cat(2,min_mov_hist,min(min_move));
            
            
            [old_agg_pos, displacement] = old_pole_agg_position(dataCell.spotCen{i});
            if median(dataCell.YFPCV{i}(ind)) > 1.3
                avg_dist_old = cat(2,avg_dist_old, mean(displacement));
            else
                avg_dist_small = cat(2,avg_dist_small, mean(displacement));
            end
            peak_over_ave_hist = cat(2,peak_over_ave_hist,median(dataCell.YFPPeak{i}(ind)./dataCell.YFPAvg{i}(ind)));
            cv_hist = cat(2,cv_hist, median(dataCell.YFPCV{i}(ind)));
        end
    end
end
 min_mov_hist(min_mov_hist>300) = 300;
figure;histogram(min_mov_hist*0.1031746/18)
xlim([0 8*0.1031746])
xlabel({'Minimum aggregate displacement', 'averaged over 3 gen (\mum)'})
ylabel('Number of cells')
ylims =get(gca,'YLim');
h=line([2.5*0.1031746 2.5*0.1031746],[ylims(1) ylims(2)]);
set(h,'linestyle','--','color','k')


figure;histogram(peak_over_ave_hist)
xlim([0 10])
xlabel('Median peak/average fluorescence')
ylabel('Number of cells')
ylims=ylim;
h=line([5 5],[ylims(1) ylims(2)]);
set(h,'linestyle','--','color','k')

figure;histogram(cv_hist)
% xlim([0 2])
xlabel('Median CV of pixel values')
ylabel('Number of cells')
ylims=ylim;
h=line([1.3 1.3],[ylims(1) ylims(2)]);
set(h,'linestyle','--','color','k')

figure;histogram(avg_dist_old*0.1031746)
hold on;
histogram(avg_dist_small*0.1031746)
xlim([0 35*0.1031746])
xlabel({'Mean displacement of topmost','(old pole) aggregate (\mum / 8 min)'})
ylabel('Number of cells')
ylims=ylim;
legend('Old-pole aggregate cells','Small aggregates cells')



%% FIG S4b
%b - comparison of two methods 
%small aggregates 
data = Ch_expt1_CV;

startTime = 40;
figure;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
sum(lostTimeNormal>0) + length(keptTimeNormal)

plot(frame_index_small *8/60, loss_curve_small)

xlabel('Time (h)')
ylabel('Fraction of cells with prions')


[lostTime,keptTime,lostCell,lostIndex,keptCell,keptIndex,rejectedCell, rejectedIndex,neverCell,neverIndex,frame_index_small, loss_curve_small] =prionLossDistribution_peak(data,0,0,startTime);

sum(lostTime>0) + length(keptTime)

hold on;plot(frame_index_small *8/60, loss_curve_small)
xlim([0 15])

[lostTime,keptTime,lostCell,lostIndex,keptCell,keptIndex,rejectedCell, rejectedIndex,neverCell,neverIndex,frame_index_small, loss_curve_small] =prionLossDistribution_CV(data,0,0,startTime);

sum(lostTime>0) + length(keptTime)

hold on;plot(frame_index_small *8/60, loss_curve_small)
xlim([0 15])
legend('Length immobile','Peak over average', 'CV')

%big aggregates 
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1,startTime);
figure;
plot(frame_index_big *8/60, loss_curve_big)
hold on;
sum(lostTimeBig>0) + length(keptTimeBig)

xlabel('Time (h)')
ylabel('Fraction of cells with prions')



[lostTime,keptTime,lostCell,lostIndex,keptCell,keptIndex,rejectedCell, rejectedIndex,neverCell,neverIndex,frame_index_big, loss_curve_big] =prionLossDistribution_peak(data,0,1, startTime);
plot(frame_index_big * 8/60, loss_curve_big)
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
xlim([0 15])
ylim([0 1])
sum(lostTime>0) + length(keptTime)


[lostTime,keptTime,lostCell,lostIndex,keptCell,keptIndex,rejectedCell, rejectedIndex,neverCell,neverIndex,frame_index_big, loss_curve_big] =prionLossDistribution_CV(data,0,1, startTime);
plot(frame_index_big * 8/60, loss_curve_big)
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
xlim([0 15])
ylim([0 1])
sum(lostTime>0) + length(keptTime)
legend('Length immobile','Peak over average', 'CV')


%%
%Figure S5 - absolute partitioning errors 
data = Ch_expt1;
startTime = 1 ;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig] =prionLossDistribution(data,0,1, startTime);

neverIndex = neverIndexNormal;
lossIndex =  [lostIndexNormal lostIndexBig ];
lossTime =  [lostTimeNormal lostTimeBig ];
autofluo = 0; %
dataQ = prionPartErrorSignTot_keepOutliers(data,lossIndex,lossTime, neverIndex, autofluo);

data = dataQ;
fieldToAlign = 'frameNo';
lossIndex = neverIndexNormal;
lossTime = ones(1,length(lossIndex))*50;
minLossTime = 1;
[alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));
frameNoNever = alignedTraces;

fieldToAlign = 'Q';
[alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

growthRateNever = alignedTraces(frameNoNever>30);

fieldToAlign = 'frameNo';
lossIndex = lostIndexNormal;
lossTime = lostTimeNormal;
[alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

frameNoLost = alignedTraces;

fieldToAlign = 'Q';
[alignedTraces] = alignFieldAtLossBeforeAfter(data,fieldToAlign,lossTime(lossTime>minLossTime),lossIndex(lossTime>minLossTime));

delayAfterLoss = 20;

afterLost = alignedTraces(:,1001+delayAfterLoss:2000);
% Only use data where frameNo >30
afterLost = afterLost(frameNoLost(:,1001+delayAfterLoss:2000)>30);
beforeLost = alignedTraces(:,1:1000);
beforeLost = beforeLost(frameNoLost(:,1:1000)>30);

afterLost = afterLost(abs(afterLost)<1);
beforeLost = beforeLost(abs(beforeLost)<1);


figure; violinplot([ abs(afterLost); abs(beforeLost)], [ones(1,length(afterLost)) zeros(1,length(beforeLost))], 'Violincolor',[0.7490 ,   0.6863  ,  0.8000])
ylim([0 1])
ylabel('Absolute part. error')
set(gca,'Xticklabel',{'With prion','Without prion'})

%%
%Figure S6 - SSB orthologs
%a - fraction of cells with small aggregates 
data = ch_pool_start_1;
startTime = 1;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1, startTime);

 [cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    f_small_Ch] = align_trace_part_div_cen(data);

f_small_Ch_std_bootstrap =bootstrap_it(100,17,@align_trace_part_div_cen,data);

data = a9_pool_start_1;
startTime = 1;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1, startTime);

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    f_small_A9] = align_trace_part_div_cen(data);
 
f_small_A9_std_bootstrap =bootstrap_it(100,17,@align_trace_part_div_cen,data);

data = a6_pool_start1;
startTime = 1;
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1, startTime);

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    f_small_A6] = align_trace_part_div_cen(data);
 
f_small_A6_std_bootstrap =bootstrap_it(100,17,@align_trace_part_div_cen,data);

hold on

figure;barwitherr(2*[f_small_Ch_std_bootstrap, f_small_A9_std_bootstrap, f_small_A6_std_bootstrap],[f_small_Ch f_small_A9 f_small_A6])
set(gca,'XTickLabel',{'Ch','Ml','Lh'});
ylabel('Fraction of cells with small aggregates')

hold off
%%
%b - Partitioning errors 
data = a6_pool_start1;
[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste] = align_trace_part_div_cen(data);

figure;
plot(x_part,part_error_in_div);

hold on
mean_traces = part_error_in_div;
x = x_part;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*part_error_in_div_ste;
hi = mean_traces + 2*part_error_in_div_ste;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

data = a9_pool_start_1;
[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste] = align_trace_part_div_cen(data);

hold on
plot(x_part,part_error_in_div);
mean_traces = part_error_in_div;
x = x_part;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*part_error_in_div_ste;
hi = mean_traces + 2*part_error_in_div_ste;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

xlim([-4 4])
ylabel('Partitioning errors')
xlabel('Divisions relative to loss')

xlims= xlim;
ylims = ylim;

plot([0 0], ylims, 'k--')
plot(xlims, [0 0], 'k--')


legend('Lh','Ml')
%%
%Figure S7

%b - Ch SSB PrD Stable lineage
%Fraction of losses at cell division  %Fraction of cells with small aggregates 
data = Ch_ultra_R2;
 startTime = 40;

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div_ultra,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    f_small_ultra] = align_trace_part_div_cen(data);

[fraction_at_div_ultra_std_bootstrap, f_small_ultra_std_bootstrap] =bootstrap_it(100,[6 17],@align_trace_part_div_cen,data);

data = Ch_ultra_R3;
 startTime = 40;

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div_ultra,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    f_small_ultra_3] = align_trace_part_div_cen(data);

[fraction_at_div_ultra_std_bootstrap_3, f_small_ultra_std_bootstrap_3] =bootstrap_it(100,[6 17],@align_trace_part_div_cen,data);

data = Ch_ultra_R4;
 startTime = 40;

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div_ultra,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    f_small_ultra_4] = align_trace_part_div_cen(data);

[fraction_at_div_ultra_std_bootstrap_4, f_small_ultra_std_bootstrap_4] =bootstrap_it(100,[6 17],@align_trace_part_div_cen,data);

 data = Ch_expt1;

 startTime = 40;

[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste, ...
    x_part_at_div, part_error_in_div_at_div, part_error_in_div_at_div_ste, ...
    x_part_not_at_div, part_error_in_div_not_at_div, part_error_in_div_not_at_div_ste, ...
    abs_part_error_in_div, abs_part_error_in_div_ste, ...
    f_small] = align_trace_part_div_cen(data);

[fraction_at_div_std_bootstrap, f_small_std_bootstrap] =bootstrap_it(100,[6 17],@align_trace_part_div_cen,data);

figure;
barwitherr(2*[f_small_std_bootstrap f_small_ultra_std_bootstrap f_small_ultra_std_bootstrap_3 f_small_ultra_std_bootstrap_4],[f_small f_small_ultra f_small_ultra_3 f_small_ultra_4])
set(gca,'XTickLabel',{'Typical lineage','Ultra stable lineage R2','Ultra stable lineage R3','Ultra stable lineage R4'},'XTickLabelRotation',45);
ylabel({'Fraction of cells', 'with small aggregates'})
ylim([0 1])
figure;
barwitherr(2*[fraction_at_div_std_bootstrap fraction_at_div_ultra_std_bootstrap],[fraction_at_div fraction_at_div_ultra])
set(gca,'XTickLabel',{'Typical lineage','Ultra stable lineage'},'XTickLabelRotation',45);
ylabel({'Fraction of losses', 'at cell division'})
ylim([0 1])



%%
%Figure S8
% a - Doubling time 
[doubling_time_no_prion, doubling_time_with_prion]= calc_growth_rate_prion({  Ch_expt2 Ch_expt3 Ch_expt1});
[doubling_time_no_prion_a6, doubling_time_with_prion_a6]= calc_growth_rate_prion({  Lh_expt1 Lh_expt2 Lh_expt3});

doubling_times_no_prion_cat = [doubling_time_no_prion; doubling_time_no_prion_a6];
prion_norm = {doubling_times_no_prion_cat/mean(doubling_times_no_prion_cat) doubling_time_with_prion./doubling_time_no_prion doubling_time_with_prion_a6./doubling_time_no_prion_a6};
ste_2 = cellfun(@std,prion_norm)./sqrt([6 3 3])*2;
figure;barwitherr(ste_2,cellfun(@mean,prion_norm))
set(gca,'xticklabel', {'No prion', 'Ch', 'Lh'})
ylabel('Normalized doubling time')

%%

ch_pool_full = pool_experiments({Ch_expt1_no_curation Ch_expt2_no_curation Ch_expt3_no_curation}, [1 1 1]);          
a6_pool_full = pool_experiments({Lh_expt4_no_curation Lh_expt5_no_curation Lh_expt6_no_curation}, [1 1 1]);          
%%
%b - death rate 
data = ch_pool_full;

data = identify_death_time(data);

[death_rate_never, death_rate_combined, death_rate_never_dist, death_rate_combined_dist] = calculate_death_rate_prion_combined(data, 1);

data = a6_pool_full;
data = identify_death_time(data);

[~, death_rate_combined_a6, ~, death_rate_combined_dist_a6] = calculate_death_rate_prion_combined(data, 1);

figure;
hold on
barwitherr([std(death_rate_never_dist) std(death_rate_combined_dist) std(death_rate_combined_dist_a6)]/8*60,[1 2 3],[death_rate_never death_rate_combined death_rate_combined_a6]/8*60 )

axis = gca;
axis.XTick =[ 1 2 3];
axis.XTickLabel = {'No prion','Ch','Lh'};
ylabel('Death rate (h^{-1})')
axis.XTickLabelRotation = 45;
axis.Box = 0;





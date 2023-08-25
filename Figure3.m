%%Figure 3

%%
%Experiments used

data_dir = '';

load([data_dir '/Ch_expt1.mat'])
load([data_dir '/Ch_expt2.mat'])
load([data_dir '/Ch_expt3.mat'])
%%
load([data_dir '/Ml_expt3.mat'])
load([data_dir '/Ml_expt2.mat'])
load([data_dir '/Ml_expt1.mat'])

load([data_dir '/Lh_expt1.mat'])
Lh_expt1 = ensure_divTime_in_absolute_time(Lh_expt1);
load([data_dir '/Lh_expt2.mat'])
Lh_expt2 = ensure_divTime_in_absolute_time(Lh_expt2);
load([data_dir '/Lh_expt3.mat'])
Lh_expt3 = ensure_divTime_in_absolute_time(Lh_expt3);
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
%Figure 3a - prion loss curves SSB orthologs
data = a6_pool;

figure;
startTime = 1;
[lostTime,~,~,~,~,~,~, ~,~,~,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
n = sum(lostTime>0)
plot(frame_index_small*8/60, loss_curve_small)
showfit('2^(-x/tau)+b')
hold on;

std_bootstrap =bootstrap_it(100,12,@prionLossDistribution,data,0,0,startTime);
mean_traces = loss_curve_small(1:length(std_bootstrap));
x = frame_index_small*8/60;
x=x(1:length(std_bootstrap));
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - std_bootstrap;
hi = mean_traces + std_bootstrap;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

data = a9_pool;

startTime = 1;
[lostTime,~,~,~,~,~,~, ~,~,~,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
n = sum(lostTime>0)

plot(frame_index_small*8/60, loss_curve_small)
showfit('2^(-x/tau)+b')
std_bootstrap =bootstrap_it(100,12,@prionLossDistribution,data,0,0,startTime);
mean_traces = loss_curve_small(1:length(std_bootstrap));
x = frame_index_small*8/60;
x=x(1:length(std_bootstrap));
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*std_bootstrap;
hi = mean_traces + 2*std_bootstrap;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)


data = ch_pool; 
startTime = 1;
[lostTime,~,~,~,~,~,~, ~,~,~,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
n = sum(lostTime>0)

plot(frame_index_small*8/60, loss_curve_small)

std_bootstrap =bootstrap_it(100,12,@prionLossDistribution,data,0,0,startTime);
mean_traces = loss_curve_small(1:length(std_bootstrap));
x = frame_index_small*8/60;
x=x(1:length(std_bootstrap));
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*std_bootstrap;
hi = mean_traces + 2*std_bootstrap;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

hold off


box off
xlim([0 15])
ylim([0 1])
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
box off

%%
%Figure 3d - Fraction of losses at cell division 
data = ch_pool_start_1;
[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_Ch] = align_trace_part_div_cen(data);
std_bootstrap_Ch =bootstrap_it(100,6,@align_trace_part_div_cen,data);

data = a6_pool_start1;
[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_A6] = align_trace_part_div_cen(data);

std_bootstrap_A6 =bootstrap_it(100,6,@align_trace_part_div_cen,data);

data = a9_pool_start_1;
[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_A9] = align_trace_part_div_cen(data);

std_bootstrap_A9 =bootstrap_it(100,6,@align_trace_part_div_cen,data);

%%
figure;barwitherr(2*[std_bootstrap_Ch std_bootstrap_A6, std_bootstrap_A9],[fraction_Ch fraction_A6 fraction_A9])
set(gca,'XTickLabel',{'Ch','Lh', 'Ml'})
ylabel('Fraction of losses at cell division')
xlims=xlim;
hold on;
plot(xlims,[1/Ch_expt1.avgDivDuration 1/Ch_expt1.avgDivDuration], '--')
%%
%Figure 3e - Mean y position of aggregates for orthologs
data = a6_pool_start1;
[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste] = align_trace_part_div_cen(data);

figure;

plot(time_min_part,-spot_pos*0.1031746);

hold on
mean_traces = -spot_pos*0.1031746;
x = time_min_part;
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*spot_pos_ste*(0.1031746);
hi = mean_traces + 2*spot_pos_ste*(0.1031746);
valid_indices = ~isnan(spot_pos_ste) & time_min_part<0;
x=x(valid_indices);
lo=lo(valid_indices);
hi=hi(valid_indices);
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)


data = a9_pool_start_1;
[cell_cycle_at_loss,x_part, part_error_in_div, time_min_part, spot_pos,fraction_at_div,part_error_in_div_ste,spot_pos_ste] = align_trace_part_div_cen(data);

hold on
plot(time_min_part,-spot_pos*0.1031746);


hold on
mean_traces = -spot_pos*0.1031746;
valid_indices = ~isnan(spot_pos_ste) & time_min_part<0;
x = time_min_part;
x=x(valid_indices);
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*spot_pos_ste*(0.1031746);
hi = mean_traces + 2*spot_pos_ste*(0.1031746);
lo=lo(valid_indices);
hi=hi(valid_indices);
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

xlim([-100 10])
ylabel('Mean y position of aggregates (\mum)')
xlabel('Time relative to loss (min)')
legend('Lh','Ml')


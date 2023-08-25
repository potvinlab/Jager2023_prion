%%Figure 1

%%
%%Experiments used 

data_dir = '';

load([data_dir '/Ch_expt1.mat'])
load([data_dir '/Ch_expt2.mat'])
load([data_dir '/Ch_expt3.mat'])
%%
%Figure 1d - Kymographs Ch SSB 

%% pool experiments

ch_pool= pool_experiments({  Ch_expt2 Ch_expt3 Ch_expt1 },[40 40 40]);

%%
% Figure 1e - Prion loss curve Ch SSB
data = ch_pool;
startTime = 1;

[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal, frame_index, loss_curve] =prionLossDistribution_mixed(data,0,startTime);
figure;
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,frame_index*8/60, loss_curve)
xlim([0 1000/60])
ylim([0 1])
hold on;

std_bootstrap =bootstrap_it(100,12,@prionLossDistribution_mixed,data,0,startTime);

mean_traces = loss_curve(1:length(std_bootstrap));
x = frame_index*8/60;
x=x(1:length(std_bootstrap));
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*std_bootstrap;
hi = mean_traces + 2*std_bootstrap;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

xlabel('Time (h)')
ylabel('Fraction of cells with prions')

disp(["Number of cells: " num2str(sum(lostTimeNormal>0) + sum(keptTimeNormal>0) )])

ax2 = axes(t);
p= plot(ax2, frame_index/Ch_expt1.avgDivDuration, loss_curve);
axes(ax2)
xlim([0 1000/60/Ch_expt1.avgDivDuration/8*60])
ylim([0 1])
xlabel('Time (generation)')

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';

ax1.Box = 'off';
ax2.Box = 'off';
p.LineStyle = 'none';


ax2.YTickLabel = {};
axes(ax2);
h = gca();
h.XAxis.Label.Color=[0 0 0];
h.XAxis.Label.Visible='on';


%%
%Figure 1f - Small aggregates loss curve

[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1, startTime);

std_bootstrap =bootstrap_it(100,12,@prionLossDistribution,data,0,0,startTime);

figure;plot(frame_index_small*8/60, loss_curve_small)
hold on;
mean_traces = loss_curve_small(1:length(std_bootstrap));
x = frame_index_small*8/60;
x=x(1:length(std_bootstrap));
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*std_bootstrap;
hi = mean_traces + 2*std_bootstrap;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

xlim([0 1000/60])
ylim([0 1])
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
box off
showfit('2^(-x/tau)+b')
disp(["Number of cells: " num2str(sum(lostTimeNormal>0) + sum(keptTimeNormal>0) )])
%%
%Figure 1g - Old pole aggregates loss curve
 
[lostTimeNormal,keptTimeNormal,~,lostIndexNormal,~,keptIndexNormal,~, ~,~,neverIndexNormal,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
[lostTimeBig,keptTimeBig,~,lostIndexBig,~,keptIndexBig,~, ~,~,neverIndexBig,frame_index_big, loss_curve_big] =prionLossDistribution(data,0,1, startTime);

figure;plot(frame_index_big*8/60, loss_curve_big)

std_bootstrap =bootstrap_it(100,12,@prionLossDistribution,data,0,1,startTime);

box off
hold on;
mean_traces = loss_curve_big(1:length(std_bootstrap));
x = frame_index_big*8/60;
x=x(1:length(std_bootstrap));
% Draw patch of shaded error bars, set color, transparency
lo = mean_traces - 2*std_bootstrap;
hi = mean_traces + 2*std_bootstrap;
error = fill([x,fliplr(x)],[lo,fliplr(hi)],[0.7490 ,   0.6863  ,  0.8000],'HandleVisibility','off');
set(error,'edgecolor','none','facealpha',0.2)

xlim([0 1000/60])
ylim([0 1])
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
box off
%%
disp(["Number of cells: " num2str(sum(lostTimeBig>0) + sum(keptTimeBig>0) )])
disp(["Average rate of loss (h^-1): " num2str(calc_loss_rate_per_frame(ch_pool,1,1)*60/8)])
disp(["1/Average rate of loss (h^-1): " num2str(1/calc_loss_rate_per_frame(ch_pool,1,1)*8/60)])


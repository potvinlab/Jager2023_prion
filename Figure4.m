%% Figure 4 
%%
%Experiments used
data_dir = '';
load([data_dir '/Ch_expt1.mat'])

load([data_dir '/Lh_colA_expt1.mat'])
load([data_dir '/Lh_colA_expt2.mat'])
load([data_dir '/Lh_colA_expt3.mat'])
load([data_dir '/Lh_colA_expt4.mat'])

load([data_dir '/Ch_ultra_R2.mat'])
load([data_dir '/Ch_ultra_R3.mat'])
load([data_dir '/Ch_ultra_R4.mat'])

%%
%Figure 4a - Same Lh colony
startTime = 32;
[~,~,~,~,~,~,~, ~,~,~,frame_index_small_1,loss_curve_small_1] =prionLossDistribution(Lh_colA_expt1,0,0,startTime);
startTime = 42;
[~,~,~,~,~,~,~, ~,~,~,frame_index_small_2,loss_curve_small_2] =prionLossDistribution(Lh_colA_expt3,0,0,startTime);
startTime = 26;
[~,~,~,~,~,~,~, ~,~,~,frame_index_small_3,loss_curve_small_3] =prionLossDistribution(Lh_colA_expt4,0,0,startTime);
startTime = 40;
[~,~,~,~,~,~,~, ~,~,~,frame_index_small_4,loss_curve_small_4] =prionLossDistribution(Lh_colA_expt2,0,0,startTime);

frame_index = 0:(min([frame_index_small_1(end) frame_index_small_2(end) frame_index_small_3(end) frame_index_small_4(end)]));
length_trace = length(frame_index);

mean_loss_curve1 = mean([loss_curve_small_1(1:length_trace); ...
                        loss_curve_small_2(1:length_trace); ...
                        loss_curve_small_3(1:length_trace); ...
                        loss_curve_small_4(1:length_trace)], 1);


frame_index = frame_index *8/60;

figure;p_mean1 = plot(frame_index,mean_loss_curve1,'-'); 

hold on
p_ind = plot(frame_index, loss_curve_small_1(1:length_trace),'linewidth',0.5,'Color',[0 0.4470 0.7410]);
plot(frame_index, loss_curve_small_2(1:length_trace),'linewidth',0.5,'Color',[0 0.4470 0.7410])
plot(frame_index, loss_curve_small_3(1:length_trace),'linewidth',0.5,'Color',[0 0.4470 0.7410])
plot(frame_index, loss_curve_small_4(1:length_trace),'linewidth',0.5,'Color',[0 0.4470 0.7410])

xlabel('Time (h)')
ylabel('Fraction of cells in prion state')
xlim([0 10])
legend([p_mean1 p_ind],{'Mean of colony A','Diff. experiments A'}, 'location', 'best')
%% ultra stable loss rate
ultra_pool = pool_experiments({Ch_ultra_R2, Ch_ultra_R3, Ch_ultra_R4 }, [32 25 25]);
1/(calc_loss_rate_per_frame(ultra_pool,1,0)/8*60)
%%
%Figure 4b - Stable Ch SSB lineage over time
data = Ch_ultra_R2;
startTime = 32;
[lostTime,~,~,~,~,~,~, ~,~,~,frame_index_small,loss_curve_small] =prionLossDistribution(data,0,0,startTime);
n = sum(lostTime>0)
figure;plot(frame_index_small*8/60, loss_curve_small)
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




data = Ch_ultra_R3;
startTime = 25;
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

data = Ch_ultra_R4;
startTime = 25;
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

data = Ch_expt1;
startTime = 40;
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

box off
xlim([0 10])
ylim([0 1])
xlabel('Time (h)')
ylabel('Fraction of cells with prions')
box off
legend('Round 1', 'Round 2','Round 3','Diff. lineage')


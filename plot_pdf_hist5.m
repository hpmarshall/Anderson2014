% plot_pdf_hist
% HPM 03/09/05
% this function plots both a histogram and nonparametric pdf,
% now with two separate y-axes (07/10/05)
% INPUT: x=data sample
%     [binx]=[50 equal-spaced] x-locations of histogram bins
%        [h]=[2*binsize] smoothing parameter, 
%     [xval]=[100 equal-spaced] locations to be sampled
% OUTPUT: xval = locations density was calculated
%         density = density at these locations
% SNTX: [xval,density,N,xout] = plot_pdf_hist5(x,binx,h,xval)


function [xval,density,N,xout] = plot_pdf_hist5(x,binx,h,xval)

if nargin<2
    [N,xout]=hist(x,50);
    wide=mean(diff(xout))*2;
    tdx=max(xout)-min(xout);x1=min(xout)-0.1*tdx; 
    x2=max(xout)+0.1*tdx; dx=(x2-x1)/300; % default evaluate at 300 points
    xi=x1:dx:x2;
    [density,xval]=ksdensity(x,xi,'width',wide);
elseif nargin<3
   [N,xout]=hist(x,binx);
    wide=mean(diff(xout))*2;
    tdx=max(xout)-min(xout);x1=min(xout)-0.1*tdx; 
    x2=max(xout)+0.1*tdx; dx=(x2-x1)/300; % default evaluate at 300 points
    xi=x1:dx:x2;
    [density,xval]=ksdensity(x,xi,'width',wide);
elseif nargin<4
    [density,xval]=ksdensity(x,'width',h);
    [N,xout]=hist(x,binx);
    wide=mean(diff(xout))*2;
    tdx=max(xout)-min(xout);x1=min(xout)-0.1*tdx; 
    x2=max(xout)+0.1*tdx; dx=(x2-x1)/300; % default evaluate at 300 points
    xi=x1:dx:x2;
    [density,xval]=ksdensity(x,'width',h);
else
    [density,xval]=ksdensity(x,xval,'width',h);
    [N,xout]=hist(x,binx);
end

% first make boxplot:
htemp=boxplot(x,'orientation','horizontal','position',1.1*max(N),'width',0.15*max(N),'notch','on');
set(htemp(isfinite(htemp)),'LineWidth',2)
%htemp=boxplot(x,'orientation','horizontal','position',1.1*max(N),'width',mean(diff(xout)),'notch','on');
h(1)=htemp(1); hold on
% first plot histogram
h(2)=bar(xout,N,1); hold on
% now lets scale the pdf and plot on top:
density=density*max(N)/max(density);
h(3)=plot(xval,density,'r','LineWidth',3);
axis([min(xval) max(xval) 0 1.2*max(N)]);
h(4)=plot([mean(x) mean(x)], [0 max(N)],'g','LineWidth',3);
%hl=legend(h,'boxplot','histogram','non-par. pdf','mean');
%set(hl,'FontSize',10); set(hl,'Location',[)
set(gca,'YTick',[(0:0.1*max(N):1.2*max(N))]);
set(gca,'YTickLabel',num2str((0:0.1*max(N):1.2*max(N))'));
set(gca,'YTickMode','auto'); % make the labeling auto again
set(gca,'YTickLabelMode','auto'); %

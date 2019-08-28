% DryCreekVariogram_hpm
% this code is an edited version of Brian Anderson's thesis work on
% variograms at Dry Creek
% HPM 08/14/11

Fs=16; Fts=20;

%clear all;close all
MagnaFiles={'TL12_23_10new.txt','TL1_7_10.txt','TL1_19_10.txt','TL2_10_10.txt','TL3_3_10.txt','316reformat.txt',...
    'LDP1_7_10.txt','LDP1_22_10.txt','LDP2_5_10.txt','LDP2_26_10.txt','LDP_3_12_10.txt','LDP4_1_10.txt',...
    'LDP4_15_10.txt','UDC1_8_10.txt','UDC1_29_10.txt','UDC_3_5_10.txt','UDC4_14_10.txt','UDC5_2_10.txt'}; % load datafiles
Date={'12/23/09','1/7/10','1/19/10','2/10/10','3/3/10','3/16/10','1/7/10','1/22/10','2/5/10','2/26/10',...
    '3/12/10','4/1/10','4/15/10','1/8/10','1/29/10','3/5/10','4/14/10','5/2/10'};
nbins=50; % number of equally spaced variogram bins
nAbin=1; % omnidirectional variograms
figure(1);clf
col=[3 6 6 6 6 3 3 3 3 3 7 7 7 6 6 6 7 6];
x0=569142;
y0=4842241;
G={[60 100 25],[50 30 20],[120 40 60],[400 40 50],[400 200 50],[550 100 50],[160 65 20],[300 80 40],[400 60 50],...
    [400 30 50],[400 30 45],[500 20 55],[400 30 60],[150 20 60],[500 200 100],[350 10 60],[600 30 90],[600 20 60]}; % inital guess of sill, range, nugget, from visual observation of data
% calculating semivariograms and fitting models
for n=1:length(MagnaFiles)
    T1=load(MagnaFiles{n});
    depth{n}=T1(:,col(n));
    if strcmp(MagnaFiles{n},'TL12_23_10new.txt') % fix the non-MagnaProbed coordinates
        x1=T1(1,2); y1=T1(1,1); % get starting point
        x2=T1(end,2); y2=T1(end,1); % get ending point
        Easting{n}=(x1+linspace(0,(x2-x1),299))'; % assume equally spaced points
        Northing{n}=(y1+linspace(0,(y2-y1),299))'; % assume equally spaced points
    else
        Easting{n}=T1(:,2);
        Northing{n}=T1(:,1);
    end
    % use only points along 300m traverse at treeline
    if strcmp(MagnaFiles{n}(1),'T')
        I2=find(diff(Easting{n})<-50);
        if isempty(I2)
            I2=length(Easting{n});
        end
        Easting{n}=Easting{n}(1:I2);
        Northing{n}=Northing{n}(1:I2);
        n
        depth{n}=depth{n}(1:I2);
    end
    % use only coincident points at LDP
    if strcmp(MagnaFiles{n}(1),'L')
        % get points we want to remove that were only measured in 1 survey
        I3=(Northing{n}>(4842973+74) & Easting{n}<(570691+20)) | Northing{n}<(4842973+14);
        I4=find(~I3); % get only those points that don't meet the above conditions.
        Easting{n}=Easting{n}(I4); % remove the bad points
        Northing{n}=Northing{n}(I4); % remove the bad points
        depth{n}=depth{n}(I4);
    end
    % use only nearly coincident points at UDC
    if strcmp(MagnaFiles{n}(1),'U')
        x0=572412.5; y0=4844464.5;
        I5=(Northing{n}>(y0+289) & Northing{n}<(y0+317) & Easting{n}<(x0+172));
        I6=(Northing{n}>(y0+96) & Northing{n}<(y0+101) & Easting{n}>(x0+184));
        I7=(Northing{n}>(y0+126.5) & Northing{n}<(y0+127.5) & Easting{n}>(x0-5) & Easting{n}<(x0+5));
        I8=(Northing{n}>(y0+62) & Northing{n}<(y0+63) & Easting{n}>(x0+16) & Easting{n}<(x0+20));
        I9=find(~(I5+I6+I7+I8));
        Easting{n}=Easting{n}(I9); % remove the bad points
        Northing{n}=Northing{n}(I9); % remove the bad points
        depth{n}=depth{n}(I9);
    end
    %pos=([Easting{n} Northing{n}]); % check order with Brian - this is flipped, but looks more likely
    %[h{n},Vexp1{n},npairs{n}]=semivar_exp(pos,depth{n},nbins); % calculate variogram of snow depth
    r(n) = variogram2D(Easting{n},Northing{n},depth{n},nbins,nAbin);
    MALL(n,1)=nanmean(depth{n});
    SALL(n,1)=nanstd(depth{n});
    %h2=log(h{n});
    %V2=log(Vexp1{n});
    %fh=@(p)model_variogram_error(h2,V2,npairs{n}(:),p(1),p(2),p(3),'L'); % make function handle for minimization
    %[pbest(n,:),fval(n)]=fminsearch(fh,log(G{n})) % simplex minimization, with initial guess
    fh=@(p)model_variogram_error(r(n).L,r(n).V,r(n).npairs,p(1),p(2),p(3),'S'); % make function handle for minimization
    [pbest(n,:),fval(n)]=fminsearch(fh,G{n}); % simplex minimization, with initial guess
    Vmod{n}=model_variogram(r(n).L,pbest(n,1),pbest(n,2),pbest(n,3),'S');
end


%% now plot
figure(1);clf
figure(2);clf
figure(3);clf
% first lets get the indexes for each site:
TL=[];
LDP=[];
UDC=[];
for n=1:length(MagnaFiles)
    switch MagnaFiles{n}(1)
        case 'T'
            TL=[TL;n];
         case 'L'
            LDP=[LDP;n];
        case 'U'
            UDC=[UDC;n];
    end
    
end
% lets make local origins for each site
x0TL=min(vertcat(Easting{TL}));
x0LDP=min(vertcat(Easting{LDP}));
x0UDC=min(vertcat(Easting{UDC}));
y0TL=min(vertcat(Northing{TL}));
y0LDP=min(vertcat(Northing{LDP}));
y0UDC=min(vertcat(Northing{UDC}));
xmTL=max(vertcat(Easting{TL}));
xmLDP=max(vertcat(Easting{LDP}));
xmUDC=max(vertcat(Easting{UDC}));
ymTL=max(vertcat(Northing{TL}));
ymLDP=max(vertcat(Northing{LDP}));
ymUDC=max(vertcat(Northing{UDC}));

%% now lets make plots for each site:
close all
figure(1); clf; set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
figure(6); clf; set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
marker2={'sb','ok','+r','xg','^m','vc','db','sy'};
marker={'ok','+r','xg','^m','vc','db','sy'};
for m=1:length(TL)
    n=TL(m);
    figure(1); 
    subplot(1,3,1)
    hp(n)=loglog(r(n).L,r(n).V,marker2{m},'markersize',8,'LineWidth',3);
    hold on
    plot(r(n).L,Vmod{n},marker2{m}(2),'LineWidth',5)
    figure(2); subplot(1,3,1)
    hp2(n)=plot(Easting{n}-x0TL,Northing{n}-y0TL,[marker2{m}],'markersize',4,'LineWidth',2);
    hold on
    figure(4); % temporary figure for pdf output
    [xval,density,N,xout] = plot_pdf_hist5(depth{n});
    figure(3); subplot(3,1,1)
    hp3(n)=plot(xval,density,marker2{m}(2),'LineWidth',4);
    hold on
    figure(5); subplot(1,3,1)
    hp4(n)=plot(r(n).L,r(n).Vgr,marker2{m},'markersize',8,'LineWidth',3);
    hold on
    figure(6); 
    subplot(1,3,1)
    hp5(n)=plot(r(n).L,100*sqrt(2*r(n).Vpr),marker2{m},'markersize',8,'LineWidth',3);
    hold on
end
figure(1); subplot(1,3,1)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp(1:5),Date{1:5},'Location','SouthEast')
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('\gamma [cm^2]','FontSize',Fts,'FontWeight','bold')
title('Treeline','FontSize',Fts,'FontWeight','bold')
axis([5 140 10 1000])
set(gca,'XTick',[10 20 50 100])
set(gca,'XTickLabel',{'10','20','50','100'})
grid on
figure(2); subplot(1,3,1)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp2(1:5),Date{1:5},'Location','Best')
xlabel('Easting [m]','FontSize',Fts,'FontWeight','bold')
ylabel('Northing [m]','FontSize',Fts,'FontWeight','bold')
title('Treeline','FontSize',Fts,'FontWeight','bold')
axis([0 xmTL-x0TL 0 ymTL-y0TL])
figure(3); subplot(3,1,1)
set(gca,'FontSize',Fs,'FontWeight','bold')
xlabel('depth [cm]','FontSize',Fts,'FontWeight','bold')
ylabel('pdf','FontSize',Fts,'FontWeight','bold')
title('Treeline','FontSize',Fts,'FontWeight','bold')
legend(hp3(TL),Date{TL},'Location','Best')
axis([0 200 0 100])
figure(5); subplot(1,3,1)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp4(1:5),Date{1:5},'Location','Best')
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('\gamma_{GR}','FontSize',Fts,'FontWeight','bold')
title('Treeline','FontSize',Fts,'FontWeight','bold')
figure(6); subplot(1,3,1)
set(gca,'FontSize',Fs,'FontWeight','bold')
%legend(hp5(1:5),Date{1:5},2)
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('% variability','FontSize',Fts,'FontWeight','bold')
title('Treeline','FontSize',Fts,'FontWeight','bold')
axis([0 150 0 160])
grid on
%
for m=1:length(LDP)
    n=LDP(m);
    figure(1); subplot(1,3,2)
    hp(n)=loglog(r(n).L,r(n).V,marker{m},'markersize',8,'LineWidth',3);
    hold on
    plot(r(n).L,Vmod{n},marker{m}(2),'LineWidth',5)
    figure(2); subplot(1,3,2)
    hp2(n)=plot(Easting{n}-x0LDP,Northing{n}-y0LDP,[marker{m}],'markersize',4,'LineWidth',2);
    hold on
    figure(4)
    [xval,density,N,xout] = plot_pdf_hist5(depth{n});
    figure(3); subplot(3,1,2)
    hp3(n)=plot(xval,density,marker{m}(2),'LineWidth',4);
    hold on
    figure(5); subplot(1,3,2)
    hp4(n)=plot(r(n).L,r(n).Vgr,marker{m},'markersize',8,'LineWidth',3);
    hold on
    figure(6); subplot(1,3,2)
    hp5(n)=plot(r(n).L,100*sqrt(2*r(n).Vpr),marker{m},'markersize',8,'LineWidth',3);
    hold on
end
%
figure(1); subplot(1,3,2)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp(LDP),Date{LDP},'Location','SouthEast')
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('\gamma [cm^2]','FontSize',Fts,'FontWeight','bold')
title('Lower Deer Point','FontSize',Fts,'FontWeight','bold')
axis([5 160 10 1000])
set(gca,'XTick',[10 20 50 100])
set(gca,'XTickLabel',{'10','20','50','100'})
grid on
figure(2); subplot(1,3,2)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp2(LDP),Date{LDP},'Location','Best')
xlabel('Easting [m]','FontSize',Fts,'FontWeight','bold')
ylabel('Northing [m]','FontSize',Fts,'FontWeight','bold')
title('Lower Deer Point','FontSize',Fts,'FontWeight','bold')
axis([0 xmLDP-x0LDP 0 ymLDP-y0LDP])
figure(3); subplot(3,1,2)
set(gca,'FontSize',Fs,'FontWeight','bold')
xlabel('depth [cm]','FontSize',Fts,'FontWeight','bold')
ylabel('pdf','FontSize',Fts,'FontWeight','bold')
title('Lower Deer Point','FontSize',Fts,'FontWeight','bold')
legend(hp3(LDP),Date{LDP},'Location','Best')
axis([0 200 0 60])
figure(5); subplot(1,3,2)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp4(LDP),Date{LDP},'Location','Best')
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('\gamma_{GR}','FontSize',Fts,'FontWeight','bold')
title('Lower Deer Point','FontSize',Fts,'FontWeight','bold')
figure(6); subplot(1,3,2)
set(gca,'FontSize',Fs,'FontWeight','bold')
%legend(hp5(LDP),Date{LDP},2)
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('% variability','FontSize',Fts,'FontWeight','bold')
title('Lower Deer Point','FontSize',Fts,'FontWeight','bold')
axis([0 150 0 160])
grid on
%
for m=1:length(UDC)
    n=UDC(m);
    figure(1); subplot(1,3,3)
    hp(n)=loglog(r(n).L,r(n).V,marker{m},'markersize',8,'LineWidth',3);
    hold on
    plot(r(n).L,Vmod{n},marker{m}(2),'LineWidth',5)
    figure(2); subplot(1,3,3)
    hp2(n)=plot(Easting{n}-x0UDC,Northing{n}-y0UDC,[marker{m}],'markersize',4,'LineWidth',2);
    hold on
    figure(4)
    [xval,density,N,xout] = plot_pdf_hist5(depth{n});
    figure(3); subplot(3,1,3)
    hp3(n)=plot(xval,density,marker{m}(2),'LineWidth',4);
    hold on
    figure(5); subplot(1,3,3)
    hp4(n)=plot(r(n).L,r(n).Vgr,marker{m},'markersize',8,'LineWidth',3);
    hold on
    figure(6); subplot(1,3,3)
    hp5(n)=plot(r(n).L,100*sqrt(2*r(n).Vpr),marker{m},'markersize',8,'LineWidth',3);
    hold on
end
figure(1); subplot(1,3,3)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp(UDC),Date{UDC},'Location','SouthEast')
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('\gamma [cm^2]','FontSize',Fts,'FontWeight','bold')
title('Upper Dry Creek','FontSize',Fts,'FontWeight','bold')
axis([5 240 10 1000])
set(gca,'XTick',[10 20 50 100])
set(gca,'XTickLabel',{'10','20','50','100'})
grid on
figure(2); subplot(1,3,3)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp2(UDC),Date{UDC},'Location','Best')
xlabel('Easting [m]','FontSize',Fts,'FontWeight','bold')
ylabel('Northing [m]','FontSize',Fts,'FontWeight','bold')
title('Upper Dry Creek','FontSize',Fts,'FontWeight','bold')
axis([0 xmUDC-x0UDC 0 ymUDC-y0UDC])
figure(3); subplot(3,1,3)
set(gca,'FontSize',Fs,'FontWeight','bold')
xlabel('depth [cm]','FontSize',Fts,'FontWeight','bold')
ylabel('pdf','FontSize',Fts,'FontWeight','bold')
title('Upper Dry Creek','FontSize',Fts,'FontWeight','bold')
legend(hp3(UDC),Date{UDC},'Location','Best')
axis([0 200 0 60])
figure(5); subplot(1,3,3)
set(gca,'FontSize',Fs,'FontWeight','bold')
legend(hp4(UDC),Date{UDC},'Location','Best')
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('\gamma_{GR}','FontSize',Fts,'FontWeight','bold')
title('Upper Dry Creek','FontSize',Fts,'FontWeight','bold')
figure(6); subplot(1,3,3)
set(gca,'FontSize',Fs,'FontWeight','bold')
%legend(hp5(UDC),Date{UDC},2)
xlabel('lag [m]','FontSize',Fts,'FontWeight','bold')
ylabel('% variability','FontSize',Fts,'FontWeight','bold')
title('Upper Dry Creek','FontSize',Fts,'FontWeight','bold')
axis([0 150 0 160])
grid on
% lets plot the evolution of depth distributions at the three sites
D=[];
G=[];
for n=1:length(TL)
    D=[D;depth{TL(n)}];
    G=[G;n*ones(length(depth{TL(n)}),1)];
end
figure(4);clf
subplot(1,3,1)
hb=boxplot(D,G);
set(hb,'LineWidth',3)
set(gca,'FontSize',Fs,'FontWeight','bold')
title('Treeline')
ylabel('depth [cm]')
set(findobj(gca,'Type','text'),'FontSize',12,'FontWeight','bold')
axis([0.5 5.5 0 200])

D=[];
G=[];
for n=1:length(LDP)
    D=[D;depth{LDP(n)}];
    G=[G;n*ones(length(depth{LDP(n)}),1)];
end
subplot(1,3,2)
hb=boxplot(D,G);
set(hb,'LineWidth',3)
set(gca,'FontSize',Fs,'FontWeight','bold')
title('Lower Deer Point')
ylabel('depth [cm]')
set(findobj(gca,'Type','text'),'FontSize',12,'FontWeight','bold')
axis([0.5 7.5 0 200])

D=[];
G=[];
for n=1:length(UDC)
    D=[D;depth{UDC(n)}];
    G=[G;n*ones(length(depth{UDC(n)}),1)];
end
subplot(1,3,3)
hb=boxplot(D,G);
set(hb,'LineWidth',3)
set(gca,'FontSize',Fs,'FontWeight','bold')
title('Upper Dry Creek','FontSize',Fts,'FontWeight','bold')
ylabel('depth [cm]','FontSize',Fts,'FontWeight','bold')
set(findobj(gca,'Type','text'),'FontSize',12,'FontWeight','bold')
axis([0.5 5.5 0 200])
%%
figure(1);
set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
print('-dpng','-r300','FIGURES/DryCreekVariograms.png')
%%
figure(2);
set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
print('-dpng','FIGURES/SurveyLocations.png')
figure(3);
set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
print('-dpng','FIGURES/depth_pdfs.png')
figure(4);
set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
print('-dpng','FIGURES/depth_boxplots.png')
figure(5);
set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
print('-dpng','FIGURES/GeneralRelativeVariograms.png')
%%
figure(6);
set(gcf,'PaperPositionMode','auto','position',[1 35 1280 671])
print('-dpng','-r300','FIGURES/PercentVariability.png')

%% now lets make a table










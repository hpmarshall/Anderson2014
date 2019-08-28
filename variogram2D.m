function r = variogram2D(x,y,m,lagbins,nAbin)
% HPM 08/19/11
% calculates the semivariogram on 2D data
% INPUT: x = x-coordinate
%        y = y-coordinate
%        m = measurement
%     lagbins = array of lag bin edges for the variogram [0,5,10], or number of equal spaced bins
%   nAbin = number of angle bins (1,2, or 4)
% OUTPUT: r = results structure array
%         r.L = vector of bin centers
%         r.A = vector of direction angle bin centers [90,270 (E/W) 0,180(N/S)]
%         r.mHT = mean measurement value [heads tails] in bin
%         r.npairs = number of point pairs in each lag/distance bin
%         r.G = [nlags,50] = 50 random samples of the differences (m_i-m_j)^2 in each lag bin (for speed, keep it same size)
%         r.V = semivariogram = 1/(2*N) sum[(m_i-m_j)^2]
%         r.Vgr = general relative semivariogram = r.V./mean(r.mHT)
%         r.Vpr = pairwise relative semivariogram (Goovaerts, 1997, p.85)

ndata=length(x); % number of data points
npt = sum(1:(length(x)-1)); % total number of pairs of points
% initialize vectors
head=zeros(npt,3); % location and starting value for each pair
tail=zeros(npt,3); % location and ending value for each pair
dr=zeros(npt,2); % delta distance vector 
lag=zeros(npt,1); % lags
gamma=zeros(npt,1); % variance of each pair (Vi-Vj)^2
ang=zeros(npt,1); % angle of each difference vector
i=0; % initialize index
for i1=1:(ndata-1) % loop over all heads
    for i2=(i1+1):ndata % loop over all tails
        i=i+1; % update index
        head(i,:)=[x(i1) y(i1) m(i1)]; % position, value of head
        tail(i,:)=[x(i2) y(i2) m(i2)]; % poistion, value of tail
        dr(i,:)=head(i,1:2)-tail(i,1:2); % difference vector
        lag(i)=sqrt(sum(dr(i,1).^2+dr(i,2).^2)); % distance between head and tail
        gamma(i)=(head(i,3)-tail(i,3)).^2; % squared difference in value between points
        ang(i)=180/pi*atan(dr(i,2)./dr(i,1)); % direction, -90=S,0=E,90=N,180=W
    end
end
% now bin the data and calculate statistics
if length(lagbins)==1 % if number of lagbins specified
    lagbins=linspace(0,(max(lag)/2),lagbins);
end
switch nAbin
    case 1 % if omnidirectional variogram
        mHT=zeros(length(lagbins)-1,2); % mean value for heads and tails
        npairs=zeros(length(lagbins)-1,1); % number of point pairs
        G=zeros(length(lagbins)-1,50); % random sample of difference
        V=zeros(length(lagbins)-1,1); % semivariance
        Vgr=zeros(length(lagbins)-1,1); % general relative semivariance
        Vpr=zeros(length(lagbins)-1,1); % pairwise relative semivariance
        A=NaN; % no angle bins
        for i3=1:(length(lagbins)-1) % loop over lag bins
            ind=find(lag>=lagbins(i3) & lag<lagbins(i3+1)); % index to lagbins within current edges
            mHT(i3,:)=[mean(head(ind,3)) mean(tail(ind,3))]; % mean of heads and tails
            npairs(i3)=length(ind); % number of pairs
            G(i3,:)=randsample(gamma(ind),50,true); % sample of squared differences
            V(i3)=0.5*mean(gamma(ind)); % semivariance
            Vgr(i3)=V(i3)/mean(mHT(i3,:)); % general relative semivariance
            HT=(head(ind,3)+tail(ind,3))./2; % mean value for each pair
            I2=find(HT>0); % only use non-zero HT values
            Vpr(i3)=0.5*mean(gamma(ind(I2))./(HT(I2).^2)); % pairwise relative semivar (Goovaerts, 97, p.85)
        end
    case 2 % if bi-directional variogram
        mHT=zeros(length(lagbins)-1,2,2); % mean value for heads and tails
        npairs=zeros(length(lagbins)-1,2); % number of point pairs
        G=zeros(length(lagbins)-1,2,50); % random sample of difference
        V=zeros(length(lagbins)-1,2); % semivariance
        Vgr=zeros(length(lagbins)-1,2); % general relative semivariance
        Vpr=zeros(length(lagbins)-1,2); % pairwise relative semivariance
        A={'N/S','E/W'}; % angle bins
        for i3=1:(length(lagbins)-1) % loop over lag bins
            ind1=find((lag>=lagbins(i3) & lag<lagbins(i3+1)) & (ang>=-45 & ang<45)); % index to E/W and current lagbin
            ind2=find((lag>=lagbins(i3) & lag<lagbins(i3+1)) & (ang>=45 | ang<-45)); % index to N/S and current lagbin
            mHT(i3,1,:)=[mean(head(ind2,3)) mean(tail(ind2,3))]; % mean of heads and tails in N/S
            mHT(i3,2,:)=[mean(head(ind1,3)) mean(tail(ind1,3))]; % mean of heads and tails in E/W
            npairs(i3,:)=[length(ind2) length(ind1)]; % number of pairs in N/S and E/W
            G(i3,1,:)=randsample(gamma(ind2),50,true); % sample of squared differences in N/S
            G(i3,2,:)=randsample(gamma(ind1),50,true); % sample of squared differences in E/W
            V(i3,:)=0.5*[mean(gamma(ind2)) mean(gamma(ind1))]; % semivariance
            Vgr(i3,:)=V(i3,:)./[mean(mHT(i3,1,:)) mean(mHT(i3,2,:))]; % general relative semivariance
            HT1=(head(ind2,3)+tail(ind2,3))./2; % mean value for each pair
            HT2=(head(ind1,3)+tail(ind1,3))./2; % mean value for each pair
            I2=find(HT1>0); I4=find(HT2>0);% only use non-zero HT values
            Vpr(i3,:)=0.5*[mean(gamma(ind2(I2))./(HT1(I2).^2)) mean(gamma(ind1(I4))./(HT2(I4).^2))]; % pairwise relative semivariance
        end
    case 4 % if quad-directional variogram
        mHT=zeros(length(lagbins)-1,4,2); % mean value for heads and tails
        npairs=zeros(length(lagbins)-1,4); % number of point pairs
        G=zeros(length(lagbins)-1,4,50); % random sample of difference
        V=zeros(length(lagbins)-1,4); % semivariance
        Vgr=zeros(length(lagbins)-1,4); % general semivariance
        Vpr=zeros(length(lagbins)-1,4); % pairwise semivariance
        A={'N/S','E/W','NE/SW','SE/NW'}; % angle bins
        for i3=1:(length(lagbins)-1) % loop over lag bins
            ind1=find((lag>=lagbins(i3) & lag<lagbins(i3+1)) & (ang>=-22.5 & ang<22.5)); % index to E/W and current lagbin
            ind2=find((lag>=lagbins(i3) & lag<lagbins(i3+1)) & (ang>=67.5 | ang<-67.5)); % index to N/S and current lagbin
            ind3=find((lag>=lagbins(i3) & lag<lagbins(i3+1)) & (ang>=22.5 & ang<67.5)); % index to NE/SW and current lagbin
            ind4=find((lag>=lagbins(i3) & lag<lagbins(i3+1)) & (ang>=-67.5 | ang<-22.5)); % index to SE/NW and current lagbin
            mHT(i3,1,:)=[mean(head(ind2,3)) mean(tail(ind2,3))]; % mean of heads and tails in N/S
            mHT(i3,2,:)=[mean(head(ind1,3)) mean(tail(ind1,3))]; % mean of heads and tails in E/W
            mHT(i3,3,:)=[mean(head(ind3,3)) mean(tail(ind3,3))]; % mean of heads and tails in N/S
            mHT(i3,4,:)=[mean(head(ind4,3)) mean(tail(ind4,3))]; % mean of heads and tails in E/W            
            npairs(i3,:)=[length(ind2) length(ind1) length(ind3) length(ind4)]; % number of pairs in all four directions
            G(i3,1,:)=randsample(gamma(ind2),50,true); % sample of squared differences in N/S
            G(i3,2,:)=randsample(gamma(ind1),50,true); % sample of squared differences in E/W
            G(i3,3,:)=randsample(gamma(ind3),50,true); % sample of squared differences in NE/SW
            G(i3,4,:)=randsample(gamma(ind4),50,true); % sample of squared differences in SE/NW
            V(i3,:)=0.5*[mean(gamma(ind2)) mean(gamma(ind1)) mean(gamma(ind3)) mean(gamma(ind4))]; % semivariance
            Vgr(i3,:)=V(i3,:)./[mean(mHT(i3,1,:)) mean(mHT(i3,2,:)) mean(mHT(i3,3,:)) mean(mHT(i3,4,:))]; % generalized relative semivariance
            HT1=(head(ind2,3)+tail(ind2,3))./2; % mean value for each pair
            HT2=(head(ind1,3)+tail(ind1,3))./2; % mean value for each pair
            HT3=(head(ind3,3)+tail(ind3,3))./2; % mean value for each pair
            HT4=(head(ind4,3)+tail(ind4,3))./2; % mean value for each pair
            I2=find(HT1>0); I4=find(HT2>0);% only use non-zero HT values
            I5=find(HT3>0); I6=find(HT4>0);% only use non-zero HT values
            Vpr(i3,:)=0.5*[mean(gamma(ind2(I2))./(HT1(I2).^2)) mean(gamma(ind1(I4))./(HT2(I4).^2)) mean(gamma(ind3(I5))./(HT3(I5).^2)) mean(gamma(ind4(I6))./(HT4(I6).^2))]; % pairwise relative semivariance

        end
    otherwise
        disp('only omnidirectional, bidirectional, and quad-directional currently implemented')
        A=NaN; mHT=NaN; npairs=NaN; G=NaN;
end
% store results in cell array
lag=lagbins(1:(end-1))+(diff(lagbins)/2); % get bin centers from bin edges
r.L=lag; r.A=A; r.mHT=mHT; r.npairs=npairs; r.G=G; r.V=V; r.Vgr=Vgr; r.Vpr=Vpr;
%r=struct('L',lag,'A',A,'mHT',mHT,'npairs',npairs,'G',G,'V',V,'Vgr',Vgr,'Vpr',Vpr);
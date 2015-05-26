function [ZP,PP]=FlipAngleMinVar(Nq,lb,ub,IX,IY)
% loading measurement data, Note rhodata is calculated by dividing the
% corresponding data from two different flip angle, for detail see the
% following:
% MATLAB/projects/mfgreScripts/ooMFGRE
load rhodata    % data directory: /mnt/FUS4/data2/madankan/FlipAngle

% number of uncertain parameters (Nd=2 here)
Nd=length(lb);
% assigning the quadrature points and their weights for the uncertain parameters
for ind=1:Nd
    [xtmp,wtmp]=lgwt(Nq,lb(ind),ub(ind));
    x(:,ind)=xtmp;
    wtmp=wtmp/(ub(ind)-lb(ind));
end
index=GenerateIndex(Nd,Nq*ones(1,Nd));
para=zeros(size(index));
weight=ones(Nq^Nd,1);
for j=1:Nd
    para(:,j)=x(index(:,j),j);
    weight=weight.*wtmp(index(:,j));    
end

% calculating prior statistics of uncertain parameters
for j=1:Nd
    zm(j,1)=sum(para(:,j).*weight); % prior mean
    Pm(j,j)=sum((para(:,j)-zm(j)).^2.*weight);% prior covariance
end

Temp_ref=0;
T1ref=para(:,1);
slope=para(:,2);

% defining the parameters
theta=pi/6;
beta=pi/12;
Tr=45;

% performing minimum variance estimation for each spatial location over
% time:
for ix=IX;% assign the x coordinate
    for iy=IY;% assign the y coordinate
        Ez=[];
        Rho=[];
        ZP=zeros(60,2);
        PP=cell(60,1);
        for tind=1:60 % loop thru the time steps
            T1=T1ref+slope.*T1ref.*(Temp(ix,iy,tind)-Temp_ref); % quadrature values of T1
            rho=sin(theta).*(1-exp(-Tr./T1))./sin(beta)./(1-cos(theta)*exp(-Tr./T1)); % model predictions of rho
            Exp_rho=sum(rho.*weight); % expected value of model predicted rho
            
            Ez=[Ez; Exp_rho];
            Rho=[Rho;rho'];
            
            % required statistics for minimum variance estimation:
            Prr=zeros(tind);
            Prp=zeros(Nd,tind);
            for i=1:tind
               for j=1:tind
                   Prr(i,j)=sum((Rho(i,:)-Ez(i)).*(Rho(j,:)-Ez(j)).*weight');
               end
               for pind=1:Nd
                   Prp(pind,i)=sum((Rho(i,:)-Ez(i))'.*(para(:,pind)-zm(pind)).*weight);
               end
            end

            % measurement noise
            R=diag(0.01*squeeze(rho_meas(IX,IY,1:tind))); 
            K= Prp/(Prr+R); % Kalman Gain!
            zp=zm+K*(squeeze(rho_meas(ix,iy,1:tind))-Ez); % posterior mean
            Pp=Pm-K*Prr*K';% posterior covariance

            ZP(tind,:)=zp;
            PP{tind}=Pp;
        end
    end
end
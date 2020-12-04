function [ c, dC_drho, dC_dF, dC_dXg ] = constrain_fcn(p, x, volfrac, nelx, nely, U, edofMat,...
    KE, freedofs, S, Sl, iK, jK, K, F, Xg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% xD = zeros(nely*nelx,1);
% x = reshape(x,nely,nelx);
% x(passive==1) = 0;
% x(passive==2) = 0;
% x(passive==3) = 1;
% xD(eFree) = x(1:xsize);
% xF = x(xsize+1:end);
% x = zeros(nely,nelx);
% xTilde = x;
% if ft == 1
%     x(:) = xD;
% elseif ft == 2
%     x(:) = (H*xD(:))./Hs;
% elseif ft == 3
%     xTilde(:) = (H*xD(:))./Hs;
%     x(:) = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
% end
% x(end,25:end) = 1;
% x(passive==1) = 0;
% x(passive==2) = 0;
% x(passive==3) = 1;

%% Volume Fraction
c(1) = 1*(mean(x(:))-volfrac);
dC1=ones(size(x))/length(x(:));
%% Stress Constraint
P=4;
mse = reshape(sqrt(abs(sum((U(edofMat)*S).*U(edofMat),2))),nely,nelx); %microscopic Von Mises Stress
rcv=x(:).*(mse(:)/Sl-1);
sS=reshape(S(:)*(x(:)'./mse(:)'/Sl.*exp(P*(rcv(:)'-max(rcv(:)))))/sum(exp(P*(rcv(:)-max(rcv(:))))),64*nelx*nely,1);
S0=sparse(iK,jK,sS); S0=(S0+S0')/2;
dG_ksdu=S0*U;
Lambda=U;
Lambda(freedofs) = K(freedofs,freedofs)\dG_ksdu(freedofs);
dG_ksdx=(mse(:)/Sl-1).*exp(P*(rcv(:)-max(rcv(:))))/sum(exp(P*(rcv(:)-max(rcv(:)))));
% se=mse.*(Emin+reshape(x,nely,nelx).^penal*(p.E0-Emin))/E0;
c(2)=max(rcv(:))+1/P*log(sum(exp(P*(rcv(:)-max(rcv(:))))))-log(length(rcv(:)))/P;
dc_du = reshape((sum((U(edofMat)*KE).*Lambda(edofMat),2)),nely,nelx);
dc_du = -p.penalty*(p.E0-p.Emin)*x(:).^(p.penalty-1).*dc_du(:);
dC2_dx = dc_du(:)+dG_ksdx;
% dC2_dF = Lambda;
dC2_df = Lambda;
%% Compliance Constraint
c0 = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
c(3) = sum(sum((p.Emin+x.^p.penalty*(p.E0-p.Emin)).*c0))-20;
lambda2 = 2*U;
dC3_dx = -p.penalty*(p.E0-p.Emin)*x(:).^(p.penalty-1).*c0(:);
dC3_df = lambda2';
%% Force component position and length constraint
c(4) = Xg(2)^2 - nely*Xg(2) + nely/2*Xg(3) - Xg(3)^2/4;
dC_dXg = sparse([2;3],4,[2*Xg(2)-nely; nely - Xg(3)/2],length(Xg),4);

% Fmin = 1.5*Fmin;
% 
% c(4) = 1-sum(F)/Fmin;
% dC4_dx = -sum(df_drho)'/Fmin;
% dC4_df = -sum(df_fi(:,forcedofs))/Fmin;

% c(4) = Fmin/sum(F)-1;
% dC4_dx = -Fmin./sum(F).^2*sum(df_drho)';
% dC4_df = -Fmin./sum(F).^2*sum(df_fi(:,forcedofs));

% c(4) = Fmin/sum(xF)-1;
% dC4_dx = zeros(size(xD));
% dC4_df = -Fmin/sum(xF)^2*ones(size(xF));

% c(4) = Fmin/sum(wP)-1;
% dC4_dx = zeros(size(xD));
% dC4_df = -Fmin/sum(wP)^2*diag(df_fi(forcedofs,forcedofs));

% if ft == 1
% elseif ft == 2
%     dC1 = H*(dC1(:)./Hs);
%     dC2_dx = H*(dC2_dx./Hs);
%       dC3_dx = H*(dC3_dx./Hs);
%       dC4_dx = H*(dC4_dx./Hs);
% elseif ft == 3
%     dx = beta*exp(-beta*xTilde)+exp(-beta);
%     dC1 = H*(dC1(:).*dx(:)./Hs);
%     dC2_dx = H*(dC2_dx.*dx(:)./Hs);
%      dC3_dx = H*(dC3_dx.*dx(:)./Hs);
% end
%  dC = [dC1 dC2_dx ; zeros(length(x)-length(xD),1) dC2_df'];
%dC = [dC1 dC2_dx dC3_dx; zeros(length(x)-length(xD),1) dC2_df' dC3_df'];
dC_drho = [dC1(:) dC2_dx dC3_dx sparse(numel(x),1)];

dC_dF = [sparse(length(F),1) dC2_df dC3_df' sparse(length(F),1)];

% c(3) = 1*c(3);
% dC(:,3) = 1*dC(:,3);

c(1) = 1*c(1);
dC_drho(:,1) = 1*dC_drho(:,1);
dC_dF(:,1) = 1*dC_dF(:,1);

c(2) = 1*c(2);
dC_drho(:,2) = 1*dC_drho(:,2);
dC_dF(:,2) = 1*dC_dF(:,2);

c(3) = .01*c(3);
dC_drho(:,3) = .01*dC_drho(:,3);
dC_dF(:,3) = 1*dC_dF(:,3);

c(4) = 1*c(4);
dC_drho(:,4) = 1*dC_drho(:,4);
dC_dF(:,4) = 1*dC_dF(:,4);
dC_dXg(:,4) = 1*dC_dXg(:,4);
end



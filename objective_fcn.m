function [o, do, U, do_df, K] = objective_fcn(p, x, KE, iK, jK, u_obj, nelx, nely, freedofs, edofMat,L, F)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global loop
loop=loop+1;
%% FE-ANALYSIS
% xD = zeros(nely*nelx,1);
% xD(eFree) = x(1:xsize);
% x = reshape(x,nely,nelx);
% x(passive==1) = 0;
% x(passive==2) = 0;
% x(passive==3) = 1;
% x(end,25:end) = 1;

% xPhys = zeros(nely,nelx);
% xTilde = xPhys;
% if ft == 1
%     xPhys(:) = xD;
% elseif ft == 2
%     xPhys(:) = (H*xD(:))./Hs;
% elseif ft == 3
%     xTilde(:) = (H*xD(:))./Hs;
%     xPhys(:) = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
% end
% % xPhys(end,25:end) = 1;
% xPhys(passive==1) = 0;
% xPhys(passive==2) = 0;
% xPhys(passive==3) = 1;

U = zeros(2*(nely+1)*(nelx+1),1);
sK = reshape(KE(:)*(p.Emin+x(:)'.^p.penalty*(p.E0-p.Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
% F = sparse(forcedofs,ones(size(wP)), wP.*xD(1:length(forcedofs))'.^(1/p.penalty), length(U),1);
U(freedofs) = K(freedofs,freedofs)\F(freedofs);

% xPhys_support = [zeros(1,nelx+2); zeros(nely,1) xPhys zeros(nely,1); zeros(1,nelx+2)];
% 
% xPhys_ul = xPhys_support(1:end-1,1:end-1);
% xPhys_ur = xPhys_support(1:end-1,2:end);
% xPhys_dl = xPhys_support(2:end,1:end-1);
% xPhys_dr = xPhys_support(2:end,2:end);
% 
% xPhys_ul(1,:) = 2*xPhys_ul(1,:);
% xPhys_ur(1,:) = 2*xPhys_ur(1,:);
% xPhys_dl(1,:) = 2*xPhys_dl(1,:);
% xPhys_dr(1,:) = 2*xPhys_dr(1,:);
% 
% xPhys_ul(end,:) = 2*xPhys_ul(end,:);
% xPhys_ur(end,:) = 2*xPhys_ur(end,:);
% xPhys_dl(end,:) = 2*xPhys_dl(end,:);
% xPhys_dr(end,:) = 2*xPhys_dr(end,:);
% 
% xPhys_ul(:,1) = 2*xPhys_ul(:,1);
% xPhys_ur(:,1) = 2*xPhys_ur(:,1);
% xPhys_dl(:,1) = 2*xPhys_dl(:,1);
% xPhys_dr(:,1) = 2*xPhys_dr(:,1);
% 
% xPhys_ul(:,end) = 2*xPhys_ul(:,end);
% xPhys_ur(:,end) = 2*xPhys_ur(:,end);
% xPhys_dl(:,end) = 2*xPhys_dl(:,end);
% xPhys_dr(:,end) = 2*xPhys_dr(:,end);
% 
% rho_nodes = (xPhys_ul+xPhys_ur+xPhys_dl+xPhys_dr)/4;
% rho2 = [rho_nodes(:) rho_nodes(:)]';
% rho_n = reshape(rho2,length(U),1);
% Ux = U.*rho_n;
% Kx = K*spdiags(1./rho_n,0,2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));

o = norm(L*U-u_obj)^2;
do_dUf = 2*U'*(L')*L-2*u_obj'*L;
% o = 100 - U(2*(nelx+1)*(nely+1)-nely);
% do_dUf = sparse(1,2*(nelx+1)*(nely+1)-nely,-1,1,2*(nelx+1)*(nely+1));
% do_dUf = do_dUf(freedofs);
lambda = U';
lambda(freedofs) = (-K(freedofs,freedofs)\do_dUf(freedofs)')';
do_df = -lambda;
% o = norm(L*Ux-u_obj)^2;
% do_dUf = (2*Ux'*L'*L-2*u_obj'*L).*rho_n;
% do_dUf = do_dUf(freedofs)';
% lambda = U';
% lambda(freedofs) = -do_dUf/K(freedofs,freedofs);

% xx = xD(:);
% df_fi = zeros(size(U));
% df_fi(forcedofs) = Fmma*s/atan(s)./(s^2*xF.^2+1).*xD(1:length(forcedofs))'.^(1/p.penalty);
%df_fi(forcedofs) = 1;
% df_fi = spdiags(df_fi,0,length(U),length(U));
% df_drho = sparse(forcedofs,1:length(forcedofs),1/p.penalty*Fmma*s/atan(s)./(s^2*xF.^2+1).*x(1:length(forcedofs))'.^(1/p.penalty-1),length(U),nelx*nely);
% sK_rho1 = reshape((U(edofMat)*KE).*repmat(xPhys(:).^(penal-1)*penal*(E0-Emin),1,8),8*nelx*nely,1);
% iKK = reshape(edofMat,8*nelx*nely,1);
% jKK = reshape(repmat([1:nelx*nely]',1,8),8*nelx*nely,1);
% iKK = edofMat;
% sK_rho1 = (U(edofMat)*KE).*repmat(xPhys(:).^(penal-1)*penal*(E0-Emin),1,8);
% jKK = repmat([1:nelx*nely]',1,8);
% dKU_drho = sparse(iKK,jKK,sK_rho1);
c00 = reshape(sum((lambda(edofMat)*KE).*U(edofMat),2),nely,nelx) ; %initial lambda compliance
do =p.penalty*(p.E0-p.Emin)*x.^(p.penalty-1).*c00;
do = reshape(do,1,nelx*nely);
%do =lambda*dKU_drho(freedofs,:);
% if ft == 1
%     do(:) = H*(x(:).*do(:))./Hs./max(1e-3,x(:));
% elseif ft == 2
%     do(:) = H*(do(:)./Hs);
% elseif ft == 3
%     dx = beta*exp(-beta*xTilde)+exp(-beta);
%     do(:) = H*(do(:).*dx(:)./Hs);
% end
% do = [do, -lambda(freedofs)*df_fi(freedofs,forcedofs)];

% figure(5)
% hold on
% plot(loop,o,'bo','MarkerFaceColor','b')
% plot(loop,c(1)*100,'ro','MarkerFaceColor','r')
% plot(loop,c(2)*100,'ro','MarkerFaceColor','r')
% % plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
% title(['Convergence volfrac = ',num2str(mean(xPhys(:))*100),', G_{KS}^l =',num2str(c*100),'%, iter = ', num2str(loop)])
% grid on
% legend('Obj %','C1 %','C2 %')
% xlabel('iter')
end



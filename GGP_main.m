%% Generalized Geometry Projection
%-------------------------------------------------------------
%    This is the file GGP_main.m you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation; either version 3 of 
%    the License, or (at your option) any later version.
%    
%    This code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    
%    You should have received a copy of the GNU General Public License
%    (file COPYING) along with this file.  If not, see 
%    <http://www.gnu.org/licenses/>.
%
%    Version Nov 2019.
%    Simone Coniglio <simone.coniglio@airbus.com>
%    Propulsion Airframe Stress Transverse,
%    31300 Toulouse, France.
%
% This is an introduction to a Matlab implementation of Generalized Geometry 
% Projection approach for topology optimization.
% 
% In this approach geometric primitives are projected on a Finite Element 
% Mesh and assembled together to build the solution. 
% Author Simone Coniglio,12/09/2019
%% Problem set-up
function GGP_main(nelx,nely,volfrac,BC,method)
% In this section of the Matlab code we define several *parameters* needed for 
% the *Generalized Geometry Projection*.
% GGP parameters
tic
switch BC 
    case 'RIB'
        afoil = importdata('e420.dat');
        ax = afoil.data(2:end,1);
        ay = afoil.data(2:end,2);
        xmaxi = find(ax == max(ax));
        ax = [ax(1:xmaxi(1));flipud(ax(xmaxi(1)+1:xmaxi(2)))];
        ay = [ay(1:xmaxi(1));flipud(ay(xmaxi(1)+1:xmaxi(2)))];
        ay = ay+abs(min(ay));
        ind = find(ax>=0.15 & ax<=0.7);
        rx = ceil(500*ax(ind));
        ry = ceil(500*ay(ind));
        if find(min(rx)) < length(rx)/2
            rx(end) = min(rx);
        else
            rx(1) = min(rx);
        end
        rxmaxi = find(rx==max(rx));
        if (max(rx) - rx(rxmaxi+1)) < (max(rx) - rx(rxmaxi-1))
            rx(rxmaxi+1) = max(rx);
        else
            rx(rxmaxi-1) = max(rx);
        end
        rx = rx - min(rx);
        ry = ry - min(ry);
        nelx = ceil(max(rx));
        nely = ceil(max(ry));
end
stopping_criteria='kktnorm'; %stopping criteria of the optimization algorithm either change or KKT norm
% nelx=100;nely=100; 
% BC='Compliant';%L-shape %Short_Cantilever%MBB
p.method=method;%MMC%MNA %GP this change the function employed for the evaluation of local volume fraction
switch method
    case 'MMC'
        q=3;
    otherwise
        q=1;
end
% q=3;
p.zp=1 ;% parameter for p-norm/mean regularization
p.alp=1; %parameter for MMC
p.epsi=0.866;% parameter for MMC
p.bet=1e-3; %parameter for MMC
p.deltamin=1e-6; %parameter for GP
p.r=1.5;%parameter for GP
minh=3;% minimal bar thickness
p.sigma=3;%parameter for MNA
p.gammav=1;%parameter for GP
p.gammac=3;%parameter for GP
p.penalty=3;%parameter for MNA
p.aggregation='KSl'; %parameter for the aggregation function to be used
% IE= Induced Exponential % KS= KS function %KSl= lowerbound KS function
% p-norm %p-mean
p.ka=10; % parameter for the aggregation constant
p.saturation=true; % switch for saturation
ncx=2; % number of components in the x direction
ncy=2; % number of components in the y direction
Ngp=2; % number of Gauss point per sampling window
R=0.5; % radius of the sampling window (infty norm)
initial_d=0.5; % initial mass variable adopted for MNA and GP
mult = 1; % Volfrac modifier in case of emptyelts != 0
%% Generate a *folder* and a prefix to save images *optimization history*:

rs=replace(num2str(R,'%3.2f'),'.','_');
vf = replace(num2str(volfrac,'%3.2f'),'.','_');
folder_name=['Optimization_history_',BC,'_',p.method,'_Volfrac_',vf,'_nelx_',num2str(nelx),...
    '_nely_',num2str(nely),'_ncx&y_',num2str(ncx),'_',num2str(ncy),'_R_',rs,'_Ngp_',num2str(Ngp),'_SC_',stopping_criteria];
image_prefix=[BC,'_',p.method,'_Volfrac_',vf,'_nelx_',num2str(nelx),'_nely_',num2str(nely),'_R_',rs,'_Ngp_',num2str(Ngp)];
mkdir(folder_name)
Path=[folder_name,'/'];
%% MATERIAL PROPERTIES
p.E0 = 1;
p.Emin = 1e-6;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nele = nelx*nely;
nnodes = (nely+1)*(nelx+1);
ndof = 2*nnodes;
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nele,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nele,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nele,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nele,1);
% U = zeros(ndof,1);
%define the nodal coordinates
[Yy,Xx]=find(nodenrs);
Yy=nely+1-Yy;
Xx=Xx-1;
% Element connectivity
enodeMat=edofMat(:,[2,4,6,8])/2;
%% Prepare the *Generalized Geometry Projection:*
% Compute the element centroid coordinates:

xc=mean(Xx(enodeMat'));
yc=mean(Yy(enodeMat'));
centroid_coordinate=[xc(:),yc(:)];
%% DEFINE LOADS AND SUPPORTS
U = sparse(ndof,1);
fixeddofs = [];
emptyelts = [];
fullelts = [];
switch BC
    case 'MBB'
        excitation_node=1;excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=[find(Xx==min(Xx));(nelx+1)*(nely+1)];fixed_dir=[ones(nely+1,1);2];
        fixeddofs=2*(fixednodes-1)+fixed_dir;
    case 'Short_Cantilever'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=repmat(find(Xx==min(Xx)),2,1);fixed_dir=[ones(nely+1,1);2*ones(nely+1,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
    case 'L_Shape'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=repmat(find(Yy==max(Yy)),2,1);fixed_dir=[ones(nelx+1,1),2*ones(nelx+1,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
        emptyelts=find(xc>=(((max(Xx)+min(Xx))/2))&(yc>=((max(Yy)+min(Yy))/2)));
    case 'Compliant'
        in = 1; out = 2*(nely+1)*nelx+1;
        F = sparse([in;out],[1;2],[1; -1],ndof,2);
        U = zeros(ndof,2);
        fixednodes = [1:(nely+1):(nelx+1)*(nely+1) nely+1 nely+1]';
        fixed_dir = [2*ones(nelx+1,1);1;2];
        fixeddofs = 2*(fixednodes-1)+fixed_dir(:);
%         fixedx = [2:2*(ne1ly+1):ndof]';
%         fixeddofs = union(fixedx,[2*(nely+1);2*(nely+1)-1]);
%         fixednodes = fixedx/(2*(nely+1));
%         amplitude = 0;
%         fixed_dir=[ones(nelx+1,1),2*ones(nelx+1,1)];
    case 'X'
        e0 = eye(3);
        ufixed = zeros(8,3);
        U = zeros(ndof,3);
        alldofs = (1:ndof);
        n1 = [nodenrs(end,[1,end]),nodenrs(1,[end,1])];
        d1 = reshape([(2*n1-1);2*n1],1,8);
        n3 = [nodenrs(2:end-1,1)',nodenrs(end,2:end-1)];
        d3 = reshape([(2*n3-1);2*n3],1,2*(nelx+nely-2));
        n4 = [nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];
        d4 = reshape([(2*n4-1);2*n4],1,2*(nelx+nely-2));
        d2 = setdiff(alldofs,[d1,d3,d4]);
        for j = 1:3
          ufixed(3:4,j) =[e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[nelx;0];
          ufixed(7:8,j) = [e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[0;nely];
          ufixed(5:6,j) = ufixed(3:4,j)+ufixed(7:8,j);
        end
        wfixed = [repmat(ufixed(3:4,:),nely-1,1); repmat(ufixed(7:8,:),nelx-1,1)];
        qe = cell(3,3);
        Q = zeros(3,3);
        dQ = cell(3,3);
    case 'LW'
        fixednodes = union([2*(nely+1):(nely+1):10*(nely+1)],[2*(nely+1)+1:(nely+1):10*(nely+1)+1]);
        fixednodes = [1*(nely+1):(nely+1):125*(nely+1)];
        fixedy =fixednodes;
        fixedx = fixednodes(ceil(length(fixednodes)/2):end);
        fixeddofs = union([2*fixednodes, 2*fixednodes-1],[2*fixednodes+1,2*fixednodes+2])';
        fixeddofs = [2*fixednodes, 2*fixednodes-1];
        fixeddofs = union(2*fixedx-1, 2*fixedy);
        alldofs = [1:2*(nely+1)*(nelx+1)]';
%         freedofs = setdiff(alldofs,fixeddofs);
        forcedofs = [1:2:2*(nely+1)-2]';
        
        passive = zeros(nely,nelx);
        for i = nelx/2:nelx
            for j = 1:nely 
                if sqrt((j-25*nely)^2+(i-.5*nelx)^2) > 25*nely
                    passive(j,i) = 1;
%               elseif sqrt((j+25*nely)^2+(i-.5*nelx)^2) > 26*nely+1
%                   passive(j,i) = 1;
                end
            end
        end
        for i = nelx/2:nelx
            for j = 1:nely/2+1 
                if sqrt((j-25*nely)^2+(i-.5*nelx)^2) < 25*nely && sqrt((j-25*nely)^2+(i-.5*nelx)^2) > 25*nely-1 && passive(j,i)==0 && passive(j+1,i)==0
                    passive(j,i) = 2;
%               elseif sqrt((j+25*nely)^2+(i-.5*nelx)^2) < 26*nely+1 && sqrt((j+25*nely)^2+(i-.5*nelx)^2) > 26*nely && passive(j,i)==0
%                   passive(j,i) = 3;
                end
            end
        end

        for i = nelx/2:nelx
            for j = 2:nely/2+1 
                if passive(j-1,i)==2
                    passive(j,i) = 0;
%               elseif sqrt((j+25*nely)^2+(i-.5*nelx)^2) < 26*nely+1 && sqrt((j+25*nely)^2+(i-.5*nelx)^2) > 26*nely && passive(j,i)==0
%                   passive(j,i) = 3;
                end
            end
        end
%       passive(end,250:500) = 3;
        passive(1,125:250) = 2;
        [rows_up, columns_up] = find(passive==2);
        rows_down = nely-rows_up+1;
        for i=1:length(rows_down)
            passive(rows_down(i),columns_up(i)) = 3;
        end
        [rows_up_hole, columns_up_hole] = find(passive==1);
        rows_down_hole = nely-rows_up_hole+1;
        for i=1:length(rows_down_hole)
            passive(rows_down_hole(i),columns_up_hole(i)) = 1;
        end
        xPhys = 0.5*ones(nely, nelx);
        xPhys(passive==1) = 0;
        xPhys(passive==2) = 0.2;
        xPhys(passive==3) = 1;

        [rows_down, columns_down] = find(passive==3);
        elements_up = rows_up+nely*(columns_up-1);
        elements_down = rows_down+nely*(columns_down-1);
        passive(1,75:124) = 2;
        
        patchx = [299;299];
        patch1y = [50;49];
%         patch2x = [299,299];
        patch2y = [0;1];
        for i = 300:nelx
            patchx = [patchx; i];
            patch1y = [patch1y; 50-find(passive(1:25,i)==1,1,'last')];
%             patch2x = [patch2x; i];
            patch2y = [patch2y; 25-find(passive(26:50,i)==1,1)+1];
        end
        patchx = [patchx;500];
        patch1y = [patch1y;50];
%         patch2x = [299,299];
        patch2y = [patch2y;0];
        
        figure(10)
        colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow; hold off;
        
        D = p.E0/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
        B = 1/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
        DB = D*B;
        Cvm = [1 -0.5 0;-0.5 1 0;0 0 3];
        Sel = DB'*Cvm*DB;
        Sl = 0.5;
        F = sparse(2*(nely+1)*(nelx+1),1);
        
        fixednodes = [1*(nely+1):(nely+1):125*(nely+1)];
        fixedy =fixednodes;
        fixedx = fixednodes(ceil(length(fixednodes)/2):end);
        fixeddofs = union(2*fixedx-1, 2*fixedy);
        forcedofs = [1:2:2*nely]';
        displ_tip = 100;
        degree = 20;
        u_obj_fcn = @(x) displ_tip*((x-0)/(nelx+1-0)).^degree;
        figure
        fplot(u_obj_fcn,[0 nelx]); axis equal; hold off;
        u_obj = u_obj_fcn((400:nelx).');
        ddl_up = edofMat(elements_up,8);
        ddl_down = edofMat(elements_down,2);
        sel = repmat([1:length(u_obj)]',1,2);
        L = sparse(sel,[ddl_down(end-length(u_obj)+1:end) ddl_up(end-length(u_obj)+1:end)],.5,length(u_obj),2*(nelx+1)*(nely+1));
        figure
        spy(L)
        sen = 5;
        Fmin = -1;
        Fmma = -Fmin;
        eFree = find(passive==0);
        xsize = length(eFree);
    case 'RIB'
        patch1x = [rx(1:find(rx==max(rx),1,'first'));nelx;0];
        patch1y = [ry(1:find(rx==max(rx),1,'first'));nely;nely];
        patch2x = [rx(find(rx==max(rx),1,'last'):end);nelx];
        patch2y = [ry(find(rx==max(rx),1,'last'):end);0];
%         elenrs = reshape(1:nele,nely,nelx);
        in = reshape(inpolygon(xc,yc,rx,ry),nely,nelx);
%         spy(in);
        ind = find(in);
        emptyelts = setdiff(1:nele,ind);
        % Hole
        hx1 = 90; hy1 = 50; hr1 = 20; % x,y and radius of hole
        hx2 = 180; hy2 = 49; hr2 = 10;
        h1 = zeros(1,nele);
        hd1 = sqrt((xc - hx1).^2 + (yc - hy1).^2);
        hd2 = sqrt((xc - hx2).^2 + (yc - hy2).^2);
        in(hd1 <= hr1) = 0;
        in(hd2 <= hr2) = 0;
%         spy(in)
        emptyhole = [find(hd1 <= hr1) find(hd2 <= hr2)];
        emptyelts = union(emptyelts,emptyhole);
        fixednodes = union(find(in(:,1)),nelx*(nely+1)+find(in(:,end)));
        fixednodes = repmat(union(fixednodes, fixednodes+1),2,1);
        fixed_dir = [ones(numel(fixednodes)/2,1);2*ones(numel(fixednodes)/2,1)];
        fixeddofs = 2*(fixednodes-1) + fixed_dir;
        forcetop = [];
        forcebot = [];
        for i = 2:nelx-1
            forcetop = [forcetop ; (i-1)*nely+find(in(:,i),1,'first')];
            forcebot = [forcebot ; (i-1)*nely+find(in(:,i),1,'last')];
        end
        fnodetop = floor(forcetop/nely)*(nely+1)+mod(forcetop,nely);
        fnodetop = [fnodetop;fnodetop(end)+nely+1];
        fnodebot = floor((forcebot-1)/nely)*(nely+1)+mod((forcebot-1),nely)+2;
        fnodebot = [fnodebot;fnodebot(end)+nely+1];
        excitation_direction=2;
        fdoftop = fnodetop*2;
        fdofbot = fnodebot*2;
%         F = sparse([fdoftop;fdofbot],1,[-ones(size(fdoftop));...
%             ones(size(fdofbot))]/numel([fdoftop;fdofbot]),ndof,1);
%         fullelts = union(fixedel, union(forcetop,forcebot));
        F = sparse(2*([fnodetop;fnodebot]-1)+excitation_direction,1,[-ones(size(fnodetop));...
            ones(size(fnodebot))]/numel([fnodetop;fnodebot]),ndof,1);
    otherwise
        error('BC! String should be a valid entry: ''MBB'',''L-Shape'',''Short_Cantilever''')
end
mult = 1-numel(emptyelts)/nele;
volfrac = volfrac*mult;
alldofs = [1:ndof]';
freedofs = setdiff(alldofs,fixeddofs);
% 
% Compute *Gauss point coordinates and weights *in a squared sampling window 
% $[2R\times2R]$ centred in the origin

a=-R;
b=R;
[gpc,wc]=lgwt(Ngp,a,b);
[gpcx,gpcy]=meshgrid(gpc,gpc);
gauss_weight=wc*wc';
%% 
% Repeat the value ones for each element in the mesh

gpcx=reshape((repmat(gpcx(:),1,size(centroid_coordinate,1)))',[],1);
gpcy=reshape((repmat(gpcy(:),1,size(centroid_coordinate,1)))',[],1);
gauss_weight=reshape((repmat(gauss_weight(:),1,size(centroid_coordinate,1)))',[],1);
%% 
% translate the sampling window Gauss points of the element centroid coordinates

cc=repmat(centroid_coordinate,Ngp^2,1);
gauss_point=cc+[gpcx,gpcy];
%% 
% Avoid to evaluate repeated value of sampling window gauss point coordinates:

[ugp,~,idgp]=unique(gauss_point,'rows');
%% Initialize design variable vector:
% The initial design is composed of couples of crossed components regularly 
% disposed in the mesh. 

xp=linspace(min(Xx),max(Xx),ncx+2);
yp=linspace(min(Yy),max(Yy),ncy+2); 
[xx,yy]=meshgrid(xp,yp);
%% Design variable for TopX
switch BC
% Composed of elements arranged to mimic SIMP w/TopX or in box and cross
% shapes (uncomment)
    case 'X'
        t = [linspace(0,2*pi-pi/8,16)]';
%         t = [0:pi/4:2*pi-0.1]';
        Xc = [nelx/2 + cos(t)*nelx/4; 0; 0; nelx/2 - nelx/10; nelx/2 + nelx/10; nelx; nelx;...
            nelx/2 + nelx/10; nelx/2 - nelx/10];
        Yc = [nely/2 + sin(t)*nely/4; nely/2 - nely/10; nely/2 + nely/10; nely; nely; nely/2 + nely/10;...
            nely/2 - nely/10; 0; 0];
        Lc = ones(numel(Xc),1)*sqrt(nelx^2+nely^2)/2;%ones(8,1)*sqrt(nelx^2+nely^2)/2];
        Tc = [t + pi/2; 0;0; pi/2; pi/2; 0;0; pi/2; pi/2;];
%         Xc = [0;0;0;nelx/2;nelx/2;nelx;nelx;nelx;0;0;nelx/2;nelx/2;nelx/2;nelx/2;nelx;nelx;0;0;0;nelx/2;nelx/2;nelx;nelx;nelx;];
%         Yc = [zeros(8,1);nely/2*(ones(8,1));nely*(ones(8,1))];
%         Lc = [nely;sqrt(nelx^2+nely^2);nelx;nely;nelx;nelx;sqrt(nelx^2+nely^2);nely;nely;nelx;nelx;nely;sqrt(nelx^2+nely^2); ...
%         sqrt(nelx^2+nely^2);nelx;nely;nelx;nely;sqrt(nelx^2+nely^2);nelx;nely;nelx;sqrt(nelx^2+nely^2);nely];
%         Tc = [pi/2;atan2(nely,nelx);0;pi/2;0;0;atan2(-nely,nelx);pi/2;pi/2;0;0;pi/2;atan2(nely,nelx);atan2(-nely,nelx);0;pi/2;0;pi/2;atan2(nely,nelx); ...
%         0;pi/2;0;atan2(nely,nelx);pi/2];
    otherwise
        Xc=repmat(xx(:),2,1); %component center X
        Yc=repmat(yy(:),2,1); %component center Y
        Lc=2*sqrt((nelx/(ncx+2))^2+(nely/(ncy+2))^2)*ones(size(Xc)); %component length L
        Tc=atan2(nely/ncy,nelx/ncx)*[ones(length(Xc)/2,1);-ones(length(Xc)/2,1)];% component orientation angle theta
        % Tc = -pi + (2*pi)*rand(length(Xc),1);
%         Tc = [zeros(length(Xc)/2,1); pi/2*ones(length(Xc)/2,1)]; % For RIB test case only!
end
%%
hc=2*ones(length(Xc),1); % component h
Mc=initial_d*ones(size(Xc)); % component mass (For MNA and GP)
Xg=reshape([Xc,Yc,Lc,hc,Tc,Mc]',[],1);
switch BC
    case 'LW'
        Xg = [0;nely/2;nely;1.5;pi/2;1;Xg];
        dF_dX = sparse(ndof,numel(Xg)/6);
        dF_dT = dF_dX;
        dF_dh = dF_dX;
        dF_dL = dF_dX;
        dF_dY = dF_dX;
        dF_dm = dF_dX;
end
%% Build upper and lower bounds of the design problem

Xl=min(Xx-1)*ones(size(Xc));Xu=max(Xx+1)*ones(size(Xc));
Yl=min(Yy-1)*ones(size(Xc));Yu=max(Yy+1)*ones(size(Xc));
Ll=0*ones(size(Xc));Lu=sqrt(nelx^2+nely^2)*ones(size(Xc));
hl=minh*ones(size(Xc));hu=sqrt(nelx^2+nely^2)*ones(size(Xc));
Tl=-2*pi*ones(size(Xc));Tu=2*pi*ones(size(Xc));
Ml=0*ones(size(Xc));Mu=ones(size(Xc));
lower_bound=reshape([Xl,Yl,Ll,hl,Tl,Ml]',[],1);
upper_bound=reshape([Xu,Yu,Lu,hu,Tu,Mu]',[],1);
switch BC 
    case 'LW'
        lower_bound = [0;min(Yy-1);0;1.1;pi/2-0.01;0.99;lower_bound];
        upper_bound = [0;max(Yy+1);sqrt(nelx^2+nely^2);0.9;pi/2+0.01;1;upper_bound];
end
%% *Scale* the *design variable vector* accordingly  $(X\in[0,1])$:

X=(Xg-lower_bound)./(upper_bound-lower_bound);
%% MMA initialization:
pcinc = 1.03;
pen = 3;
loop = 0;
switch BC
    case 'LW'
        m = 4;
    otherwise
        m = 1;
end
n = length(X(:));
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = X(:);
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmax    = eeen;
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 2000;
kkttol  =0.001;
changetol=0.001;
kktnorm = kkttol+10;
outit = 0;
change=1;
%% 
% choose the allowable *volfrac:*
% volfrac=.4;
%% 
% Prepare plots and quantity storage:
cvec=zeros(maxoutit,1);
vvec=cvec;kvec=cvec;          %gvec=cvec;pvec=cvec;
plot_rate=5;
%initialize variables for plot
tt=0:0.005:(2*pi);tt=repmat(tt,length(Xg)/6,1);
cc=cos(tt);ss=sin(tt);
%% 
% Initialize the stopping criterion
stop_cond = 1;
%% Write logs in file for tracking command window outputs
fname = fullfile(folder_name, image_prefix);
if exist([fname,'.txt']) == 2
    f1 = fopen([fname,'.txt'],'at+');
else
    f1 = fopen([fname,'.txt'],'wt+');
end
if exist([image_prefix,'.mat']) == 2
        load([image_prefix,'.mat']);
end
storeXg = zeros(length(Xg),maxoutit);
%% Start the design loop:
while stop_cond
    outit   = outit+1;
    outeriter = outeriter+1;
    %Compute the smooth characteristic functions and gradients for each component 
    % on each sampling window Gauss point (Can support GPU)
    switch BC 
        case 'LW'
            [WF,~,dWF_dY,~,dWF_dL,~]=Wgp(Xx,Yy,Xg(1:6),p); % dimensions due to just one force component
%             F = Fmin*WF/Xg(3);
            F(1:2:end) = Fmin*WF*sin(Xg(5))/Xg(3);
            F(2:2:end) = Fmin*WF*cos(Xg(5))/Xg(3);
            dF_dY(1:2:end,1) = Fmin*dWF_dY*sin(Xg(5))/Xg(3);
            dF_dY(2:2:end,1) = Fmin*dWF_dY*cos(Xg(5))/Xg(3);
            dF_dL(1:2:end,1) = -Fmin*dWF_dL*sin(Xg(5))/Xg(3)^2;
            dF_dL(2:2:end,1) = -Fmin*dWF_dL*cos(Xg(5))/Xg(3)^2;
    end
    [W,dW_dX,dW_dY,dW_dT,dW_dL,dW_dh]=Wgp(ugp(:,1),ugp(:,2),Xg,p); % dimensions due to all component
    
    %Compute local volume fractions and gradients using generalized projection
    % delta is for densities, deltac for Young modulus
    delta=sum(reshape(W(:,idgp).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
    ddelta_dX=sum(reshape(dW_dX(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_dY=sum(reshape(dW_dY(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_dT=sum(reshape(dW_dT(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_dL=sum(reshape(dW_dL(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_dh=sum(reshape(dW_dh(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    delta_c=sum(reshape(W(:,idgp).^q.*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
    ddelta_c_dX=sum(reshape(q*dW_dX(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dY=sum(reshape(q*dW_dY(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dT=sum(reshape(q*dW_dT(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dL=sum(reshape(q*dW_dL(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dh=sum(reshape(q*dW_dh(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3)...
        ./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    % model update 
    % compute young modulus and gradients
    [E,dE,dE_dm]=model_updateM(delta_c,p,X);
    dE_dX=dE.*ddelta_c_dX;
    dE_dY=dE.*ddelta_c_dY;
    dE_dT=dE.*ddelta_c_dT;
    dE_dL=dE.*ddelta_c_dL;
    dE_dh=dE.*ddelta_c_dh;
    E=full(reshape(E(:),nely,nelx));
    %compute densities
    [rho,drho_ddelta,drho_dm]=model_updateV(delta,p,X);
    drho_dX=drho_ddelta.*ddelta_dX;
    drho_dY=drho_ddelta.*ddelta_dY;
    drho_dT=drho_ddelta.*ddelta_dT;
    drho_dL=drho_ddelta.*ddelta_dL;
    drho_dh=drho_ddelta.*ddelta_dh;
    xPhys=full(reshape(rho(:),nely,nelx));
    %Take in account passive elements
    switch BC
        case 'LW'
            xPhys(passive==1) = 0;
            xPhys(passive==2) = 0;
            xPhys(passive==3) = 1;
        otherwise
            xPhys(emptyelts) = 0;
            xPhys(fullelts) = 1;
            E(emptyelts) = p.Emin;
            E(fullelts) = p.E0;
            % FE-ANALYSIS
            sK = reshape(KE(:)*(E(:)'),64*nelx*nely,1);
            K = sparse(iK,jK,sK); K = (K+K')/2;
    end
    switch BC
        case 'LW'
            
        case 'Compliant'
            K(in,in) = K(in,in) + 0.1;
            K(out,out) = K(out,out) + 0.1;
        case 'X'
            Kr = [K(d2,d2), K(d2,d3)+K(d2,d4); K(d3,d2)+K(d4,d2), K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
    end
    switch BC
        case 'X'
            U(d1,:) = ufixed;
            U([d2,d3],:) = Kr\(-[K(d2,d1); K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4); K(d3,d4)+K(d4,d4)]*wfixed);
            U(d4,:) = U(d3,:)+wfixed;
        case 'LW'
            
        otherwise
            U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    end
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    switch BC
        case 'LW'
            [c, dc_drho, U, dc_df, K] = objective_fcn(p, xPhys, KE, iK, jK,...
                u_obj, nelx, nely, freedofs, edofMat,L, F);
        case 'Compliant'
            U1 = U(:,1); U2 = U(:,2);
            ce = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),[nely,nelx]);
            c  = U(out,1);
            dc_dE = (xPhys.^p.penalty).*ce; %penal*(E0-Emin)*xPhys.^(penal-1).*ce;
        case 'X'
            for i = 1:3
                for j = 1:3
                    U1 = U(:,i); U2 = U(:,j);
                    qe{i,j} = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx)/nele;
                    Q(i,j) = sum(sum((p.Emin+xPhys.^p.penalty*(p.E0-p.Emin)).*qe{i,j}));
                    dQ{i,j} =  xPhys.^(p.penalty).*qe{i,j};
                end
            end
            c = -(Q(1,1)+Q(2,2)+Q(1,2)+Q(2,1));
            dc_dE = -(dQ{1,1}+dQ{2,2}+dQ{1,2}+dQ{2,1});
        otherwise
            ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
            c = sum(sum((E).*ce));
            dc_dE = -ce;
    end
    switch BC
        case 'LW'
            dc_drho(emptyelts) = 0;
            dc_drho(fullelts) = 0;
            dc_dX = drho_dX*dc_drho' + [dc_df*dF_dX]';
%             dc_dX(1,1) = dc_dX(1,1) + dc_df*dF_dX;
            dc_dY = drho_dY*dc_drho' + [dc_df*dF_dY]';
            dc_dL = drho_dL*dc_drho' + [dc_df*dF_dL]'; 
            dc_dh = drho_dh*dc_drho' + [dc_df*dF_dh]';
            dc_dT = drho_dT*dc_drho' + [dc_df*dF_dT]';
            dc_dm = drho_dm*dc_drho' + [dc_df*dF_dm]';
        otherwise
            dc_dE(emptyelts) = 0;
            dc_dE(fullelts) = 0;
            dc_dX=dE_dX*dc_dE(:);
            dc_dY=dE_dY*dc_dE(:);
            dc_dL=dE_dL*dc_dE(:);
            dc_dh=dE_dh*dc_dE(:);
            dc_dT=dE_dT*dc_dE(:);
            dc_dm=dE_dm*dc_dE(:);
    end
    dc=zeros(size(X));
    dc(1:6:end)=dc_dX;
    dc(2:6:end)=dc_dY;
    dc(3:6:end)=dc_dL;
    dc(4:6:end)=dc_dh;
    dc(5:6:end)=dc_dT;
    dc(6:6:end)=dc_dm;
    switch BC
        case 'LW'
            [ v, dv_drho, dv_dF, dv_dXg ] = constrain_fcn(p, xPhys, volfrac,...
                nelx, nely, U, edofMat, KE, freedofs, Sel, Sl, iK, jK, K, F, Xg);
            dv_dX = drho_dX*dv_drho + dF_dX'*dv_dF;
            dv_dY = drho_dY*dv_drho + dF_dY'*dv_dF;
            dv_dL = drho_dL*dv_drho + dF_dL'*dv_dF;
            dv_dh = drho_dh*dv_drho + dF_dh'*dv_dF;
            dv_dT = drho_dT*dv_drho + dF_dT'*dv_dF;
            dv_dm = drho_dm*dv_drho + dF_dm'*dv_dF;
            dv = zeros(n,m);
        otherwise
            v=mean(xPhys(:));
            dv_dxPhys = ones(nely,nelx)/nele;
            dv_dxPhys(emptyelts) = 0;
            dv_dxPhys(fullelts) = 0;
            dv_dX=drho_dX*dv_dxPhys(:);
            dv_dY=drho_dY*dv_dxPhys(:);
            dv_dL=drho_dL*dv_dxPhys(:);
            dv_dh=drho_dh*dv_dxPhys(:);
            dv_dT=drho_dT*dv_dxPhys(:);
            dv_dm=drho_dm*dv_dxPhys(:);
            dv=zeros(size(X));
    end
    dv(1:6:end,:)=dv_dX;
    dv(2:6:end,:)=dv_dY;
    dv(3:6:end,:)=dv_dL;
    dv(4:6:end,:)=dv_dh;
    dv(5:6:end,:)=dv_dT;
    dv(6:6:end,:)=dv_dm;
    % store the output for plot
    cvec(outit)=c;vvec(outit)=v(1);
    %% Greyness Level
    gl = 4/nele*sum(xPhys(:).*(1-xPhys(:)));
    %% PRINT RESULTS
    fprintf(f1,' It.:%5i Obj.:%4.3e Vol.:%7.3f kktnorm.:%7.3f ch.:%7.3f GL.: %5.3f T: %10.3f \n',outit-1,full(c), ...
        mean(xPhys(:)),kktnorm,change,gl, toc);
    storeXg(:,outit) = Xg;
    % pass scaled objective and constraint function and sensitivities to MMA
    f0val=log(c+1);
    switch BC
        case 'LW'
            dv = dv + dv_dXg;
            fval = v';
            dfdx = dv'.*repmat((upper_bound(:)-lower_bound(:))',m,1);
        otherwise
            fval=[(v-volfrac)/volfrac]*100;
            dfdx = [dv(:)'/volfrac]*100.*(upper_bound(:)-lower_bound(:))';
    end
    df0dx = (dc(:)/(c+1).*(upper_bound(:)-lower_bound(:)));
    %plot every plot_rate iterations
    if outit > 101
        plot_rate = 20;
    end
    if or(mod(outit,plot_rate) == 0,outit == 1)
        %convergence plot
        figure(3)
        subplot(2,1,1)
        plot(1:outit,cvec(1:outit),'bo','MarkerFaceColor','b')
        grid on
        hold on
        scatter(outit,c,'k','fill')
        hold off
        text(outit,c,['C =',num2str(c,'%4.2f'),' at iteration ', num2str(outit)],...
            'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
         xlabel('iter')
        ylabel('C')
        subplot(2,1,2)
        plot(1:outit,vvec(1:outit)*100,'ro','MarkerFaceColor','r')
        grid on
        hold on
        scatter(outit,mean(xPhys(:))*100,'k','fill')
        hold off
        text(outit,mean(xPhys(:))*100,['V = ',num2str(mean(xPhys(:))*100,'%4.2f'),'% at iteration ', num2str(outit)],...
            'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
        xlabel('iter')
        ylabel('V [%]')
        print([Path,image_prefix,'convergence'],'-dpng')
        % KKT PLOT
        figure(4)
        plot(1:outit,kvec(1:outit),'bo','MarkerFaceColor','b');
        grid on; hold on;
        scatter(outit,c,'k','fill')
        hold off
        xlabel('iter')
        ylabel('kktnorm')
        %% PLOT DENSITIES
        figure(1)
        map=colormap(gray);
        map=map(end:-1:1,:);
        caxis([0 1])
        switch BC
            case 'RIB'
                patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],...
            'FaceColor','flat','EdgeColor','none');
            otherwise
                patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],...
            'FaceColor','flat','EdgeColor','none');
        end
        axis equal; axis off;
        hold on
        fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
        switch BC
            case 'LW'
%             case 'Compliant'
            otherwise
                scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
                scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
                scal=10;
        end
%         quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
%             -(excitation_direction==2), scal,'r','Linewidth',2)
        colormap(map)
        colorbar
        drawnow
        hold off
        axis([min(Xx),max(Xx),min(Yy),max(Yy)])
        print([Path,'density_',num2str(outit-1,'%04d')],'-dpng')
        %% Component Plot
        figure(2)
        Xc=Xg(1:6:end);
        Yc=Xg(2:6:end);
        Lc=Xg(3:6:end);
        hc=Xg(4:6:end);
        Tc=Xg(5:6:end);
        Mc=Xg(6:6:end);
        C0=repmat(cos(Tc),1,size(cc,2));S0=repmat(sin(Tc),1,size(cc,2));
        xxx=repmat(Xc(:),1,size(cc,2))+cc;
        yyy=repmat(Yc(:),1,size(cc,2))+ss;
        xi=C0.*(xxx-Xc)+S0.*(yyy-Yc);
        Eta=-S0.*(xxx-Xc)+C0.*(yyy-Yc);
        [dd]=norato_bar(xi,Eta,repmat(Lc(:),1,size(cc,2)),repmat(hc(:),1,size(cc,2)));
        xn=repmat(Xc,1,size(cc,2))+dd.*cc;
        yn=repmat(Yc,1,size(cc,2))+dd.*ss;
        tolshow=0.1;
        Shown_compo=find(Mc>tolshow);
        fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
        hold on
        fill(xn(Shown_compo,:)',yn(Shown_compo,:)',Mc(Shown_compo),'FaceAlpha',0.5)
        if strcmp(BC,'L_Shape')
            fill([fix((min(Xx)+max(Xx))/2),max(Xx),max(Xx),fix((min(Xx)+max(Xx))/2)],[fix((min(Yy)+max(Yy))/2),...
                fix((min(Yy)+max(Yy))/2),max(Yy),max(Yy)],'w')
        end
        caxis([0,1])
        colormap 'jet'
        switch BC
            case 'LW'
                patch(patchx, patch1y, 'w')
                patch(patchx, patch2y, 'w')
            case 'RIB'
                patch(patch1x, patch1y, 'w')
                patch(patch2x, patch2y, 'w')
                rectangle('Position',[hx1-hr1 hy1-hr1 2*hr1 2*hr1],'Curvature',[1 1],'FaceColor','w')
                rectangle('Position',[hx2-hr2 hy2-hr2 2*hr2 2*hr2],'Curvature',[1 1],'FaceColor','w')
        end
        axis equal; axis off;
        hold on
%         scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
%         scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
        scal=10;
%         quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
%             -(excitation_direction==2),scal,'r','Linewidth',2)
        colorbar
        axis([min(Xx),max(Xx),min(Yy),max(Yy)])
        hold off
        print([Path,'component_',num2str(outit-1,'%04d')],'-dpng')        
    end
    %% MMA code optimization
        [X,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
            mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
            f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
        xold2 = xold1;
        xold1 = xval;
        xval  = X;
        Xg=lower_bound+(upper_bound-lower_bound).*X;
        change=norm(xval-xold1);
        %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,d);
    % update the stopping criterion
    switch stopping_criteria
        case 'kktnorm'
            stop_cond = (outit < maxoutit && kktnorm>kkttol);
            kvec(outit) = kktnorm;
        case 'change'
            stop_cond = (outit < maxoutit &&change>changetol);
    end
    if toc/3600 > 23
        save([image_prefix,'.mat']);
    end
end
% Make the plot of the solution
% convergence plot
figure(3)
subplot(2,1,1)
plot(1:outit,cvec(1:outit),'bo','MarkerFaceColor','b')
grid on
hold on
scatter(outit,c,'k','fill')
hold off
text(outit,c,['C =',num2str(c,'%4.2f'),' at iteration ', num2str(outit)],...
    'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
xlabel('iter')
ylabel('C')
subplot(2,1,2)
plot(1:outit,vvec(1:outit)*100,'ro','MarkerFaceColor','r')
grid on
hold on
scatter(outit,mean(xPhys(:))*100,'k','fill')
hold off
text(outit,mean(xPhys(:))*100,['V = ',num2str(mean(xPhys(:))*100,'%4.2f'),'% at iteration ', num2str(outit)],...
    'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
xlabel('iter')
ylabel('V [%]')
print([Path,image_prefix,'convergence'],'-dpng')
%% PLOT DENSITIES
figure(1)
map=colormap(gray);
map=map(end:-1:1,:);
caxis([0 1])
patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],...
    'FaceColor','flat','EdgeColor','none'); axis equal; axis off; hold on
hold on
fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
% scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
% scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
% quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
%     -(excitation_direction==2), scal,'r','Linewidth',2)
colormap(map)
colorbar
drawnow
hold off
axis([min(Xx),max(Xx),min(Yy),max(Yy)])
print([Path,'density_',num2str(outit-1,'%03d')],'-dpng')
%% Component Plot
figure(2)
Xc=Xg(1:6:end);
Yc=Xg(2:6:end);
Lc=Xg(3:6:end);
hc=Xg(4:6:end);
Tc=Xg(5:6:end) ;
Mc=Xg(6:6:end) ;
C0=repmat(cos(Tc),1,size(cc,2));S0=repmat(sin(Tc),1,size(cc,2));
xxx=repmat(Xc(:),1,size(cc,2))+cc;
yyy=repmat(Yc(:),1,size(cc,2))+ss;
xi=C0.*(xxx-Xc)+S0.*(yyy-Yc);
Eta=-S0.*(xxx-Xc)+C0.*(yyy-Yc);
[dd]=norato_bar(xi,Eta,repmat(Lc(:),1,size(cc,2)),repmat(hc(:),1,size(cc,2)));
xn=repmat(Xc,1,size(cc,2))+dd.*cc;
yn=repmat(Yc,1,size(cc,2))+dd.*ss;
tolshow=0.1;
Shown_compo=find(Mc>tolshow);
fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
hold on
fill(xn(Shown_compo,:)',yn(Shown_compo,:)',Mc(Shown_compo),'FaceAlpha',0.5)
if strcmp(BC,'L-shape')
    fill([fix((min(Xx)+max(Xx))/2),max(Xx),max(Xx),fix((min(Xx)+max(Xx))/2)],...
        [fix((min(Yy)+max(Yy))/2),fix((min(Yy)+max(Yy))/2),max(Yy),max(Yy)],'w')
end
caxis([0,1])
colormap 'jet'
axis equal; axis off;
hold on
% scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
% scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
% quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
%     -(excitation_direction==2),scal,'r','Linewidth',2)
colorbar
axis([min(Xx),max(Xx),min(Yy),max(Yy)])
print([Path,'component_',num2str(outit-1,'%03d')],'-dpng')
hold off
fprintf(f1,'Elapsed Time: %10.3f',toc);
fclose(f1);
matname = fullfile(folder_name, 'desvar.mat');
save(matname,'storeXg','xPhys');
end
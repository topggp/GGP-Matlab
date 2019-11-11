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
% In this section of the Matlab code we define several *parameters* needed for 
% the *Generalized Geometry Projection*.

% GGP parameters
stopping_criteria='change'; %stopping criteria of the optimization algorithm either change or KKT norm
nelx=60;nely=30; 
BC='Short_Cantilever';%L-shape %Short_Cantilever%MBB
p.method='GP';%MMC%MNA %GP this change the function employed for the evaluation of local volume fraction
q=1;%q=1
p.zp=1 ;% parameter for p-norm/mean regularization
p.alp=1; %parameter for MMC
p.epsi=0.866;% parameter for MMC
p.bet=1e-3; %parameter for MMC
p.deltamin=1e-6; %parameter for GP
p.r=.5;%parameter for GP
minh=1;% minimal bar thickness
p.sigma=1;%parameter for MNA
p.gammav=1;%parameter for GP
p.gammac=3;%parameter for GP
p.penalty=3;%parameter for MNA
p.aggregation='KSl'; %parameter for the aggregation function to be used
% IE= Induced Exponential % KS= KS function %KSl= lowerbound KS function
% p-norm %p-mean
p.ka=10; % parameter for the aggregation constant
p.saturation=true; % switch for saturation
ncx=1; % number of components in the x direction
ncy=1; % number of components in the y direction
Ngp=2; % number of Gauss point per sampling window
R=0.5; % radius of the sampling window (infty norm)
initial_d=0.5; % initial mass variable adopted for MNA and GP
%% 
% Generate a *folder* and a prefix to save images *optimization history*:

rs=replace(num2str(R,'%3.2f'),'.','_');
folder_name=['Optimization_history_',BC,p.method,'nelx_',num2str(nelx),...
    'nely_',num2str(nely),'_R_',rs,'_Ngp_',num2str(Ngp),'_SC_',stopping_criteria];
image_prefix=[BC,p.method,'nelx_',num2str(nelx),'nely_',num2str(nely),'_R_',rs,'_Ngp_',num2str(Ngp)];
mkdir(folder_name)
Path=[folder_name,'/'];
%% 
% Define *Material properties*:

% MATERIAL PROPERTIES
p.E0 = 1;
p.Emin = 1e-6;
nu = 0.3;
%% 
% Prepare *finite element analysis*

% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
U = zeros(2*(nely+1)*(nelx+1),1);
%define the nodal coordinates
[Yy,Xx]=find(nodenrs);
Yy=nely+1-Yy;
Xx=Xx-1;
% Element connectivity
enodeMat=edofMat(:,[2,4,6,8])/2;
% DEFINE LOADS AND SUPPORTS 
switch BC
    case 'MBB'
        excitation_node=1;excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=[find(Xx==min(Xx));(nelx+1)*(nely+1)];fixed_dir=[ones(nely+1,1);2];
        fixeddofs=2*(fixednodes-1)+fixed_dir;
        emptyelts=[]; fullelts = [];
    case 'Short_Cantilever'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=repmat(find(Xx==min(Xx)),2,1);fixed_dir=[ones(nely+1,1);2*ones(nely+1,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
        emptyelts=[]; fullelts = [];
    case 'L-shape'
        excitation_node=find((Xx==max(Xx))&(Yy==fix(0.5*min(Yy)+0.5*max(Yy))));excitation_direction=2;
        amplitude=-1;
        F = sparse(2*(excitation_node-1)+excitation_direction,1,amplitude,2*(nely+1)*(nelx+1),1);
        fixednodes=repmat(find(Yy==max(Yy)),2,1);fixed_dir=[ones(nelx+1,1),2*ones(nelx+1,1)];
        fixeddofs=2*(fixednodes-1)+fixed_dir(:);
        emptyelts=find(xc>=(((max(Xx)+min(Xx))/2))&(yc>=((max(Yy)+min(Yy))/2)));
        fullelts = [];
    otherwise
        error('BC string should be a valid entry: ''MBB'',''L-Shape'',''Short_Cantilever''')
end
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% Prepare the *Generalized Geometry Projection:*
% Compute the element centroid coordinates:

xc=mean(Xx(enodeMat'));
yc=mean(Yy(enodeMat'));
centroid_coordinate=[xc(:),yc(:)];
%% 
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
Xc=repmat(xx(:),2,1); %component center X
Yc=repmat(yy(:),2,1); %component center Y
Lc=2*sqrt((nelx/(ncx+2))^2+(nely/(ncy+2))^2)*ones(size(Xc)); %component length L
Tc=atan2(nely/ncy,nelx/ncx)*[ones(length(Xc)/2,1);-ones(length(Xc)/2,1)];% component orientation angle tetha
hc=2*ones(length(Xc),1); % component h
Mc=initial_d*ones(size(Xc)); % component mass (For MNA and GP)
Xg=reshape([Xc,Yc,Lc,hc,Tc,Mc]',[],1);
%% Build upper and lower bounds of the design problem

Xl=min(Xx-1)*ones(size(Xc));Xu=max(Xx+1)*ones(size(Xc));
Yl=min(Yy-1)*ones(size(Xc));Yu=max(Yy+1)*ones(size(Xc));
Ll=0*ones(size(Xc));Lu=sqrt(nelx^2+nely^2)*ones(size(Xc));
hl=minh*ones(size(Xc));hu=sqrt(nelx^2+nely^2)*ones(size(Xc));
Tl=-2*pi*ones(size(Xc));Tu=2*pi*ones(size(Xc));
Ml=0*ones(size(Xc));Mu=ones(size(Xc));
lower_bound=reshape([Xl,Yl,Ll,hl,Tl,Ml]',[],1);
upper_bound=reshape([Xu,Yu,Lu,hu,Tu,Mu]',[],1);
%% 
% *Scale* the *design variable vector* accordingly  $(X\in[0,1])$:

X=(Xg-lower_bound)./(upper_bound-lower_bound);
%% MMA initialization:

loop = 0;
m = 1;
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
volfrac=.4;
%% 
% Prepare plots and quantity storage:
cvec=zeros(maxoutit,1);
vvec=cvec;ovvec=cvec;gvec=cvec;pvec=cvec;
plot_rate=10;
%initialize variables for plot
tt=0:0.005:(2*pi);tt=repmat(tt,length(Xc),1);
cc=cos(tt);ss=sin(tt);
%% 
% Initialize the stopping criterion
switch stopping_criteria
    case 'kktnorm'
        stop_cond=outit < maxoutit && kktnorm>kkttol;
    case 'change'
        stop_cond=outit < maxoutit &&change>changetol;
end
%% Start the design loop:
while  stop_cond
    outit   = outit+1;
    outeriter = outeriter+1;
    %Compute the smooth characteristic functions and gradients for each component 
    % on each sampling window Gauss point (Can support GPU)
    [W,dW_dX,dW_dY,dW_dT,dW_dL,dW_dh]=Wgp(ugp(:,1),ugp(:,2),Xg,p);
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
    xPhys(emptyelts) = 0;
    xPhys(fullelts) = 1;
    E(emptyelts) = p.Emin;
    E(fullelts) = p.E0;
    % FE-ANALYSIS
    sK = reshape(KE(:)*(E(:)'),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((E).*ce));
    v=mean(xPhys(:));
    dc_dE = -ce;
    dc_dE(emptyelts) = 0;
    dc_dE(fullelts) = 0;
    dc_dX=dE_dX*dc_dE(:);
    dc_dY=dE_dY*dc_dE(:);
    dc_dL=dE_dL*dc_dE(:);
    dc_dh=dE_dh*dc_dE(:);
    dc_dT=dE_dT*dc_dE(:);
    dc_dm=dE_dm*dc_dE(:);
    dc=zeros(size(X));
    dc(1:6:end)=dc_dX;
    dc(2:6:end)=dc_dY;
    dc(3:6:end)=dc_dL;
    dc(4:6:end)=dc_dh;
    dc(5:6:end)=dc_dT;
    dc(6:6:end)=dc_dm;
    dv_dxPhys = ones(nely,nelx)/nelx/nely;
    dv_dxPhys(emptyelts) = 0;
    dv_dxPhys(fullelts) = 0;
    dv_dX=drho_dX*dv_dxPhys(:);
    dv_dY=drho_dY*dv_dxPhys(:);
    dv_dL=drho_dL*dv_dxPhys(:);
    dv_dh=drho_dh*dv_dxPhys(:);
    dv_dT=drho_dT*dv_dxPhys(:);
    dv_dm=drho_dm*dv_dxPhys(:);
    dv=zeros(size(X));
    dv(1:6:end)=dv_dX;
    dv(2:6:end)=dv_dY;
    dv(3:6:end)=dv_dL;
    dv(4:6:end)=dv_dh;
    dv(5:6:end)=dv_dT;
    dv(6:6:end)=dv_dm;
    % store the output for plot
    cvec(outit)=c;vvec(outit)=v;
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%4.3e Vol.:%7.3f kktnorm.:%7.3f ch.:%7.3f\n',outit,c, ...
        mean(xPhys(:)),kktnorm,change);
    % pass scaled objective and constraint function and sensitivities to MMA
    f0val=log(c+1);
    fval=[(v-volfrac)/volfrac]*100;
    df0dx=(dc(:)/(c+1).*(upper_bound(:)-lower_bound(:)));
    dfdx=[dv(:)'/volfrac]*100.*(upper_bound(:)-lower_bound(:))';
    %plot every plot_rate iterations
    if rem(outit,plot_rate)==0
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
        %% PLOT DENSITIES
        figure(1)
        map=colormap(gray);
        map=map(end:-1:1,:);
        caxis([0 1])
        patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],...
            'FaceColor','flat','EdgeColor','none'); axis equal; axis off; hold on
        hold on
        fill([min(Xx),max(Xx),max(Xx),min(Xx)],[min(Yy),min(Yy),max(Yy),max(Yy)],'w','FaceAlpha',0.)
        scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
        scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
        scal=10;
        quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
            -(excitation_direction==2), scal,'r','Linewidth',2)
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
            fill([fix((min(Xx)+max(Xx))/2),max(Xx),max(Xx),fix((min(Xx)+max(Xx))/2)],[fix((min(Yy)+max(Yy))/2),...
                fix((min(Yy)+max(Yy))/2),max(Yy),max(Yy)],'w')
        end
        caxis([0,1])
        colormap 'jet'
        axis equal; axis off;
        hold on
        scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
        scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
        scal=10;
        quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
            -(excitation_direction==2),scal,'r','Linewidth',2)
        colorbar
        axis([min(Xx),max(Xx),min(Yy),max(Yy)])
        print([Path,'component_',num2str(outit-1,'%03d')],'-dpng')
        hold off
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
            stop_cond=outit < maxoutit && kktnorm>kkttol;
        case 'change'
            stop_cond=outit < maxoutit &&change>changetol;
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
scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
    -(excitation_direction==2), scal,'r','Linewidth',2)
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
scatter(Xx(fixednodes(fixed_dir==1)),Yy(fixednodes(fixed_dir==1)),'>b','filled')
scatter(Xx(fixednodes(fixed_dir==2)),Yy(fixednodes(fixed_dir==2)),'^b','filled')
scal=10;
quiver(Xx(excitation_node),Yy(excitation_node)+scal*(excitation_direction==2),excitation_direction==1,...
    -(excitation_direction==2),scal,'r','Linewidth',2)
colorbar
axis([min(Xx),max(Xx),min(Yy),max(Yy)])
print([Path,'component_',num2str(outit-1,'%03d')],'-dpng')
hold off
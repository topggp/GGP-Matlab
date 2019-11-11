%-------------------------------------------------------------
%    This is the file Wgp.m
%
%    Version Nov 2019.
%    Simone Coniglio <simone.coniglio@airbus.com>
%    Propulsion Airframe Stress Transverse,
%    31300 Toulouse, France.
%
function [W,dW_dX,dW_dY,dW_dT,dW_dL,dW_dh]=Wgp(x,y,Xc,p)
%  Evaluate characteristic function in each Gauss point
ii=1:numel(x);
X=Xc(1:6:end);
Y=Xc(2:6:end);
L=Xc(3:6:end);
h=Xc(4:6:end);
T=Xc(5:6:end);
jj=1:numel(X);
[I,J]=meshgrid(ii,jj);
xi=reshape(x(I),size(I));
yi=reshape(y(I),size(I));
rho=sqrt((X(J)-xi).^2+(Y(J)-yi).^2);
drho_dX=(X(J)-xi)./(rho+(rho==0));
drho_dY=(Y(J)-yi)./(rho+(rho==0));
phi=atan2(-Y(J)+yi,-(X(J)-xi))-T(J);
dphi_dX=((-Y(J)+yi)./(rho.^2+(rho==0)));
dphi_dY=(X(J)-xi)./(rho.^2+(rho==0));
dphi_dT=-ones(size(J));
upsi=sqrt(rho.^2+L(J).^2/4-rho.*L(J).*abs(cos(phi))).*(((rho.*cos(phi)).^2)>=(L(J).^2/4))+~(((rho.*cos(phi)).^2)>=(L(J).^2/4)).*abs(rho.*sin(phi));
dupsi_drho=(2*rho-L(J).*abs(cos(phi)))/2./(upsi+(upsi==0)).*((((rho.*cos(phi)).^2)>=(L(J).^2/4)))+~(((rho.*cos(phi)).^2)>=(L(J).^2/4)).*abs(sin(phi));
dupsi_dphi=(L(J).*rho.*sign(cos(phi)).*sin(phi))/2./(upsi+(upsi==0)).*((((rho.*cos(phi)).^2)>=(L(J).^2/4)))+~(((rho.*cos(phi)).^2)>=(L(J).^2/4)).*rho.*sign(sin(phi)).*cos(phi);
dupsi_dL=(L(J)/2-rho.*abs(cos(phi)))./2./(upsi+(upsi==0)).*((((rho.*cos(phi)).^2)>=(L(J).^2/4))&upsi~=0);
switch p.method
    case 'MMC'
        alp=p.alp;
        epsi=p.epsi;
        bet=p.bet;
        chi0=1-(4*upsi.^2./h(J).^2).^alp;
        dchi0_dh=8*alp*upsi.^2.*(4*upsi.^2./h(J).^2).^(alp-1)./h.^3;
        dchi0_dupsi=-8*alp*upsi.*(4*upsi.^2./h(J).^2).^(alp-1)./h.^2;
        [chi,dchi]=Aggregation_Pi(chi0,p);
        dchi_dh=(dchi0_dh.*dchi);
        dchi_dupsi=(dchi0_dupsi.*dchi);
        chi(chi<=-1e6)=-1e6;
        W=(chi>epsi)+(chi<=epsi&chi>=-epsi).*(3/4*(1-bet)*(chi/epsi-chi.^3/3/epsi^3)+(1+bet)/2)+(chi<-epsi)*bet;
        dW_dchi=-3/4*(1/epsi-chi.^2/epsi^3).*(bet-1).*(abs(chi)<epsi);
        dW_dupsi=repmat(dW_dchi,size(dchi_dh,1),1).*dchi_dupsi;
        dW_dh=repmat(dW_dchi,size(dchi_dh,1),1).*dchi_dh;
        dW_dX=dW_dupsi.*(dupsi_dphi.*dphi_dX+dupsi_drho.*drho_dX);
        dW_dY=dW_dupsi.*(dupsi_dphi.*dphi_dY+dupsi_drho.*drho_dY);
        dW_dL=dW_dupsi.*dupsi_dL;
        dW_dT=dW_dupsi.*dupsi_dphi.*dphi_dT;
    case 'GP'
        deltamin=p.deltamin;
        r=p.r;
        zetavar=upsi-h(J)/2;
        dzetavar_dupsi=ones(size(upsi));
        dzetavar_dh=-0.5*ones(size(J));
        deltaiel=(1/pi/r^2*(r^2*acos(zetavar/r)-zetavar.*sqrt(r^2-zetavar.^2))).*(abs(zetavar)<=r)+((zetavar<-r));
        ddetlaiel_dzetavar=(-2*sqrt(r^2-zetavar.^2)/pi/r^2).*(abs(zetavar)<=r);
        W=deltamin+(1-deltamin)*deltaiel;
        dW_ddeltaiel=(1-deltamin);
        dW_dh=dW_ddeltaiel*ddetlaiel_dzetavar.*dzetavar_dh;
        dW_dupsi=dW_ddeltaiel*ddetlaiel_dzetavar.*dzetavar_dupsi;
        dW_dX=dW_dupsi.*(dupsi_dphi.*dphi_dX+dupsi_drho.*drho_dX);
        dW_dY=dW_dupsi.*(dupsi_dphi.*dphi_dY+dupsi_drho.*drho_dY);
        dW_dL=dW_dupsi.*dupsi_dL;
        dW_dT=dW_dupsi.*dupsi_dphi.*dphi_dT;
    case 'MNA'
        epsi=p.sigma;
        ds=upsi;
        d=abs(upsi);
        l=h(J)/2-epsi/2;
        u=h(J)/2+epsi/2;
        a3= -2./((l - u).*(l.^2 - 2*l.*u + u.^2));
        a2=   (3*(l + u))./((l - u).*(l.^2 - 2*l.*u + u.^2));
        a1=    -(6*l.*u)./((l - u).*(l.^2 - 2*l.*u + u.^2));
        a0=(u.*(- u.^2 + 3*l.*u))./((l - u).*(l.^2 - 2*l.*u + u.^2));
        W=1*(d<=l)+(a3.*d.^3+a2.*d.^2+a1.*d+a0).*(d<=u&d>l);
        dW_dupsi=sign(ds).*(3*a3.*d.^2+2*a2.*d+a1).*(d<=u&d>l);
        da3_du=- 2./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) - (2*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        da2_du=3./((l - u).*(l.^2 - 2*l.*u + u.^2)) + (3*(l + u))./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) + (3*(l + u).*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        da1_du=- (6*l)./((l - u).*(l.^2 - 2*l.*u + u.^2)) - (6*l.*u)./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) - (6*l.*u.*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        da0_du=(- u.^2 + 3*l.*u)./((l - u).*(l.^2 - 2*l.*u + u.^2)) + (u.*(- u.^2 + 3*l.*u))./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) + (u.*(3*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2)) + (u.*(- u.^2 + 3*l.*u).*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        dWf_du=(da3_du.*d.^3+da2_du.*d.^2+da1_du.*d+da0_du).*(d<=u&d>l);
        da3_dl=2./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) + (2*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        da2_dl= 3./((l - u).*(l.^2 - 2*l.*u + u.^2)) - (3*(l + u))./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) - (3*(l + u).*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        da1_dl=    (6*l.*u)./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) - (6*u)./((l - u).*(l.^2 - 2*l.*u + u.^2)) + (6*l.*u.*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        da0_dl= (3*u.^2)./((l - u).*(l.^2 - 2*l.*u + u.^2)) - (u.*(- u.^2 + 3*l.*u))./((l - u).^2.*(l.^2 - 2*l.*u + u.^2)) - (u.*(- u.^2 + 3*l.*u).*(2*l - 2*u))./((l - u).*(l.^2 - 2*l.*u + u.^2).^2);
        dWf_dl=(da3_dl.*d.^3+da2_dl.*d.^2+da1_dl.*d+da0_dl).*(d<=u&d>l);
        dW_dh=0.5*sign(ds).*(dWf_du+dWf_dl);
        dW_dX=dW_dupsi.*(dupsi_dphi.*dphi_dX+dupsi_drho.*drho_dX);
        dW_dY=dW_dupsi.*(dupsi_dphi.*dphi_dY+dupsi_drho.*drho_dY);
        dW_dL=dW_dupsi.*dupsi_dL;
        dW_dT=dW_dupsi.*dupsi_dphi.*dphi_dT;
end
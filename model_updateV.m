function [rho,drho_ddelta,drho_dm]=model_updateV(delta,p,X)
% update local densities
m=X(6:6:end);
nc=length(m);
m=repmat(m(:),1,size(delta,2));
switch p.method
    case 'MMC'
        rho=delta;
        drho_ddelta=ones(size(delta));
        drho_ddelta=repmat(drho_ddelta,size(m,1),1);
        drho_dm=0*m;
    case 'GP'
        hatdelta=delta.*m.^p.gammav;
        [rho,drho_dhatdelta]=Aggregation_Pi(hatdelta,p);
        if p.saturation
        [rho,ds]=smooth_sat(rho,p,nc);
        drho_dhatdelta=ds.*drho_dhatdelta;
        end
        dhatdelta_ddelta=m.^p.gammav;
        dhatdelta_dm=p.gammav*delta.*m.^(p.gammav-1);
        drho_ddelta=dhatdelta_ddelta.*drho_dhatdelta;
        drho_dm=drho_dhatdelta.*dhatdelta_dm;
    case 'MNA'
        hatdelta=delta.*m.^p.gammav;
        [rho,drho_dhatdelta]=Aggregation_Pi(hatdelta,p);
        if p.saturation
        [rho,ds]=smooth_sat(rho,p,nc);
        drho_dhatdelta=ds.*drho_dhatdelta;
        end
        dhatdelta_ddelta=m.^p.gammav;
        dhatdelta_dm=p.gammav*delta.*m.^(p.gammav-1);
        drho_ddelta=dhatdelta_ddelta.*drho_dhatdelta;
        drho_dm=drho_dhatdelta.*dhatdelta_dm;
end

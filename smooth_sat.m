function [s,ds]=smooth_sat(y,p,nc)
switch p.aggregation
     case 'p-norm'
       xt=1;
    case 'p-mean'
       xt=(((nc-1)*p.zp^p.ka+(1+p.zp)^p.ka)/nc)^(1/p.ka)-p.zp;
    case 'KS'
         xt=1;
    case 'KSl'
         xt=1+1/p.ka*log((1+(nc-1)*exp(-p.ka))/nc);
    case 'IE'
        xt=1+1/p.ka*log((1+(nc-1)*exp(-p.ka))/nc);
end
pp=100;
s0=-log(exp(-pp)+1.0./(exp((pp.*0)./xt)+1.0))./pp;
s=@(xs,a,pa)(-log(exp(-pa)+1.0./(exp((pa.*xs)./a)+1.0))./pa-s0)/(1-s0);
ds=@(xs,a,pa)(exp((pa.*xs)./a).*1.0./(exp((pa.*xs)./a)+1.0).^2)./(a.*(exp(-pa)+1.0./(exp((pa.*xs)./a)+1.0)))/(1-s0);
% syms a xs
s=s(y,xt,pp);
ds=ds(y,xt,pp);
% s=((xt-2)/xt^3*y.^3+(3-2*xt)/xt^2*y.^2+y).*(y<=xt)+~(y<=xt);
% ds=(3*(xt-2)/xt^3*y.^2+2*(3-2*xt)/xt^2*y+ones(size(y))).*(y<=xt);
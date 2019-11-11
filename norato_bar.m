%-------------------------------------------------------------
%    This is the file norato_bar.m
%
%    Version Nov 2019.
%    Simone Coniglio <simone.coniglio@airbus.com>
%    Propulsion Airframe Stress Transverse,
%    31300 Toulouse, France.
%
function [d]=norato_bar(xi,eta,L,h)
% to be used for plot with fill function 
d=((L/2.*sqrt(xi.^2./(xi.^2+eta.^2))+sqrt(h.^2/4-L.^2/4.*eta.^2./(xi.^2+eta.^2)))...
    .*(xi.^2./(xi.^2+eta.^2)>=(L.^2./(h.^2+L.^2)))+h./2.*sqrt(1+xi.^2./(eta.^2+(eta==0)))...
    .*(~(xi.^2./(xi.^2+eta.^2)>=(L.^2./(h.^2+L.^2)))))...
    .*(xi~=0|eta~=0)+sqrt(2)/2*h.*(xi==0&eta==0);

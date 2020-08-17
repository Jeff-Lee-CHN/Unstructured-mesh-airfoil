
function [NablaT_fluxface]= gardientOnFluxface(  NablaT_node,  fluxface,T_old,Length_edge, tau )
    
NablaT_ave(:,1)=( NablaT_node(fluxface(:,1),1)+NablaT_node(fluxface(:,2),1) )/2;
NablaT_ave(:,2)=( NablaT_node(fluxface(:,1),2)+NablaT_node(fluxface(:,2),2) )/2;

dT_dl=( T_old(fluxface(:,2),1)-T_old(fluxface(:,1),1) )./ Length_edge;

NablaT_fluxface(:,1)=NablaT_ave(:,1)-(  NablaT_ave(:,1).*tau(:,1)+ NablaT_ave(:,2).*tau(:,2) -dT_dl ).*tau(:,1);
NablaT_fluxface(:,2)=NablaT_ave(:,2)-(  NablaT_ave(:,1).*tau(:,1)+ NablaT_ave(:,2).*tau(:,2) -dT_dl ).*tau(:,2);
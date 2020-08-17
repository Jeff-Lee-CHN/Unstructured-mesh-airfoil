function   [ rho_u ]=convectionOnFluxface( fluxfaceList,Ux,Uy,varargin  )
[N,~]=size(fluxfaceList);
 rho_u1=ones(N,1);
 rho_u2=ones(N,1); 
 rho_u=zeros(N,2);
 Ngetin=nargin-3;
for i=1: Ngetin
 rho_u1= rho_u1.* varargin{1,i}(fluxfaceList(:,1),1);
 rho_u2= rho_u2.*varargin{1,i}(fluxfaceList(:,2),1);
end
rho_u(:,1)=( rho_u1.*Ux(fluxfaceList(:,1),1)+ rho_u2.*Ux(fluxfaceList(:,2),1) )/2;
rho_u(:,2)=(  rho_u1.*Uy(fluxfaceList(:,1),1)+ rho_u2.*Uy(fluxfaceList(:,2),1) )/2;
    
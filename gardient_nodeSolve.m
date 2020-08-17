
function  [ dT ] =  gardient_nodeSolve(Nnode,fluxfaceList,fluxfaceList1,T)
[N,~]=size(fluxfaceList);
dT=zeros(Nnode,2);
for i=1:N
    dT( fluxfaceList(i,1) ,1)=dT( fluxfaceList(i,1) ,1)+( T(fluxfaceList(i,2))- T(fluxfaceList(i,1)))*fluxfaceList1(i,1)*(-fluxfaceList1(i,4))/2;
        dT( fluxfaceList(i,2) ,1)=dT( fluxfaceList(i,2) ,1)+( T(fluxfaceList(i,1))- T(fluxfaceList(i,2)))*fluxfaceList1(i,1)*(fluxfaceList1(i,4))/2;
        dT( fluxfaceList(i,1) ,2)=dT( fluxfaceList(i,1) ,2)+( T(fluxfaceList(i,2))- T(fluxfaceList(i,1)))*fluxfaceList1(i,1)*(-fluxfaceList1(i,5))/2;
        dT( fluxfaceList(i,2) ,2)=dT( fluxfaceList(i,2) ,2)+( T(fluxfaceList(i,1))- T(fluxfaceList(i,2)))*fluxfaceList1(i,1)*(fluxfaceList1(i,5))/2;
end
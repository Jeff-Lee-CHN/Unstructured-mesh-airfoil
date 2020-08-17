% function Conduction2Dmain
%% 1. Imposed the mesh data
clear % clear all parameters
clc % clear command window
tic
Filename='naca0012whk.msh';
[Nnode,~,nodeList,elemList]=readmsh(Filename);
elemList0=elemList;
elemList0(:,3:5)=[]; 
%% 2. Preproces
%barycenter of each element
[barycenter]=trianglecenter(elemList0,nodeList);
%edge and neighbour element of fluxface
[fluxfaceList]=assocoatedgrid1(elemList0);%第一二列为相邻单元，三四列为相邻节点
%求解伴生网格面信息
[fluxfaceList1]=fluxface(fluxfaceList,barycenter,nodeList);%第一列为长度，第二三列为对应面积，第四五列为法向量坐标
%节点信息列表包含体积信息
[Nfluxface,~]=size(fluxfaceList);
nodeList1=zeros(Nnode,1);%第一列为节点的体积
for i=1:Nfluxface
    nodeList1(fluxfaceList(i,3),1)=nodeList1(fluxfaceList(i,3),1)+fluxfaceList1(i,2);
    nodeList1(fluxfaceList(i,4),1)=nodeList1(fluxfaceList(i,4),1)+fluxfaceList1(i,3);
end
V=diag(nodeList1(:,1));
V_1=diag((1./ nodeList1(:,1)));
%edge的长度
Length_edge=sqrt( (nodeList(fluxfaceList(:,3),2)-nodeList(fluxfaceList(:,4),2)  ).^2+...
        (nodeList(fluxfaceList(:,3),3)-nodeList(fluxfaceList(:,4),3)  ).^2 );
% Vector tau
tau=zeros(Nfluxface,2);
tau(:,1)=( nodeList(fluxfaceList(:,4),2)-nodeList(fluxfaceList(:,3),2) )./ Length_edge;
tau(:,2)=( nodeList(fluxfaceList(:,4),3)-nodeList(fluxfaceList(:,3),3) )./Length_edge;
%% boundaryList
bouelem_out=find(  (elemList(:,2)==1)&(elemList(:,5)==2));
bouelem_int=find(  (elemList(:,2)==1)&(elemList(:,5)==1));
[BoundaryNode_int]=SolveBoundaryList(bouelem_int,elemList0,nodeList,barycenter);
[BoundaryNode_out]=SolveBoundaryList(bouelem_out,elemList0,nodeList,barycenter);
toc
%% 3. Iterative solve
%% 常数
 mu=1.79*1e-5;
R=287.1;
gamma=1.4;
Pr=0.72;
Ma=0.85;%Ma number
alpha=0.0;%迎角
Cv=R/(0.4);
Cp=Cv*gamma;
T_tot=288.15;
P_tot=101325;
%%
T_field=T_tot/(1+(gamma-1)/2*Ma^2);
P_field=P_tot/(1+(gamma-1)/2*Ma^2)^(gamma/(gamma-1));
rho_field=P_field/R/T_field;
Ux_field=Ma*sqrt(gamma*R*T_field)*cos(alpha);
Uy_field=Ma*sqrt(gamma*R*T_field)*sin(alpha);
%%
% Using the green function to solve the gradient of each node, and solving
% the gradient vector of fluxface by weighted average method.
% The time step
deltT=1e-6;
% Initization and boundary conditions
Ux=Ux_field*ones(Nnode,1);%velocity of x coordinate
Uy=Uy_field*ones(Nnode,1);%velocity of x coordinate
rho=rho_field*ones(Nnode,1);%density
T=T_field*ones(Nnode,1);
P=P_field*ones(Nnode,1);

%% 守恒量计算
Result_old=zeros(Nnode,4);%计算结果，第一列为连续性方程，第二列为动量方程ux，第三列为动量方程uy，第四列为能量方程
E=P./rho/0.4+(Ux.^2+Uy.^2)/2;
Result_old(:,1)=rho;
Result_old(:,2)=rho.*Ux;
Result_old(:,3)=rho.*Uy;
Result_old(:,4)=rho.*E;
%%
err_2=1;
num=0;
% Start solve
tic % start the time
% while(err_2>1e-4)
 for num=1:1000000
%
rhoH=Result_old(:,4)+P;
H=rhoH./rho;
mu=(1.45*T.^1.5)./(T+110)*(1e-6);
k=Cp*mu/Pr;
Ux(BoundaryNode_int(:,1))=0;
Uy(BoundaryNode_int(:,1))=0;
%
    [ dT ] =  gardient_nodeSolve(Nnode,fluxfaceList(:,3:4),fluxfaceList1,T);
    [ dux ] =  gardient_nodeSolve(Nnode,fluxfaceList(:,3:4),fluxfaceList1,Ux);
    [ duy ] =  gardient_nodeSolve(Nnode,fluxfaceList(:,3:4),fluxfaceList1,Uy);
    NablaT_node=V_1*dT;%NablaT in node
    NablaUx_node=V_1*dux;%NablaUx in node
    NablaUy_node=V_1*duy;%NablaUy in node
    %Correction of dissipation term in fluxface
    [NablaT_fluxface]= gardientOnFluxface(  NablaT_node,  fluxfaceList(:,3:4),T,Length_edge, tau ); % gradient on fulxface
    [NablaUx_fluxface]= gardientOnFluxface(  NablaUx_node,  fluxfaceList(:,3:4),Ux,Length_edge, tau );
    [NablaUy_fluxface]= gardientOnFluxface(  NablaUy_node,  fluxfaceList(:,3:4),Uy,Length_edge, tau );
    %fluxface上的对流项
    P_fluxface=( P(fluxfaceList(:,3),1)+P(fluxfaceList(:,4),1) )/2;
    Ux_fluxface=( Ux(fluxfaceList(:,3),1)+Ux(fluxfaceList(:,4),1) )/2;
    Uy_fluxface=( Uy(fluxfaceList(:,3),1)+Uy(fluxfaceList(:,4),1) )/2;
  %% Rho格式
  uL=Ux(fluxfaceList(:,3));      uR=Ux(fluxfaceList(:,4));
  vL=Uy(fluxfaceList(:,3));       vR=Uy(fluxfaceList(:,4));
  rhoL=rho(fluxfaceList(:,3));   rhoR=rho(fluxfaceList(:,4));
  HL=H(fluxfaceList(:,3));         HR=H(fluxfaceList(:,4));
  nx=-fluxfaceList1(:,4);          ny=-fluxfaceList1(:,5);  
  PL=P(fluxfaceList(:,3));          PR=P(fluxfaceList(:,4));
  %
  rho_bo=sqrt( rhoL.* rhoR );
  u_bo=( uL.*sqrt(rhoL)+uR.*sqrt(rhoR) )./( sqrt(rhoL)+sqrt(rhoR) );
  v_bo=( vL.*sqrt(rhoL)+vR.*sqrt(rhoR) )./( sqrt(rhoL)+sqrt(rhoR) );
  H_bo=( HL.*sqrt(rhoL)+HR.*sqrt(rhoR) )./( sqrt(rhoL)+sqrt(rhoR) );
 V_bo=u_bo.*nx+v_bo.*ny;
  q_bo2=u_bo.^2+v_bo.^2;
  c_bo=sqrt(  (gamma-1)*(H_bo-q_bo2/2) );
  %
  deltV=(uR.*nx+vR.*ny)-(uL.*nx+vL.*ny);
  deltP=PR-PL;
  deltu=uR-uL;
  deltv=vR-vL;
deltrho=rhoR-rhoL;
  deltF1=abs(V_bo-c_bo).*( ( deltP-rho_bo.*c_bo.*deltV)./( 2*(c_bo.^2)) ).*[ones(Nfluxface,1), u_bo-c_bo.*nx , v_bo-c_bo.*ny, H_bo-c_bo.*V_bo ];
  deltF234=abs(V_bo).*  (  ( deltrho -deltP./(c_bo.^2) ).*[ ones(Nfluxface,1), u_bo, v_bo, q_bo2/2  ] + rho_bo.*[ zeros(Nfluxface,1),deltu-deltV.*nx, deltv-deltV.*ny , u_bo.*deltu+v_bo.*deltv-V_bo.*deltV ]  );
  deltF5=abs( V_bo+c_bo ).*(  ( deltP+rho_bo.*c_bo.*deltV)./( 2*(c_bo.^2)) ).*[ones(Nfluxface,1),u_bo+c_bo.*nx, v_bo+c_bo.*ny, H_bo+c_bo.*V_bo];
  Rho1=(1/2)*(deltF1+deltF234+deltF5);
%% 残差
    Residual=zeros(Nnode,4);
    rofulx=(  Result_old(fluxfaceList(:,3),1)+ Result_old(fluxfaceList(:,4),1) )/2;
    roufulx=(  Result_old(fluxfaceList(:,3),2)+ Result_old(fluxfaceList(:,4),2) )/2;
    rovfulx=(  Result_old(fluxfaceList(:,3),3)+ Result_old(fluxfaceList(:,4),3) )/2; 
    roHfulx=(  Result_old(fluxfaceList(:,3),4)+P(fluxfaceList(:,3),1)+ Result_old(fluxfaceList(:,4),4)+P(fluxfaceList(:,4),1)  )/2; 
V=fluxfaceList1(:,4).*Ux_fluxface+fluxfaceList1(:,5).* Uy_fluxface;
%中心差分
% for i=1:Nfluxface
%         Residual( fluxfaceList(i,3) ,1 ) =Residual( fluxfaceList(i,3) ,1 )-(rofulx(i)*(-V(i)) )*fluxfaceList1(i,1);
%         Residual( fluxfaceList(i,4) ,1 ) =Residual( fluxfaceList(i,4) ,1 )-(rofulx(i)*V(i)  )*fluxfaceList1(i,1);
%         
%         Residual( fluxfaceList(i,3) ,2 ) =Residual( fluxfaceList(i,3) ,2 )-(roufulx(i)*(-V(i))+P_fluxface(i,1)*(-fluxfaceList1(i,4)))*fluxfaceList1(i,1);
%         Residual( fluxfaceList(i,4) ,2 ) =Residual( fluxfaceList(i,4) ,2 )-(roufulx(i)*V(i)+P_fluxface(i,1)*fluxfaceList1(i,4))*fluxfaceList1(i,1);
%         
%         Residual( fluxfaceList(i,3) ,3 ) =Residual( fluxfaceList(i,3) ,3 )-(rovfulx(i)*(-V(i))+P_fluxface(i,1)*(-fluxfaceList1(i,5)))*fluxfaceList1(i,1);
%         Residual( fluxfaceList(i,4) ,3 ) =Residual( fluxfaceList(i,4) ,3 )-(rovfulx(i)*V(i)+P_fluxface(i,1)*fluxfaceList1(i,5) )*fluxfaceList1(i,1);
%         
%         Residual( fluxfaceList(i,3) ,4 ) =Residual( fluxfaceList(i,3) ,4 )- ( roHfulx(i)*(-V(i))  )*fluxfaceList1(i,1);
%         Residual( fluxfaceList(i,4) ,4 ) =Residual( fluxfaceList(i,4) ,4 )- (  roHfulx(i)*V(i))*fluxfaceList1(i,1);
% end
%%%%%%
 [ rho_V ]=convectionOnFluxface( fluxfaceList(:,3:4),Ux,Uy, Result_old(:,1)  );
[ rhou_V ]=convectionOnFluxface( fluxfaceList(:,3:4),Ux,Uy, Result_old(:,2)  );
[ rhov_V ]=convectionOnFluxface( fluxfaceList(:,3:4),Ux,Uy, Result_old(:,3)  );
[ rhoH_V ]=convectionOnFluxface( fluxfaceList(:,3:4),Ux,Uy, rhoH(:,1)  );
for i=1:Nfluxface
        Residual( fluxfaceList(i,3) ,1 ) =Residual( fluxfaceList(i,3) ,1 )-(rho_V(i,1)*nx(i)+rho_V(i,2)*ny(i) )*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,1 ) =Residual( fluxfaceList(i,4) ,1 )-(rho_V(i,1)*(-nx(i))+rho_V(i,2)*(-ny(i))  )*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,2 ) =Residual( fluxfaceList(i,3) ,2 )-(rhou_V(i,1)*nx(i)+rhou_V(i,2)*ny(i)+P_fluxface(i,1)*nx(i))*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,2 ) =Residual( fluxfaceList(i,4) ,2 )-(rhou_V(i,1)*(-nx(i))+rhou_V(i,2)*(-ny(i))+P_fluxface(i,1)*(-nx(i)))*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,3 ) =Residual( fluxfaceList(i,3) ,3 )-(rhov_V(i,1)*nx(i)+rhov_V(i,2)*ny(i)+P_fluxface(i,1)*ny(i))*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,3 ) =Residual( fluxfaceList(i,4) ,3 )-(rhov_V(i,1)*(-nx(i))+rhov_V(i,2)*(-ny(i))+P_fluxface(i,1)*(-ny(i)) )*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,4 ) =Residual( fluxfaceList(i,3) ,4 )- ( rhoH_V(i,1)*nx(i)+rhoH_V(i,2)*ny(i) )*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,4 ) =Residual( fluxfaceList(i,4) ,4 )- ( rhoH_V(i,1)*(-nx(i))+rhoH_V(i,2)*(-ny(i)) )*fluxfaceList1(i,1);
end
%roe迎风修正
for i=1:Nfluxface
        Residual( fluxfaceList(i,3) ,1 ) =Residual( fluxfaceList(i,3) ,1 )+Rho1(i,1)*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,1 ) =Residual( fluxfaceList(i,4) ,1 )-Rho1(i,1)*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,2 ) =Residual( fluxfaceList(i,3) ,2 )+ Rho1(i,2) *fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,2 ) =Residual( fluxfaceList(i,4) ,2 )-Rho1(i,2)*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,3 ) =Residual( fluxfaceList(i,3) ,3 )+Rho1(i,3)*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,3 ) =Residual( fluxfaceList(i,4) ,3 )-Rho1(i,3)*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,4 ) =Residual( fluxfaceList(i,3) ,4 )+Rho1(i,4) *fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,4 ) =Residual( fluxfaceList(i,4) ,4 )- Rho1(i,4)*fluxfaceList1(i,1);
end
%粘性项
mu_fluxface=(mu(fluxfaceList(i,3))+mu(fluxfaceList(i,4)))/2;
k_fluxface=(k(fluxfaceList(i,3))+k(fluxfaceList(i,4)))/2;
    div_V=NablaUx_fluxface( :,1)+NablaUy_fluxface(:,2);
    tauxx=2*mu_fluxface.*( NablaUx_fluxface( :,1) -(1/3)*div_V);
    tauyy=2*mu_fluxface.*( NablaUy_fluxface( :,2) -(1/3)*div_V);
    tauxy=mu_fluxface.*(  NablaUx_fluxface(:,2)+NablaUy_fluxface(:,1) );
    tauyx=tauxy;
    thetax= Ux_fluxface(:,1).*tauxx+ Uy_fluxface(:,1).*tauxy+k_fluxface.*NablaT_fluxface(:,1);
    thetay= Ux_fluxface(:,1).*tauyx+ Uy_fluxface(:,1).*tauyy+k_fluxface.*NablaT_fluxface(:,2);
    for i=1:Nfluxface
        Residual( fluxfaceList(i,3) ,2 ) =Residual( fluxfaceList(i,3) ,2 )+ (   nx(i)*tauxx(i,1)+ny(i)*tauxy(i,1)    )*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,2 ) =Residual( fluxfaceList(i,4) ,2 )+(    (-nx(i))*tauxx(i,1)+(-ny(i))*tauxy(i,1)     )*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,3 ) =Residual( fluxfaceList(i,3) ,3 )+ (    ny(i)*tauyy(i,1) +nx(i)*tauyx(i,1)     )*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,3 ) =Residual( fluxfaceList(i,4) ,3 )+ (     (-ny(i) )*tauyy(i,1)+(-nx(i))*tauyx(i,1)     )*fluxfaceList1(i,1);
        
        Residual( fluxfaceList(i,3) ,4 ) =Residual( fluxfaceList(i,3) ,4 )+  (nx(i) *thetax(i,1)+ ny(i)*thetay(i,1) )*fluxfaceList1(i,1);
        Residual( fluxfaceList(i,4) ,4 ) =Residual( fluxfaceList(i,4) ,4 )+  ( -nx(i)*thetax(i,1)-ny(i)*thetay(i,1) )*fluxfaceList1(i,1);
    end
 %% 边界条件(无滑移边界条件)
%内边界条件
%         Residual(BoundaryNode_int(:,1),2)=Residual(BoundaryNode_int(:,1),2)-P(BoundaryNode_int(:,1),1).*BoundaryNode_int(:,2);
%         Residual(BoundaryNode_int(:,1),3)=Residual(BoundaryNode_int(:,1),3)-P(BoundaryNode_int(:,1),1).*BoundaryNode_int(:,3);
%         
%                 Residual(BoundaryNode_int(:,1),2)=Residual(BoundaryNode_int(:,1),2)+(2*mu*(NablaUx_node(BoundaryNode_int(:,1),1)-1/3*(NablaUx_node(BoundaryNode_int(:,1),1) +...
%             NablaUy_node(BoundaryNode_int(:,1),2)))).*BoundaryNode_int(:,2)+(mu*(NablaUx_node(BoundaryNode_int(:,1),2)+NablaUy_node(BoundaryNode_int(:,1),1))).*BoundaryNode_int(:,3);
%     
%         Residual(BoundaryNode_int(:,1),3)=Residual(BoundaryNode_int(:,1),3)+(2*mu*(NablaUy_node(BoundaryNode_int(:,1),2)-1/3*(NablaUx_node(BoundaryNode_int(:,1),1) +...
%             NablaUy_node(BoundaryNode_int(:,1),2)))).*BoundaryNode_int(:,3)+(mu*(NablaUx_node(BoundaryNode_int(:,1),2)+NablaUy_node(BoundaryNode_int(:,1),1))).*BoundaryNode_int(:,2);
 Residual(BoundaryNode_int(:,1),2)=0;
Residual(BoundaryNode_int(:,1),3)=0; 

           Residual(BoundaryNode_int(:,1),4)=Residual(BoundaryNode_int(:,1),4)+BoundaryNode_int(:,2).*(k(BoundaryNode_int(:,1)).*NablaT_node(BoundaryNode_int(:,1),1))+BoundaryNode_int(:,3).*(k(BoundaryNode_int(:,1)).*NablaT_node(BoundaryNode_int(:,1),2));
%外边界条件
    Residual( BoundaryNode_out(:,1),1:4 )=0; 
%% 结果更新
    Residual=deltT.*(V_1*Residual);
    Result_new=Residual+Result_old;
    
     Result_new(BoundaryNode_int(:,1),2)=0;
Result_new(BoundaryNode_int(:,1),3)=0;
    %
   rho= Result_new(:,1);
    %    
    Ux=Result_new(:,2)./rho;
    Uy=Result_new(:,3)./rho;
    P=0.4*( Result_new(:,4)-( Result_new(:,2).^2+ Result_new(:,3).^2)./rho/2);
    T=P./rho/R;
    %
    err_1=max(abs(Result_new(:,1)-Result_old(:,1)));
    err_2=max(abs(Result_new(:,2)-Result_old(:,2)));
    err_3=max(abs(Result_new(:,3)-Result_old(:,3)));
    err_4=max(abs(Result_new(:,4)-Result_old(:,4)));
    % The mxiumum error between adjacent time steps and update tempearture
%     num=num+1;
    % Print the iter steps and the maximum error
    if(mod(num,1000)==0)
        disp(['step: ',num2str(num),' Rhoerror: ',num2str( err_1), ' Uxerror: ',num2str( err_2) ,' Uyerror:  ',num2str( err_3), ' Perror: ',num2str( err_4),' Terror: ',num2str( err_4)]);
    end
    Result_old=Result_new;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

end % end solve
toc % end compute
 %% output result
resultoutput(Filename,Nnode,T,rho,P,Ux,Uy,gamma,R);
 %%
x=nodeList(:,2);
y=nodeList(:,3);
z=T;
p=P;
ux=Ux;
uy=Uy;
figure(1)
plot3(x,y,z,'bo');
title('温度')
figure(2)
plot3(x,y,p,'bo'); 
title('压力')
figure(3)
plot3(x,y,ux,'bo'); 
title('ux')
figure(4)
plot3(x,y,uy,'bo'); 
title('uy')
figure(5)
plot3(x,y,rho,'bo'); 
title('密度')
% %%%%%%%%%%%%%
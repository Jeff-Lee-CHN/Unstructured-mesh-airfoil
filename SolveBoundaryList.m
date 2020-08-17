
function [BoundaryNode]=SolveBoundaryList(bouelem,elemList,nodeList,barycenter)
[Nnode,~]=size(nodeList);
Nbou=length(bouelem);
BoundaryElem=zeros(Nbou,5);
BoundaryElem(:,1:2)=elemList(bouelem,3:4);%�߽絥Ԫ��������
numberBou=elemList(bouelem,1);%�߽絥Ԫ���
BoundaryNode=unique(BoundaryElem(:,1:2));%�߽�ڵ����
%������
intelem=find(elemList(:,2)==2);
IntElem=elemList(intelem,:);
[Nint,~]=size(IntElem);
BoundaryElem(:,5)=sqrt( ( nodeList(BoundaryElem(:,1),2 )- nodeList(BoundaryElem(:,2),2 ) ).^2+...
                                         (nodeList(BoundaryElem(:,1),3)- nodeList(BoundaryElem(:,2),3)).^2  );
x=nodeList(BoundaryElem(:,1),2)-nodeList(BoundaryElem(:,2),2);
y=nodeList(BoundaryElem(:,1),3)-nodeList(BoundaryElem(:,2),3);
for i=1:Nbou
    for j=1:Nint
%         sameElement=intersect(BoundaryElem(i,:),IntElem(j,3:5));
%         Nsame=length(sameElement);
        %%%%%%%%
        Nsame=0;
        for k=1:3
            for l=1:3
                if BoundaryElem(i,k)==IntElem(j,l+2)
                    Nsame=Nsame+1;
                end
            end
        end
        %%%%%%%%
        if Nsame==2
 vector(1,1)=barycenter(numberBou(i),1 )-barycenter(IntElem(j,1),1);
 vector(1,2)=barycenter(numberBou(i),2 )-barycenter(IntElem(j,1),2);            
        [n]=normalVector(x(i),y(i),vector);
        BoundaryElem(i,3:4)=n;
        end
    end
end
%ѭ���߽���Ϣ������
BoundaryList=zeros(Nnode,2);
for i=1:Nbou
    BoundaryList(BoundaryElem(i,1),1)= BoundaryList(BoundaryElem(i,1),1)+BoundaryElem(i,3)*BoundaryElem(i,5)/2;
    BoundaryList(BoundaryElem(i,2),1)= BoundaryList(BoundaryElem(i,2),1)+BoundaryElem(i,3)*BoundaryElem(i,5)/2;    
    %
    BoundaryList(BoundaryElem(i,1),2)= BoundaryList(BoundaryElem(i,1),2)+BoundaryElem(i,4)*BoundaryElem(i,5)/2;
    BoundaryList(BoundaryElem(i,2),2)= BoundaryList(BoundaryElem(i,2),2)+BoundaryElem(i,4)*BoundaryElem(i,5)/2;    
end
BoundaryNode(:,2:3)=BoundaryList( BoundaryNode(:,1),: );%���ձ߽緽��
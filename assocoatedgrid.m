%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ܣ����ְ�������
%���룺��Ԫ��Ϣ
%���أ�����������Ϣ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [list]=assocoatedgrid(elem)
%%%%%%%%%%%%%%%%%%
[Nelem,~]=size(elem);
%%%%%%%%%%%%%%%%%%
list=zeros(2*Nelem,4);
num=0;
for i=1:Nelem-1
    for j=i+1:Nelem
        %%%%%%%%%%%
        ifbou=0;
        if elem(j,5)==0
        else
           sameElement=intersect(elem(i,3:5),elem(j,3:5));%��intersect���ø�ʺһ��
            ifbou=length( sameElement);
        end
        %%%%%%%%%%%
        if ifbou==2
            num=num+1;
            list(num,1)=i;%���ڵ�Ԫ�����Ϊelem�е�˳�򣬷�������
            list(num,2)=j;%���ڵ�Ԫ�����Ϊelem�е�˳�򣬷�������
            list(num,3:4)=sameElement;%���ڽڵ㣬���Ϊelem�е�˳�򣬷�������
        end
        %%%%%%%%%%%
    end
end
list(num+1:end,:)=[];








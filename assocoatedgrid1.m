%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ܣ����ְ�������
%���룺��Ԫ��Ϣ
%���أ�����������Ϣ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [list1]=assocoatedgrid1(elem)
%%%%%%%%%%%%%%%%%%
[Nelem,~]=size(elem);
%%%%%%%%%%%%%%%%%%
list1=zeros(2*Nelem,4);
num=0;
for i=1:Nelem-1
    for j=i+1:Nelem
        ifbou1=0;
        ifbou=ones(1,3);%�ж��Ƿ����
        %%%%%%%%%%%%
        for k=1:3
            for l=1:3
                if elem(j,5)==0
                else
                    if elem(i,k+2)==elem(j,l+2)
                        ifbou(k)=0;
                        ifbou1=ifbou1+1;
                    end
                end
            end
        end
        %%%%%%%%%%%
        if ifbou1==2
            num=num+1;
            list1(num,1)=i;%���ڵ�Ԫ�����Ϊelem�е�˳�򣬷�������
            list1(num,2)=j;%���ڵ�Ԫ�����Ϊelem�е�˳�򣬷�������
             m=0;
            for k=1:3
                if ifbou(k)==0
                    list1(num,3+m)=elem(i,k+2);%���ڽڵ���
                    m=m+1;
                end
            end
        end
        %%%%%%%%%%%
    end
end
list1(num+1:end,:)=[];
%%%%%%%%%%%%%%%%%%%%








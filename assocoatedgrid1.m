%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：划分伴生网格
%输入：单元信息
%返回：伴生网格信息
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
        ifbou=ones(1,3);%判断是否相等
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
            list1(num,1)=i;%相邻单元，序号为elem中的顺序，方便索引
            list1(num,2)=j;%相邻单元，序号为elem中的顺序，方便索引
             m=0;
            for k=1:3
                if ifbou(k)==0
                    list1(num,3+m)=elem(i,k+2);%相邻节点编号
                    m=m+1;
                end
            end
        end
        %%%%%%%%%%%
    end
end
list1(num+1:end,:)=[];
%%%%%%%%%%%%%%%%%%%%








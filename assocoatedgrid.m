%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：划分伴生网格
%输入：单元信息
%返回：伴生网格信息
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
           sameElement=intersect(elem(i,3:5),elem(j,3:5));%此intersect慢得跟屎一样
            ifbou=length( sameElement);
        end
        %%%%%%%%%%%
        if ifbou==2
            num=num+1;
            list(num,1)=i;%相邻单元，序号为elem中的顺序，方便索引
            list(num,2)=j;%相邻单元，序号为elem中的顺序，方便索引
            list(num,3:4)=sameElement;%相邻节点，序号为elem中的顺序，方便索引
        end
        %%%%%%%%%%%
    end
end
list(num+1:end,:)=[];








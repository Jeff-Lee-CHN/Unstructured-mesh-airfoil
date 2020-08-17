%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：求每个单元重心
%输入：单元列表，节点列表
%返回：外心xye坐标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [barycenter]=trianglecenter(elem,nodeList)
[Node,~]=size(elem);
%%%%%%%%%%%%%%%%%
barycenter=zeros(Node,2);%重心
for i=1:Node
    if elem(i,2)==2
        barycenter(i,1)=(nodeList(  elem(i,3) , 2 )+nodeList(  elem(i,4) , 2 )+nodeList(  elem(i,5) , 2 )   )/3;
        barycenter(i,2)=( nodeList(  elem(i,3) , 3 )+nodeList(  elem(i,4) , 3)+ nodeList(  elem(i,5) , 3 )  )/3;
    else
        if elem(i,2)==1
            barycenter(i,1)=(nodeList(  elem(i,3) , 2 )+nodeList(  elem(i,4) , 2 ) )/2;
            barycenter(i,2)=( nodeList(  elem(i,3) , 3 )+nodeList(  elem(i,4) , 3 ) )/2;
        else
            barycenter(i,1)=nodeList(  elem(i,3) , 2 );
            barycenter(i,2)=nodeList(  elem(i,3) , 3 );
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能：导入网格文件，提取信息
%输入：网格文件名称
%返回：节点数及其列表，单元数及其列表
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Nnode,Nelem,node,elem]=readmsh(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat=textread(filename,'%s','delimiter','\n');%read file as string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
ifbreak=1;
while ifbreak==1
    dd=cell2mat(deblank(dat(i)));
    length_txt=length(dd);
    if length_txt==6
        if dd=='$Nodes'
            Nnode=str2double(dat(i+1)) ;
           break
        end
    end
    i=i+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node=zeros(Nnode,4);
for j=1:Nnode
    line_num=i+1+j;
    data_1=cell2mat(deblank(dat(line_num)))  ;
   node(j,:)=str2num(data_1);
end
i=i+Nnode+4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nelem=str2double(dat(i)) ; 
elem=zeros(Nelem,9);
for j=1:Nelem
    line_num=j+i;
data_1=cell2mat(deblank(dat(line_num)));
data_2=str2num(data_1);
[~,length_elem]=size(data_2);
elem(j,1:length_elem)=data_2;
end



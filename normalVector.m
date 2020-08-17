%����������������x��y��ֵ���Լ�����vector�����жϷ������������vector�ĵ��Ϊ��
function [n]=normalVector(x,y,vector)
N=length(x);
n=zeros(N,2);
for i=1:N
    if abs(x(i,1))<1E-9
        n(i,2)=0;
        if vector(i ,1)<0
            n(i,1)= -1;
        else
            n(i,1)=1;
        end
    else
        n(i,2)=sqrt(1/ ( 1+ y(i)^2/x(i)^2 ));
        n(i,1)=-(n(i,2)*y(i)/x(i));
        if (vector(i ,1)* n(i,1)+vector(i ,2)* n(i,2))<0
              n(i,1)= -n(i,1);
            n(i,2)= -n(i,2);
        end
    end
end
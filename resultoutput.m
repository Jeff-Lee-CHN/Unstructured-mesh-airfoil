function resultoutput(filename,Nnode,T,rho,P,Ux,Uy,gamma,R)
 %% Ma number
 Vv=sqrt( Ux.^2+Uy.^2 );
Maa=Vv./sqrt(gamma.*R.*T);
%%
tmp=0.0;
 delete('result1.msh');
 FileMesh = fopen(filename,'r') ;
FileResult=fopen('result1.msh','a+');

while ~feof(FileMesh)
    str = fgetl(FileMesh) ;
    fprintf(FileResult,'%s\n',str);
end
%
fprintf(FileResult,'$NodeData\n1\n\"T\"\n1\n0.0\n3\n0\n1\n%s\n',num2str(Nnode));
for i=1:Nnode
    fprintf(FileResult,'%s %s\n',num2str(i),num2str(T(i)));
end
fprintf(FileResult,'%s\n','$EndNodeData');
%
fprintf(FileResult,'$NodeData\n1\n\"rho\"\n1\n0.0\n3\n0\n1\n%s\n',num2str(Nnode));
for i=1:Nnode
    fprintf(FileResult,'%s %s\n',num2str(i),num2str(rho(i)));
end
fprintf(FileResult,'%s\n','$EndNodeData');
%
fprintf(FileResult,'$NodeData\n1\n\"p\"\n1\n0.0\n3\n0\n1\n%s\n',num2str(Nnode));
for i=1:Nnode
    fprintf(FileResult,'%s %s\n',num2str(i),num2str(P(i)));
end
fprintf(FileResult,'%s\n','$EndNodeData');
%
fprintf(FileResult,'$NodeData\n1\n\"vel\"\n1\n0.0\n3\n0\n3\n%s\n',num2str(Nnode));
for i=1:Nnode
    fprintf(FileResult,'%s %s %s %s\n',num2str(i),num2str(Ux(i)),num2str(Uy(i)),num2str(tmp));
end
fprintf(FileResult,'%s\n','$EndNodeData');
%
%Ma
fprintf(FileResult,'$NodeData\n1\n\"Ma\"\n1\n0.0\n3\n0\n1\n%s\n',num2str(Nnode));
for i=1:Nnode
    fprintf(FileResult,'%s %s\n',num2str(i),num2str(Maa(i)));
end
fprintf(FileResult,'%s\n','$EndNodeData');
%
 fclose(FileResult);
function [U,V]=GNMFLMI(inputdata,params)

%  Wang, M.N. et al. GNMFLMI: Graph Regularized Nonnegative Matrix Factorization   
%  for Predicting LncRNA-MiRNA Interactions.
%  2019/09/25


p=params.p;
LM_mat=inputdata.LM_mat;
SL_mat=inputdata.SL_mat;
SM_mat=inputdata.SM_mat;


Y=LM_mat; 
L_graph_mat = Graph( SL_mat , p );
M_graph_mat = Graph( SM_mat , p ); 
S_M = SM_mat.* M_graph_mat ;  
S_L = SL_mat.* L_graph_mat ;
clear p;
 
k=params.k;
iterate=params.iterate; 
beta=params.beta;   
lamda_l=params.lamda_l;
lamda_m=params.lamda_m; 
beta > 0 && lamda_l >0 && lamda_m >0;
fprintf('k=%d  maxiter=%d  beta=%d  lamda_l=%d lamda_m=%d\n', k, iterate,beta, lamda_l, lamda_m); 

[rows,cols] = size(Y);
U=abs(rand(k,rows));        
V=abs(rand(k,cols));

D_L = diag(sum(S_L,2));
D_M = diag(sum(S_M,2));
L_L=D_L-S_L;
L_M=D_M-S_M;


fid = fopen( 'RunResult.txt','wt+');
for step=1:iterate
        U1=U.*((V*Y'+lmada_l*U*S_L)./(V*V'*U+lmada_l*U*D_L+beta*U));
        V1=V.*((U1*Y+lmada_m*V*S_M)./(U1*U1'*V+lmada_m*V*D_M+beta*V));
         
        ULU = sum(diag((U1*L_L)*U1'));
        VLV = sum(diag((V1*L_M)*V1'));
        obj = sum(sum((Y-U1'*V1).^2))+beta*(sum(sum(U1.^2)) )+beta*(sum(sum(V1.^2)))+lmada_l*ULU+lmada_m*VLV; 
        
        error=max([max(sum((U1-U).^2)),max(sum((V1-V).^2))]);      
        
        fprintf(fid,'%s\n',[sprintf('step = \t'),int2str(step),...
            sprintf('\t obj = \t'),num2str(obj),...
		    sprintf('\t error = \t'), num2str(error)]);
        fprintf('step=%d  obj=%d  error=%d\n',step, obj, error);   
        if error< 10^(-4)
            fprintf('step=%d\n',step);
            break;
        end
        
        U=U1; 
        V=V1;
        
end
fclose(fid);

end

 

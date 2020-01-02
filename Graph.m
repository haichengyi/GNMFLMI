function [ graph ] = Graph( network , p )

tic    
[rows, cols] = size( network );   
        
    %------------------p nearest neighbors--------------%
    
    network= network-diag(diag(network)); 
    PNN = zeros(rows, cols);
    graph = zeros(rows, cols);
    [sort_network,idx]=sort(network,2,'descend');
    for i = 1 : rows
        PNN(i,idx(i,1:p))=sort_network(i,1:p);
    end    
    for i = 1 : rows
        idx_i=find(PNN(i,:));
        for j = 1 : rows            
            idx_j=find(PNN(j,:));
            if ismember(j,idx_i) && ismember(i,idx_j)                
                graph(i,j)=1;
            elseif ~ismember(j,idx_i) && ~ismember(i,idx_j)   
                graph(i,j)=0;
            else
                graph(i,j)=0.5;               
            end       
       end
    end
    
 toc

end




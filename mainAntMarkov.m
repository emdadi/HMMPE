function [A, E,itert]=mainAntMarkov(seqCell,stateNum,line_count,characters)
% In God we trust
% AntMarkov
% Date: 27 shahrivar 96:
% time: 7:26

%% Initialization
%type(seqCell)
%characters=unique(seqCell);

maxiter=100;      % Maximum Number of Iterations
charNum=numel(characters);
tau0=1;	% Initial Phromone
pheromone_edge_augment=randi([100,1000]);
pheromone_node_augment=rand(1);
pheromone_edge_evaporate=pheromone_edge_augment/10;
tol=1e-6;
trtol=tol;
etol=tol;

A=zeros(stateNum,stateNum);
E=zeros(stateNum,charNum);

visited_edges=zeros(stateNum,stateNum);
visited_nodes=zeros(stateNum,charNum);

e=double(1)/stateNum;
StartMatrix=e*ones(stateNum,1);


e=double(1)/stateNum;
Startorigin=e*ones(stateNum,1);

Pheromone_edge=tau0*ones(stateNum,stateNum);   % Phromone Matrix
Pheromone_node=tau0*ones(stateNum,charNum); 
start_pheromone_edge=ones(stateNum,1);


%% Initial Loop

for i=1:line_count
        %% moveAnt Initial
        Symbols=seqCell{i};
        curr_node=0;
        states=zeros(numel(Symbols),1);
        strlen=numel(Symbols);
        for k =1:strlen
                      symbol=int32(Symbols(k));
                      flag=false;
                      while ~flag
                        j=randi([1, stateNum]);
                        j=int32(j);
                       
                            flag=true;
                            %type(j)
                            %type(symbol)
                            visited_nodes(j,symbol)=visited_nodes(j,symbol)+1;
                            if k>1
                                visited_edges(curr_node,j)=visited_edges(curr_node,j)+1;
                            end
                        
                      
                            curr_node=j;
                            states(k)=curr_node;
                        
                      end
        end
end

%% pheromone update
Pheromone_edge=Pheromone_edge+(pheromone_edge_augment).*(visited_edges);
Pheromone_node=Pheromone_node+((pheromone_node_augment).*(visited_nodes));
                                                                                                
%% Evaporate
stateNum=size(Pheromone_node,1); 
Pheromone_edge=(1-pheromone_edge_evaporate).*Pheromone_edge;
A=A+Pheromone_edge;
Pheromone_node=(1-pheromone_edge_evaporate)*Pheromone_node;
E=E+Pheromone_node;
totalA=sum(A,2);
totalE=sum(E,2);
E=E./(repmat(totalE,1,charNum));
A=A./(repmat(totalA,1,stateNum));
        
%% liklihood1
likly=0;
for k=1:line_count
   S=seqCell{k};
   %% bestWay
   logS=log(StartMatrix);
   [best_likly, best_node]=max(logS);
   total_likly=best_likly+log(E(best_node,int32(S(1))));
   for i=2:numel(S)
        pre_best_node=best_node;
        logArow=log(A(pre_best_node,:));
        [best_likly, best_node]=max(logArow);
        total_likly=total_likly+best_likly+log(E(best_node,S(i)));
   end
   %% 
    
   likly=likly+total_likly;
end
loglik=likly/line_count;

%% ACO Main Loop
converged=0;
itert=0;
while and(converged==0 , itert<=(maxiter))
        itert=itert+1;
        oldLL=loglik;
        oldguessTR=A;
        oldguessE= E;
        
        %% MoveAnt
        
        for l=1:line_count
            Symbols=seqCell{l};
            curr_node=1;
            for k=1:numel(Symbols)
                if and(Symbols(k)~=' ' , Symbols(k)~='\n')
                    symbol=int32(Symbols(k));
                              
                 
                    %% ProbCalc
                    probability=zeros(1,stateNum);
                    if k>0
                            
                            probability(:)=(log(A(curr_node,:))+Pheromone_edge(curr_node,:)).*transpose((log(E(:,symbol))+Pheromone_node(:,symbol)));
                    elseif k==0
                            probability(:)=(log(StartMatrix(:))+start_pheromone_edge(:)).*transpose((log(E(:,symbol))+Pheromone_node(:,symbol)));
                    end
               
                    sumProb=sum(probability,2);
                    probability=probability./(repmat(sumProb,1,stateNum));
                    %% SelectNextNode
                    
                    [bestprob,beststate]=max(probability);
                    j=beststate;
                    visited_nodes(j,symbol)=visited_nodes(j,symbol)+1;
                    if k==0
                        start_visited_edges(j)=start_visited_edges(j)+1;
                    else  
                        visited_edges(curr_node,j)=visited_edges(curr_node,j)+1;
                    end
                    curr_node=j;
                end
            end
        end
        
        %% pheromone update
        Pheromone_edge=Pheromone_edge+((pheromone_edge_augment).*(visited_edges));
        Pheromone_node=Pheromone_node+((pheromone_node_augment).*(visited_nodes));
        
        %% Evaporate
        stateNum=size(Pheromone_node,1); 
        Pheromone_edge=(1-pheromone_edge_evaporate).*Pheromone_edge;
        A=A+Pheromone_edge;
        Pheromone_node=(1-pheromone_edge_evaporate)*Pheromone_node;
        E=E+Pheromone_node;
        totalA=sum(A,2);
        totalE=sum(E,2);
        E=E./(repmat(totalE,1,charNum));
        A=A./(repmat(totalA,1,stateNum));
        
        %% liklihood1
        likly=0;
        for k=1:line_count
            seq=seqCell{k};
            %% bestWay
            logS=log(StartMatrix);
            [best_likly, best_node]=max(logS);
            total_likly=best_likly+log(E(best_node,S(1)));
            for i=2:numel(S)
                pre_best_node=best_node;
                logArow=log(A(pre_best_node,:));
                [best_likly, best_node]=max(logArow);
                total_likly=total_likly+best_likly+log(E(best_node,S(i)));
            end
            %% 
    
              likly=likly+total_likly;
         end
loglik=likly/line_count;

        normt=(norm(A-oldguessTR,'fro')/stateNum);
        norme=(norm(E-oldguessE,'fro')/stateNum);

        if (abs(loglik-oldLL)/(1+abs(oldLL)))< tol
                    if and(normt<trtol , norme<etol)
                            converged=1;
                            break 
                    end
        end
end
size_E=size(E)
size_A=size(A)
end

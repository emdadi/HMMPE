function [A, E,itert]=mainAntMarkov(seqCell,stateNum,line_count,characters)
%mainAntMarkov 
%Calculates estimates for a hidden Markov model parameter
%by antMarkov algorithm
%
% [A, E,itert]=mainAntMarkov(seqCell,stateNum,line_count,characters);

% Inputs:      seqCell    original emission matrix 
%              stateNum   each original transition matrix 
%              line_count     the set of training sequences
%              characters      the set of test sequences
% Outputs: 
%           A        estimator for observation probability matrix
%           A    estimated transition matrix calculated by AntMarkov algorithm
%           E     estimated emission matrices calculated by AntMarkov algorithm
%           itert   the number of iterations
% Please note: the set of sequnces should be given to this code as a cell
%

maxiter=100;      % Maximum Number of Iterations
charNum=numel(characters);     % The number of uniqe charecters in sequences 
tau0=1;  	% Initial Pheromone
pheromone_edge_augment=randi([100,1000]);       % The amount of pheromone that putting on edges 
pheromone_node_augment=rand(1);                  % The amount of pheromone that putting on nodes
pheromone_edge_evaporate=pheromone_edge_augment/10;
tol=1e-6;      % Thereshold of convergence
trtol=tol;     % Thereshold of convergence for Transition matrix   
etol=tol;        % Thereshold of convergence for Emission matrix

A=zeros(stateNum,stateNum);        % Transition Matrix    stateNum:the number of states of Hidden markov model
E=zeros(stateNum,charNum);         % Emission Matrix

visited_edges=zeros(stateNum,stateNum);
visited_nodes=zeros(stateNum,charNum);

e=double(1)/stateNum;
StartMatrix=e*ones(stateNum,1);        % The matrix for first step of movement


e=double(1)/stateNum;
Startorigin=e*ones(stateNum,1);

Pheromone_edge=tau0*ones(stateNum,stateNum);   % Pheromone-edge Matrix
Pheromone_node=tau0*ones(stateNum,charNum);     % Pheromone-node Matrix
start_pheromone_edge=ones(stateNum,1);


%% Initial Loop

for i=1:line_count           % line_count is the number of sequences 
        %% moveAnt Initial
        Symbols=seqCell{i};
        curr_node=0;         %current_node
        states=zeros(numel(Symbols),1);
        strlen=numel(Symbols);       % length of sequence
        for k =1:strlen
                      symbol=int32(Symbols(k));
                      flag=false;
                      while ~flag
                        j=randi([1, stateNum]);
                        j=int32(j);
                       
                            flag=true;
                           
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
                                                                                                
%% pheromone Evaporate
stateNum=size(Pheromone_node,1); 
Pheromone_edge=(1-pheromone_edge_evaporate).*Pheromone_edge;
A=A+Pheromone_edge;
Pheromone_node=(1-pheromone_edge_evaporate)*Pheromone_node;
E=E+Pheromone_node;
totalA=sum(A,2);
totalE=sum(E,2);
E=E./(repmat(totalE,1,charNum));
A=A./(repmat(totalA,1,stateNum));
        
%% likelihood
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
  
    
   likly=likly+total_likly;
end
loglik=likly/line_count;

%% ACO Main Loop
converged=0;
itert=0;       % the number of iteration
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
                              
                 
                    %% Probability Calculate
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
        
        %% pheromone Evaporate
        stateNum=size(Pheromone_node,1); 
        Pheromone_edge=(1-pheromone_edge_evaporate).*Pheromone_edge;
        A=A+Pheromone_edge;
        Pheromone_node=(1-pheromone_edge_evaporate)*Pheromone_node;
        E=E+Pheromone_node;
        totalA=sum(A,2);
        totalE=sum(E,2);
        E=E./(repmat(totalE,1,charNum));
        A=A./(repmat(totalA,1,stateNum));
        
        %% liklihood
        likly=0;
        for k=1:line_count
            seq=seqCell{k};
            %%  selecting bestWay
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

        normt=(norm(A-oldguessTR,'fro')/stateNum);       % Frobenius norm for Transition matrix
        norme=(norm(E-oldguessE,'fro')/stateNum);         % Frobenius norm for Emission matrix

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

function [A_c, E_c]=tabu(seqCell,stateNum,line_count,characters)
% be name khoda
% Date: 7 mehr 96:
% time: 12:16
tabusize=20;
testNo=20;
maxiter=100;
charNum=numel(characters);
insertion=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=rand([stateNum,stateNum,testNo]);
E=rand([stateNum,charNum,testNo]);    
e=double(1)/stateNum;
StartMatrix=e*ones(stateNum,1);

%% Normalization
for k=1:testNo
    totala=sum(A(:,:,k),2);
    for row=1:stateNum
        for column=1:stateNum
            A(row,column,k)=A(row,column,k)/totala(row);
        end
    end
    totale=sum(E(:,:,k),2);
    for row=1:stateNum
        for column=1:charNum
            E(row,column,k)=E(row,column,k)/totale(row);
        end
    end
end

%% run tabu

swaps=zeros(testNo ,5);

tabu_list=zeros(tabusize ,5);  
liklyArray=zeros(testNo,1);

for m=1:testNo
    
    %% liklihood
    likly=0;
    for t=1:line_count
        S=seqCell{t};
        %% probGivenModel
        len=numel(S);
        p=zeros(len,stateNum);
        
        p(1,:)=StartMatrix(:).*(E(:,S(1),m));
        z=1;
        for i=2:len
            if and(S(i)~=' ' , S(i)~='\n')
                z=z+1;
                for j=1:stateNum
                    for k=1:stateNum
                        p(z,j)=p(z,j)+p(z-1,k)*A(k,j,m)*E(j,S(i),m);
                    end
                end
            end
        end
            
        total=0;
        for j=1:stateNum
            total=total+(p(len,j));
        end
        likly=likly+log(total);
    end
    liklyArray(m)=likly/line_count;
end
[v_c,bestM]=max(liklyArray);
A_c=A(:,:,bestM);
E_c=E(:,:,bestM);
    
v_b=v_c;
A_b=A_c;
E_b=E_c;

%% step 1

for iteration=1:maxiter
    swap_c=zeros(5,1);
    swap_b=zeros(5,1);
    for k=1:testNo
        A(:,:,k)=A_c;
        E(:,:,k)=E_c;
    
        matrix=rand;
        if matrix<0.5
            firstrow=randi([1, stateNum]);
            secrow=randi([1, stateNum]);
            firstCol=randi([1, stateNum]);
            secCol=randi([1, stateNum]);
            extra=A(firstrow,firstCol,k);
            A(firstrow,firstCol,k)=A(secrow,secCol,k);
            A(secrow,secCol,k)=extra;
            swaps(k,2)=firstrow;
            swaps(k,3)=secrow;
            swaps(k,4)=firstCol;
            swaps(k,5)=secCol;
            swaps(k,1)=0;
        else
           firstrow=randi([1, stateNum]);
            secrow=randi([1, stateNum]);
           firstCol=randi([1, charNum]);
           secCol=randi([1, charNum]);
           extra=E(firstrow,firstCol,k);
            E(firstrow,firstCol,k)=E(secrow,secCol,k);
            E(secrow,secCol,k)=extra;
            swaps(k,2)=firstrow;
            swaps(k,3)=secrow;
            swaps(k,4)=firstCol;
            swaps(k,5)=secCol;
           swaps(k,1)=1;
        end
    end
    
    %% Normalization
    for k=1:testNo
        totala=sum(A(:,:,k),2);
        for row=1:stateNum
            for column=1:stateNum
                A(row,column,k)=A(row,column,k)/totala(row);
            end
        end
        totale=sum(E(:,:,k),2);
        for row=1:stateNum
            for column=1:charNum
                E(row,column,k)=E(row,column,k)/totale(row);
            end
        end
    end
    liklyArray=zeros(testNo,1);
    for m=1:testNo
        
    %% liklihood
        likly=0;
        for t=1:line_count
            S=seqCell{t};
            %% probGivenModel
            len=numel(S);
            p=zeros(len,stateNum);
    
            p(1,:)=StartMatrix(:).*(E(:,S(1),m));
            z=1;
            for i=2:len
                if and(S(i)~=' ' , S(i)~='\n')
                    z=z+1;
                    for j=1:stateNum
                        for k=1:stateNum
                            p(z,j)=p(z,j)+p(z-1,k)*A(k,j,m)*E(j,S(i),m);
                        end
                    end
                end
            end
    
            total=0;
            for j=1:stateNum
                total=total+(p(z,j));
            end
            likly=likly+log(total);
        end
        liklyArray(m)=likly/line_count;
    end  
    [sort_v,ind]=sort(liklyArray,'descend');

   %% step2
   flag=false;
   k=1;
   while flag==false
         current=sort_v(k);
         index=ind(k);
         presence=false;
         p=0;
         while and(presence==false  , p<tabusize)
            p=p+1;
            if tabu_list(p,:)==swaps(index,:)
                presence=true;
            end
         end
         if or(presence==false , current>v_b)
            flag=true;
            v_c=current;
            A_c=A(:,:,index);
            E_c=E(:,:,index);
            swap_c=swaps(index,:);
         elseif and(flag==false , k==testNo)
            swaps=zeros(testNo,5);
            for k=1:testNo
                A(:,:,k)=A_c;
                E(:,:,k)=E_c;
                matrix=rand;
               
                if matrix<0.5
                   firstrow=randi([1, stateNum]);
                   secrow=randi([1, stateNum]);
                   firstCol=randi([1, stateNum]);
                   secCol=randi([1, stateNum]);
                   extra=A(firstrow,firstCol,k);
                   A(firstrow,firstCol,k)=A(secrow,secCol,k);
                   A(secrow,secCol,k)=extra;
                   swaps(k,2)=firstrow;
                   swaps(k,2)=secrow;
                   swaps(k,4)=firstCol;
                   swaps(k,5)=secCol;
                   swaps(k,1)=0;
                else
                   firstrow=randi([1, stateNum]);
                   secrow=randi([1, stateNum]);
                   firstCol=randi([1, charNum]);
                   secCol=randi([1, charNum]);
                   extra=E(firstrow,firstCol,k);
                   E(firstrow,firstCol,k)=E(secrow,secCol,k);
                   E(secrow,secCol,k)=extra;
                   swaps(k,2)=firstrow;
                   swaps(k,2)=secrow;
                   swaps(k,4)=firstCol;
                   swaps(k,5)=secCol;
                   swaps(k,1)=1;
                end
             end
            %% Normalization
            for k=1:testNo
                totala=sum(A(:,:,k),2);
                for row=1:stateNum
                    for column=1:stateNum
                        A(row,column,k)=A(row,column,k)/totala(row);
                    end
                end
                totale=sum(E(:,:,k),2);
                for row=1:stateNum
                    for column=1:charNum
                        E(row,column,k)=E(row,column,k)/totale(row);
                    end
                end
            end
            liklyArray=zeros(testNo,1);
            for m=1:testNo
        
            %% liklihood
                likly=0;
                for t=1:line_count
                    S=seqCell{t};
                    %% probGivenModel
                    len=numel(S);
                    p=zeros(len,stateNum);
    
                    p(1,:)=StartMatrix(:).*(E(:,:,S(1),m));
                    z=1;
                    for i=2:len
                        if and(S(i)~=' ' , S(i)~='\n')
                            z=z+1;
                            for j=1:stateNum
                                for k=1:stateNum
                                    p(z,j)=p(z,j)+p(z-1,k)*A(k,j,m)*E(j,S(i),m);
                                end
                            end
                        end
                    end
        
                    total=0;
                    for j=1:stateNum
                        total=total+(p(z,j));
                    end
                    likly=likly+loga(total);
                end
                liklyArray(m)=likly/line_count;
            end           
            [sort_v,ind]=sort(liklyArray,'descend');
            k=1;
         else
                k=k+1;
         end
   end
 %% step3           
   if v_c<v_b
      v_c=v_b;
      %best_iter=iteration;
      A_c=A_b;
      E_c=E_b;
      swap_c=swap_b;
   else
       v_b=v_c;
       A_b=A_c;
       E_b=E_c;
       swap_b=swap_c;
       
       
   end
    tabu_list(insertion,:)=swap_c(:);
        
    insertion=insertion+1;
    if insertion>tabusize
       insertion=1;
    end
end

end
   
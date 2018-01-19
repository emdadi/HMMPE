function testProb=testprog(A,E,testCell,line_count)
stateNum=size(A,1);
StartMatrix=(double(1)/stateNum)*ones(stateNum,1);
ehtem=0;
%line_count=numel(testSymbols);
for P=1:size(testCell,1)
    S=testCell{P};
    
    %% probGivenModel
    len=numel(S);
    
    p=zeros(len,stateNum);
    p(1,:)=StartMatrix(:).*E(:,int32(S(1)));
    z=1;
    
    for i=2:len
            if and(S(i)~=' ' , S(i)~='\n')
                z=z+1;
                for j=1:stateNum
                    for k=1:stateNum
                        p(z,j)=p(z,j)+p(z-1,k)*A(k,j)*E(j,int32(S(i)));
                    end
                end
            end
    end
    total2=sum(p(z,:));
    prob2=log(total2);
    
    ehtem=ehtem+prob2;
    
end

testProb=ehtem/(line_count);
end
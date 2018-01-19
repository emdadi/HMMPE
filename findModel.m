function b=findModel(model,i,j)
b=false;
for k=1:numel(model.tran)
    if and(model.tran(k,1)==i , model.tran(k,2)==j)
        b=true;
    end
end
end

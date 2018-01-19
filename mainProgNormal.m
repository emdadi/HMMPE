% In God we trust
% AntMarkov
% Date: 18 mehr 96:
% time: 12:17

clc;
clear;
close all;

%% open file
resultFile = fopen('resultNormal.txt','w');
baumMatrix = fopen('baumMatrixNormal.txt','w');
antMatrix = fopen('antMatrixNormal.txt','w');
SNNMFMatrix=fopen('SNNMFMatrixNormal.txt','w');
vitMatrix = fopen('vitMatrixNormal.txt','w');
tabuMatrix = fopen('tabuMatrixNormal.txt','w');
test_no=1;


fprintf(resultFile,'test_no \t AvgOrg \t AntProb \t baumProb \t vitProb  \t SNNMFProb \t tabuProb');
fprintf(resultFile,'\t antAnorm2 \t baumAnorm2 \t vitAnorm2  \t SNNMFAnorm2 \t tabuAnorm2 ');
fprintf(resultFile,'\t antEnorm2 \t baumEnorm2 \t vitEnorm2 \t SNNMFEnorm2 \t tabuEnorm2 ');
fprintf(resultFile,'\t antAnormInf \t baumAnormInf \t vitAnormInf \t SNNMFAnormInf \t tabuAnormInf ');
fprintf(resultFile,'\t antEnormInf \t baumEnormInf \t vitEnormInf  \t SNNMFEnormInf \t tabuEnormInf ');
fprintf(resultFile,'\t antAnormFro \t baumAnormFro \t vitAnormFro  \t SNNMFAnormFro \t tabuAnormFro ');
fprintf(resultFile,'\t antEnormFro \t baumEnormFro \t vitEnormFro   \t SNNMFEnormFro \t tabuEnormFro ');
fprintf(resultFile,'\t antTime \t baumTime \t vitTime  \t SNNMFTime \t tabuTime ');
fprintf(resultFile,'\t antItr \t baumItr \t vitItr \t issparse(A) \t issparse(E)\n ');

%% Problem Definition
for count1=200:432
    test_no=count1
    
    stre=strcat('emission_normal_Complete_',num2str(count1));
    strte=strcat(stre,'.txt');
       
    strt=strcat('test_normal_Complete_',num2str(count1));
    strtt=strcat(strt,'.txt');
    
    strs=strcat('seq_normal_Complete_',num2str(count1));
    strts=strcat(strs,'.txt');
    
    strtr=strcat('transition_normal_Complete_',num2str(count1));
    strttr=strcat(strtr,'.txt');
    
    
    streR=strcat('emission_random_normal_Complete_',num2str(count1));
    strteR=strcat(streR,'.txt');
    
    strtrR=strcat('transition_random_normal_Complete_',num2str(count1));
    strttrR=strcat(strtrR,'.txt');
    
    
    seqCell = dlmread(strts);
    testCell=dlmread(strtt);
    Aorigin=importdata(strttr);
    Eorigin=importdata(strte);
    guessE=importdata(strteR);
    guessTR=importdata(strttrR);
    

    characters=unique(seqCell);
    characters(characters==0)=[];
    characters = characters(~isnan(characters)) ;
    charNum=numel(characters);
    stateNum=size(Aorigin,1);
    line_count=size(seqCell,1);
    
    e=double(1)/stateNum;
    Startorigin=e*ones(stateNum,1);

    
    characterspc=unique(seqCell);
    characterspc(characterspc==0)=[];
    characterspc = characterspc(~isnan(characterspc)) ;
    charNumSpc=numel(characterspc);
    
    fprintf(resultFile,num2str(count1));
    fprintf(resultFile,'\t');

%% edit seq,test
    seqCell = arrayfun(@(x) seqCell(x,(seqCell(x,:)~=0)), 1:size(seqCell,1), 'uni', 0);

    if ~iscell(seqCell) || ischar(seqCell{1})
            [~, seqCell]  = ismember(seqCell,characters);
            %seqCell
            if any(seqCell(:)==0)
                %seqCell
                error(message('stats:hmmtrain:MissingSymbol'));
            end
    else  % now deal with a cell array of sequences
    
            numSeqs = numel(seqCell);
            newSeqs = cell(numSeqs,1);
            for count = 1:numSeqs
                [~, newSeqs{count}] = ismember(seqCell{count},characters);
                if any(newSeqs{count}(:)==0)
                    error(message('stats:hmmtrain:MissingSymbol'));
                end
            end
            for count = 1:numSeqs
                newSeqs{count}=int32(newSeqs{count});
            end
            seqCell = newSeqs;
    end
    
    
    
    testCell = arrayfun(@(x) testCell(x,(testCell(x,:)~=0)), 1:size(testCell,1), 'uni', 0);
    if ~iscell(testCell) || ischar(testCell{1})
            [~, testCell]  = ismember(testCell,characters);
            %seqCell
            if any(testCell(:)==0)
                %seqCell
                error(message('stats:hmmtrain:MissingSymbol'));
            end
    else  % now deal with a cell array of sequences
        numtest = numel(testCell);
        newtest = cell(numtest,1);
        for count = 1:numtest
            [~, newtest{count}] = ismember(testCell{count},characters);
            if any(newtest{count}(:)==0)
                error(message('stats:hmmtrain:MissingSymbol'));
            end
        end
        for count = 1:numSeqs
                newtest{count}=int32(newtest{count});
        end
        testCell = newtest;
    end
    
    
    
    
    seqCellNNMF=[];
    for i=1:size(seqCell,1)
        for j=1:numel(seqCell{i})-1
            
            a=seqCell{i}(j:j+1);
            seqCellNNMF=[seqCellNNMF, a.'];
        end
    end
%% Avg Original
    sumtest=0;
    for P=1:size(testCell,1)
        S=testCell{P};
        %% probGivenModel
        len=numel(S);
        p=zeros(len,stateNum);
        p(1,:)=Startorigin(:).*Eorigin(:,int32(S(1)));
        z=1;
        for i=2:len
            if and(S(i)~=' ' , S(i)~='\n')
               z=z+1;
               for j=1:stateNum
                   for k=1:stateNum
                       p(z,j)=p(z,j)+p(z-1,k)*Aorigin(k,j)*Eorigin(j,int32(S(i)));
                   end
               end
            end
        end
        total=0;
        for j=1:stateNum
            total=total+p(z,j);
        end
        prob1=log(total);
        sumtest=sumtest+prob1;
    end
    
    AvgOriginTest=sumtest/line_count;
    fprintf(resultFile,num2str(AvgOriginTest));
    fprintf(resultFile,'\t');
    %% calling programs
    tic
    [antA,antE,antItr]=mainAntMarkov(seqCell,stateNum,line_count,characters);
    antTime=toc;
  
    tic
    [baumA,baumE,baumItr]=hmmtrain(seqCell,guessTR,guessE);
    baumTime=toc;
    tic
    [vitA,vitE,vitItr,stop]=hmmtrain(seqCell,guessTR,guessE,'Algorithm','Viterbi');
    vitTime=toc;
    if stateNum<=charNum
        tic
        [SNNMFE,SNNMFA]=learnSNNMF(seqCellNNMF,stateNum,charNum);
        SNNMFE=SNNMFE.';
        SNNMFTime=toc;
    else
        SNNMFE=antE;
        SNNMFA=antA;
        SNNMFTime=Inf;
    end
    if mod(count1,18)==0
        tabuE=antE;
        tabuA=antA;
        tabuTime=Inf;
    else
        tic
    [tabuA,tabuE]=tabu(seqCell,stateNum,line_count,characters);
    tabuTime=toc;
        
    end
    
    %% print A,E
    fprintf(antMatrix,num2str(test_no));
    fprintf(antMatrix,'\n A= \n');
    for i=1:stateNum
         for j=1:stateNum
            fprintf(antMatrix,num2str(antA(i,j)));
            fprintf(antMatrix,'\t');
         end
         fprintf(antMatrix,'\n');
    end
    fprintf(antMatrix,'\n\n E= \n');
    for i=1:stateNum
         for j=1:charNum
            fprintf(antMatrix,num2str(antE(i,j)));
            fprintf(antMatrix,'\t');
         end
         fprintf(antMatrix,'\n');
    end
    fprintf(antMatrix,'\n*******************\n');
    
    
    fprintf(baumMatrix,num2str(test_no));
    fprintf(baumMatrix,'\n A= \n');
    for i=1:stateNum
         for j=1:stateNum
            fprintf(baumMatrix,num2str(baumA(i,j)));
            fprintf(baumMatrix,'\t');
         end
         fprintf(baumMatrix,'\n');
     end
    fprintf(baumMatrix,'\n \n E= \n');
    for i=1:stateNum
         for j=1:charNum
            fprintf(baumMatrix,num2str(baumE(i,j)));
            fprintf(baumMatrix,'\t');
         end
         fprintf(baumMatrix,'\n');
    end
    fprintf(baumMatrix,'\n*******************\n');
    
    fprintf(vitMatrix,num2str(test_no));
    fprintf(vitMatrix,'\n A= \n');
    for i=1:stateNum
         for j=1:stateNum
            fprintf(vitMatrix,num2str(vitA(i,j)));
            fprintf(vitMatrix,'\t');
         end
         fprintf(vitMatrix,'\n');
     end
    fprintf(vitMatrix,'\n \n E= \n');
    for i=1:stateNum
         for j=1:charNum
            fprintf(vitMatrix,num2str(vitE(i,j)));
            fprintf(vitMatrix,'\t');
         end
         fprintf(vitMatrix,'\n');
    end
    fprintf(vitMatrix,'\n*******************\n');
    
    
    
    fprintf(SNNMFMatrix,num2str(test_no));
    fprintf(SNNMFMatrix,'\n A= \n');
     for i=1:stateNum
         for j=1:stateNum
           fprintf(SNNMFMatrix,num2str(SNNMFA(i,j)));
                fprintf(SNNMFMatrix,'\t');
         end
         fprintf(SNNMFMatrix,'\n');
     end
    fprintf(SNNMFMatrix,'\n \n E= \n');
    for i=1:stateNum
         for j=1:charNum
            fprintf(SNNMFMatrix,num2str(SNNMFE(i,j)));
            fprintf(SNNMFMatrix,'\t');
         end
        fprintf(SNNMFMatrix,'\n');
    end
    fprintf(SNNMFMatrix,'\n*******************\n');
    
    
    fprintf(tabuMatrix,num2str(test_no));
    fprintf(tabuMatrix,'\n A= \n');
    
    tabuAA=zeros(stateNum,stateNum);
    tabuEE=zeros(stateNum,charNum);
     for i=1:stateNum
         for j=1:stateNum
            tabuAA(i,j)=tabuA(i,j,1);
            fprintf(tabuMatrix,num2str(tabuAA(i,j)));
            fprintf(tabuMatrix,'\t');
         end
         fprintf(tabuMatrix,'\n');
     end
    fprintf(tabuMatrix,'\n \n E= \n');
    for i=1:stateNum
         for j=1:charNum
             
            tabuEE(i,j)=tabuE(i,j,1);
            fprintf(tabuMatrix,num2str(tabuEE(i,j)));
            fprintf(tabuMatrix,'\t');
         end
         fprintf(tabuMatrix,'\n');
    end
    fprintf(tabuMatrix,'\n*******************\n');
%% calling test

    antProb=testprog(antA,antE,testCell,line_count);
    if stop==0
        vitProb=-Inf;
    else
    vitProb=testprog(vitA,vitE,testCell,line_count);
    end
    baumProb=testprog(baumA,baumE,testCell,line_count);
    
    if tabuTime==Inf
        tabuProb=-Inf
    else
    tabuProb=testprog(tabuA,tabuE,testCell,line_count);
    end
    
    
    if SNNMFTime==Inf
        SNNMFProb=-Inf
    else
    SNNMFProb=testprog(SNNMFA,SNNMFE,testCell,line_count);
    end
%% print liklihhod
    
    
    fprintf(resultFile,num2str(antProb));
    fprintf(resultFile,'\t');
    
    fprintf(resultFile,num2str(baumProb));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(vitProb));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFProb));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(tabuProb));
    fprintf(resultFile,'\t');   
%% calculate norm

    antC=antA-Aorigin;
    antD=antE-Eorigin;
    
    antAnorm2=norm(antC);
    antAnormInf=norm(antC,'inf');
    antAnormFro=norm(antC,'fro');
    
    antEnorm2=norm(antD);
    antEnormInf=norm(antD,'inf');
    antEnormFro=norm(antD,'fro');
    
    baumC=baumA-Aorigin;
    baumD=baumE-Eorigin;
    baumAnorm2=norm(baumC);
    baumAnormInf=norm(baumC,'inf');
    baumAnormFro=norm(baumC,'fro');

    baumEnorm2=norm(baumD);
    baumEnormInf=norm(baumD,'inf');
    baumEnormFro=norm(baumD,'fro');

    vitC=vitA-Aorigin;
    vitD=vitE-Eorigin;

    vitAnorm2=norm(vitC);
    vitAnormInf=norm(vitC,'inf');
    vitAnormFro=norm(vitC,'fro');

    vitEnorm2=norm(vitD);
    vitEnormInf=norm(vitD,'inf');
    vitEnormFro=norm(vitD,'fro');

    
    
    SNNMFC=SNNMFA-Aorigin;
    SNNMFD=SNNMFE-Eorigin;

    SNNMFAnorm2=norm(SNNMFC);
    SNNMFAnormInf=norm(SNNMFC,'inf');
    SNNMFAnormFro=norm(SNNMFC,'fro');

    SNNMFEnorm2=norm(SNNMFD);
    SNNMFEnormInf=norm(SNNMFD,'inf');
    SNNMFEnormFro=norm(SNNMFD,'fro');
    
    
    tabuC=tabuAA-Aorigin;
    tabuD=tabuEE-Eorigin;
    tabuAnorm2=norm(tabuC);
    tabuAnormInf=norm(tabuC,'inf');
    tabuAnormFro=norm(tabuC,'fro');

    tabuEnorm2=norm(tabuD);
    tabuEnormInf=norm(tabuD,'inf');
    tabuEnormFro=norm(tabuD,'fro');
    
    %% print norm2
    fprintf(resultFile,num2str(antAnorm2));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumAnorm2));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitAnorm2));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFAnorm2));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(tabuAnorm2));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(antEnorm2));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumEnorm2));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitEnorm2));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFEnorm2));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(tabuEnorm2));
    fprintf(resultFile,'\t');
%% print normInf
    fprintf(resultFile,num2str(antAnormInf));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumAnormInf));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitAnormInf));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFAnormInf));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(tabuAnormInf));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(antEnormInf));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumEnormInf));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitEnormInf));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFEnormInf));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(tabuEnormInf));
    fprintf(resultFile,'\t');
%% print normFro
    fprintf(resultFile,num2str(antAnormFro));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumAnormFro));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitAnormFro));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFAnormFro));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(tabuAnormFro));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(antEnormFro));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumEnormFro));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitEnormFro));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFEnormFro));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(tabuEnormFro));
    fprintf(resultFile,'\t');
%% print time
    fprintf(resultFile,num2str(antTime));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumTime));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitTime));
    fprintf(resultFile,'\t');

    fprintf(resultFile,num2str(SNNMFTime));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(tabuTime));
    fprintf(resultFile,'\t');
    
%%print Iteration
    fprintf(resultFile,num2str(antItr));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(baumItr));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(vitItr));
    fprintf(resultFile,'\t');
    
    fprintf(resultFile,num2str(issparse(Aorigin)));
    fprintf(resultFile,'\t');
    fprintf(resultFile,num2str(issparse(Eorigin)));
    fprintf(resultFile,'\t');
    
    fprintf(resultFile,'\n');
    
    
    %% close files
end
fclose(resultFile);
fclose(antMatrix);
fclose(baumMatrix);
fclose(vitMatrix);

fclose(tabuMatrix);
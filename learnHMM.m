function [b1hat, binfhat, Bxhat, EMIShat, PINIThat, TRANShat]...
    = learnHMM(seq,k,d)
%LEARNHMM 
%Implementation of a learning algorithm as proposed by [1].
%Calculates estimates for a hidden Markov model parameter
%representation. Furthermore the algorithm calculates estimates for 
%the emission probability matrix, the initial state distribution
%and state transition probability matrix of a hidden Markov model
%
% [b1hat, binfhat, Bxhat, EMIShat, PINIThat, TRANShat]...
%    = learnHMM(seq,k,d);
%
% Inputs:   seq     sequences of triples of observations
%           k       number of states
%           d       number of observations
% Returns:  HMM model parameterized by b^_1, b^_inf, B^_x
%                                   for all x in [n]
%           EMIShat     estimator for observation probability matrix
%           PINIThat    estimator for initial state distribution
%           TRANShat    estimator for state transition probability matrix
%
% Please note: estimators for matrices are not explicitly stated in
% the algorithm, calculation is based on a method by [2].
% Implementation as scetched in [1], Appendix C.
%
% [1] Hsu, Daniel, Sham M. Kakade, and Tong Zhang. 
% "A spectral algorithm for learning hidden Markov models." 
% Journal of Computer and System Sciences 78.5 (2012): 1460-1480.
% [2] Mossel, Elchanan, and Sébastien Roch. 
% "Learning nonsingular phylogenies and hidden Markov models." 
% Proceedings of the thirty-seventh annual 
%      ACM symposium on Theory of computing. ACM, 2005.
%
% Carl Mattfeld, 2014
% Version 2014-03-08

N=length(seq);

% =================================================
% LearnHMM Step 1 - calculating empirical averages
% =================================================
% Vector and matrix quantities for HMM representation
% [P_1]_i = Pr[x_1=i]
% [P_2,1]_ij = Pr[x_2=i,x_1=j]
% [P_3,x,1]_ij = Pr[x_3=i,x_2=x,x_1=j] for all x in [n]
% where P_1 n-vector, P_2,1 n-by-n matrix and P_3,x,1 n-by-n matrix

% Calculate P1^
% =============
P1hat=hist(seq(:,1),d)./N;
P1hat=transpose(P1hat);

% Calculate P21^
% ==============
P21hat=zeros(d);
% count occurances of sequence
for i=1:N
    P21hat(seq(i,2),seq(i,1))=P21hat(seq(i,2),seq(i,1))+1;
end
P21hat=P21hat./N;

% Calculate P3x1^
% ===============
P3x1hat=zeros(d,d,d);

% count occurances of sequence
for i=1:N
    P3x1hat(seq(i,3),seq(i,1),seq(i,2))=P3x1hat(seq(i,3),seq(i,1),seq(i,2))+1;
end
P3x1hat = P3x1hat./N;

% =================================================
% LearnHMM Step 2 - computing the SVD
% =================================================

 [Uhat,~,~]=svds(P21hat,k);

 % =================================================
% LearnHMM Step 3 - compute model paramters
% =================================================

b1hat = transpose(Uhat)*P1hat;  % calculate b^_1
binfhat = pinv(transpose(P21hat)*Uhat)*P1hat; % calculate b^_inf
% calculate B_x^ for all x in [n]
Bxhat=zeros(k,k,d);
for i=1:d
    Bxhat(:,:,i) = transpose(Uhat)*P3x1hat(:,:,i)*...
        pinv(transpose(Uhat)*P21hat);
end

% ======================================================
% LearnHMM - recover observation and transition matrices
% ======================================================
% cf HKZ '12, Appendix C

% Calculate P31^
% [P_3,1]_ij = Pr[x_3=i,x_1=j]

P31hat=zeros(d);
% count occurances of sequence
for i=1:N
    P31hat(seq(i,3),seq(i,1))=P31hat(seq(i,3),seq(i,1))+1;
end
P31hat=P31hat./N;

% for all x in [n]:
% (U^t * P_3,x,1)(U^t P_3,1)^+ 
% =(U^t O T)O_x(U^t O T)^-1
% O_x is diagonal, eigenvalues are thus the observation probabilities
% O_(r,1),...,O_(r,m) for r in [n]

% calculate matrix to be decomposed 
Decomp=zeros(k,k,d);
for i=1:d
    Decomp(:,:,i)=(transpose(Uhat)*P3x1hat(:,:,i))*...
        pinv(transpose(Uhat)*P31hat);
end

% O_x for an x in [n] corresponds to the observation probabilities
% Decompose matrices
EMIShat=zeros(k,d);
for i=1:d
    EMIShat(:,i)=eig(Decomp(:,:,i));
end

EMIShat=transpose(EMIShat);

% Calculate initial state distribution
PINIThat  =pinv(EMIShat)*P1hat;

% Calculate transition matrix
TRANShat = pinv(EMIShat)*P21hat*...
    transpose(pinv(EMIShat))*diag(PINIThat)^(-1);

end
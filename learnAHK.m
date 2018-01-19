function [EMIShat, TRANShat] = learnAHK(seq,k,d)
%LEARNAHK 
%Implementation of a learning algorithm as proposed by [1].
%Calculates numerical estimates for the emission probability
%matrix and state transition matrix of a hidden Markov model, provided
%seqences of observations (grouped in triples).
%
% [EMIShat, TRANShat] = learnAHK(seq,k,d);
%
% INPUT:    seq     sequences of triples of observations
%           k       number of states
%           d       number of observations
% OUTPUT:   EMIShat     estimate for emission probability matrix
%           TRANShat    estimate for state transisiton probability matrix
%
% 
% [1] Anandkumar, Animashree, Daniel Hsu, and Sham M. Kakade. 
%  "A method of moments for mixture models and hidden Markov models." 
%  arXiv preprint arXiv:1203.0683 (2012).
%
% Carl Mattfeld, 2014
% Version 2014-03-08

% Algorithm B - calculate M2
% ===========================
% Steps corresponding to formulation of Algorithm B in original work

N=length(seq);

% Step 1 - obtain empirical averages
% ==============================================
% by counting occurances in sequence
% P31hat_ij = Pr[x_3=i,x_1=j]
% P312hat_ijk = Pr[x_3=i,x_1=j,x_2=k]

% Calculate P31hat
P31hat=zeros(d);
for i=1:N
    
    P31hat(seq(i,3),seq(i,1))=P31hat(seq(i,3),seq(i,1))+1;
end
P31hat=P31hat./N;

% Calculate P32hat
P32hat=zeros(d);
for i=1:N
    P32hat(seq(i,3),seq(i,2))=P32hat(seq(i,3),seq(i,2))+1;
end
P32hat=P32hat./N;

% Calculate P312hat
P312hat=zeros(d,d,d);
for i=1:N
    P312hat(seq(i,3),seq(i,1),seq(i,2))=P312hat(seq(i,3),seq(i,1),seq(i,2))+1;
end
P312hat = P312hat./N;

% Step 2 - calculate singular vectors
% ==============================================
[U3hat,~,U1hat]=svds(P31hat,k);
[~,~,U2hat]=svds(P32hat,k);

% Step 3 - This is the long one
% ==============================================
% ==============================================
% "Pick an invertible matrix THETA with its i-th row denoted as
% transpose(THETA(:,i)). In the absence of any prior information about M_3 
% (M_2), a suitable choice for THETA is a random rotation matrix


% Pick random rotation matrix
[THETA, R] = qr(randn(k));
THETA = THETA*diag(sign(diag(R))); 

% exchange two columns in case det(THETA)=-1
if det(THETA)<0
    columntemp = THETA(:,1);
    THETA(:,1) = THETA(:,2);
    THETA(:,2) = columntemp;
    clear columntemp;
end

% Form the Matrix B312hat(U2hat*THETA_1)
% ======================================
% Calculate P123hat(U3hat THETA_j) for j in {1,...,m}
% all these matrices calculated now, since needed later
ARG=zeros(d,k);
P312hatARG = zeros(d,d,k);

% refer to section 2.2 in AHK '12
for j=1:k
    ARG(:,j) = U2hat*transpose(THETA(j,:));
    for i=1:d
        P312hatARG(:,:,j) = P312hatARG(:,:,j) +...
            ARG(i,j)*P312hat(:,:,i);
    end
end

% Calculate B123hat (already for all U3hat THETA_j, j in {1,..,m}, since
% needed later

B312hat = zeros(k,k,k);

for j=1:k
    B312hat(:,:,j) = (transpose(U3hat)*P312hatARG(:,:,j)*U1hat)...
        /(transpose(U3hat)*P31hat*U1hat);
end

% Compute R1hat that diagonalizes B123hat(U3hat THETA_1)
diagonals = zeros(k,k,k);
[R3hat,diagonals(:,:,1)] = eig(B312hat(:,:,1));

% Step 4 - form matrix Lhat
% ==========================
% "For each i in {2,...,m} obtain the diagonal entries
% \lambda_i1,...,\lambda_im of R3hat^(-1)*B312hat(U2 THETA_i)*R3hat and
% form the matrix Lhat whose (i,j)-th entry is \lambda_ij

Lhat = zeros(k,k);
Lhat(1,:) = diag(diagonals(:,:,1));

for j=2:k
    diagonals(:,:,j) = ((R3hat)\B312hat(:,:,j))*R3hat;
    Lhat(j,:) = diag(diagonals(:,:,j));
end

% Step 5 - return M2
% ========================
% estimator for M2 (M2 = EMIS*TRANS = OT)
M2hat = U2hat*(THETA\Lhat);

EMIShat = real(M2hat);

% Calculate TRANShat
TRANShat = (transpose(U3hat)*EMIShat)\R3hat;

% scale columns
TRANShat = real(TRANShat);
TRANShat = bsxfun(@rdivide, TRANShat, sum(TRANShat,1));

end
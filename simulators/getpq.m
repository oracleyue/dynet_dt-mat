function [Q,P,D]=getpq(A,B,p)
% Find Q and P
% [Q,P]=getpq(A,B,p)  where p is the size of the output.
% C must be = [I 0]

%clear

%sss=3;
%syst

n=length(A);
%m=min(size(B));
%p=min(size(C));
C=[eye(p) zeros(p,n-p)];

A11=A(1:p,1:p);
A12=A(1:p,p+1:n);
A21=A(p+1:n,1:p);
A22=A(p+1:n,p+1:n);

B1=B(1:p,:);
B2=B(p+1:n,:);

syms s
if n>p
  W=simplify(A11+A12*inv(s*eye(n-p)-A22)*A21);
  D=simplify(diag(diag(W)));  
  V=simplify(A12*inv(s*eye(n-p)-A22)*B2+B1);
  Q=simplify(inv(s*eye(p)-D)*(W-D));
  P=simplify(inv(s*eye(p)-D)*V);
else
  W=A11;
  D=diag(diag(W));
  V=B1;
  Q=inv(s*eye(p)-D)*(W-D);
  P=inv(s*eye(p)-D)*V;
end

  
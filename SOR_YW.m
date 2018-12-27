%SOR_Method
%A: a nXn matrix
%b: a nX1 vector
%x: a solution of Ax=b
function x=SOR_YW(A,b,omega,x0,options)
D=diag(A);
L=tril(A)-diag(D);
U=triu(A)-diag(D);
CurIter=0;
n=size(A,1);
x=x0;
while 1
    for i=1:n
        x(i)=omega*(b(i)-L(i,:)*x-U(i,:)*x0)/D(i)...
             +(1-omega)*x0(i);
    end
    CurIter=CurIter+1;
    disp([num2str(CurIter) '-th: Residual is ' num2str(norm(A*x-b))])
    if CurIter>options.MaxIter
        break
    end
    if norm(x-x0)<options.Tol
        break
    end
    if norm(A*x-b)<options.Tol
        break
    end
    x0=x;
end
%Jacobi_Method
%A: a nXn matrix
%b: a nX1 vector
%x: a solution of Ax=b
function x=Jacobi_YW(A,b,x0,options)
D=diag(A);
L=tril(A)-diag(D);
U=triu(A)-diag(D);
CurIter=0;
while 1
    x=-(L+U)*x0./D+b./D;
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
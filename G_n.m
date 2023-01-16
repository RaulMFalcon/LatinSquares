function p=G_n(x,y,n)
    [r,p]=polynomialReduce(F_n(x,n)-F_n(y,n),x-y);
end

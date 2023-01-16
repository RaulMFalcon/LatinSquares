function [List_L]=LS(P)
    n=size(P,1);
    syms x [n n];
    s=0;
    for i=1:n
        for j=1:n
            if P(i,j)==0
                s=s+1;
                I(s)=F_n(x(i,j),n);
                for k=1:n
                    if (k~=i)
                        if P(k,j)==0
                            s=s+1;
                            I(s)=G_n(x(i,j),x(k,j),n);
                        else
                            s=s+1;
                            I(s)=G_n(x(i,j),P(k,j),n);
                        end
                    end
                    if (k~=j)
                        if P(i,k)==0
                            s=s+1;
                            I(s)=G_n(x(i,j),x(i,k),n);
                        else
                            s=s+1;
                            I(s)=G_n(x(i,j),P(i,k),n);
                        end
                    end
                end
            end
        end
    end
    gb=gbasis(I);
    if gb~=1
         SOL=table2array(struct2table(solve(gb)));
         for l=1:size(SOL,1)
            Q=P;
            s=0;
            for i=1:n
                for j=1:n
                    if Q(i,j)==0
                        s=s+1;
                        Q(i,j)=SOL(l,s);
                    end
                end
            end
            List_L{l}=Q;
         end
    else
        List_L={};
    end
end

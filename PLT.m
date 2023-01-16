function [List_Q]=PLT(P,T)
    List_Q={};
    n=size(P,1);
    syms x [n n];
    s=0;
    for i=1:n
        if P(T(i,1),T(i,2))==0
            s=s+1;
            I(s)=F_n(x(T(i,1),T(i,2)),n);
            for j=1:n
                if j~=i
                    s=s+1;
                    if P(T(j,1),T(j,2))==0
                        if j>i
                            I(s)=G_n(x(T(i,1),T(i,2)),x(T(j,1),T(j,2)),n);    
                        end
                    else
                        I(s)=G_n(x(T(i,1),T(i,2)),P(T(j,1),T(j,2)),n);    
                    end
                end
            end
            for j=1:n
                if P(T(i,1),j)>0
                    s=s+1;
                    I(s)=G_n(x(T(i,1),T(i,2)),P(T(i,1),j),n);   
                end
                if P(j,T(i,2))>0
                    s=s+1;
                    I(s)=G_n(x(T(i,1),T(i,2)),P(j,T(i,2)),n);   
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
            for j=1:n
                for k=1:n
                    if T(k,1)==j
                        if Q(j,T(k,2))==0
                            s=s+1;
                            Q(j,T(k,2))=SOL(l,s);
                        end
                    end
                end
            end
            List_Q{l}=Q;
         end
    else
        List_Q={};
    end
end

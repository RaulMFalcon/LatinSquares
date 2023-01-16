function [L] = HL(P)
    L={};
    L1={};
    n=size(P,1);
    syms x [n n];
    for i=1:n
        T(i,1)=P(P(i,i),P(i,i));
        T(i,2)=P(i,i);
    end
    C=PLT(P,T);
    for i=1:size(C,2)
        L{size(L,2)+1}={C{i},T};
    end
    while size(L,2)>0
        Q=L{1}{1};
        T=L{1}{2};
        s=0;
        L(1)=[];
        for i=1:n
            T(i,1)=Q(T(i,1),P(i,i));
        end
        for i=1:n
            if Q(T(i,1),T(i,2))==0
                s=1;
                break
            end
        end
        if s==0
            L1{size(L1,2)+1}=Q;
        else
            C=PLT(Q,T);
            if size(C,2)>0
                for i=1:size(C,2)
                    L{size(L,2)+1}={C{i},T};
                end
            else
                L1{size(L1,2)+1}=Q;
            end
        end
    end
    for i=1:size(L1,2)
        L2=LS(L1{i});
        for j=1:size(L2,2)
            L{size(L,2)+1}=L2{j};
        end
    end
end

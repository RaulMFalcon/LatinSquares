function [List_L]=LS(P)
% LS   Generate the set of Latin squares containing a partial Latin square.
%    L = LS(P) returns the list of Latin squares to which the partial Latin
%    square P is completable.
%
%    [1] V. Álvarez, J.A. Armario, R.M. Falcón, M.D. Frau, F. Gudiel and
%        M.B. Güemes. A computational approach to analyze the Hadamard 
%        quasigroup product. Submitted, 2023.
% 
%    Víctor Álvarez, José Andrés Armario, Raúl M. Falcón, 
%    María Dolores Frau, Felix Gudiel and María Belén Güemes.
%    January 16, 2023
%    Dpt. Applied Mathematics I.
%    University of Seville, Spain.

% Construct an ideal associated to the partial Latin square P, whose
% algebraic set can be identified with the required list of Latin squares.
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
% Obtain the Gröbner basis of the ideal. If it is 1, then there is no
% solution. Otherwise, the list of Latin squares is obtained.
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

function L = HL(P)
% HL   Generate the set of Latin squares to which a given partial Latin 
%      square is completable and whose successive Hadamard quasigroup 
%      products preserve the Latin square property.
%    L = HL(P) returns the list of Latin squares to which the partial Latin
%    square P is completable. Empty cells in P are represented by zeros.
%
%    Example: celldisp(HL([2 0 0 0; 0 1 0 0; 0 0 4 0; 0 0 0 3]))
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

% Construct the first transversal that must be filled to get a Latin square
% satisfying the required conditions.
    L={};
    L1={};
    n=size(P,1);
    syms x [n n];
    for i=1:n
        T(i,1)=P(P(i,i),P(i,i));
        T(i,2)=P(i,i);
    end
 % Generate the list of partial Latin squares that result after filling the
 % cells of the previous transversal.
    C=PLT(P,T);
 % Associate each partial Latin square in the list to the previous 
 % transversal.
    for i=1:size(C,2)
        L{size(L,2)+1}={C{i},T};
    end
 % Generate successive localizated transversals that must be filled to get
 % the required Latin squares. Then, fill these transversals.
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
            end
        end
    end
% Generate the Latin squares to which each obtained partial Latin
% square is completable. Each one of these Latin squares satisfies the
% required conditions.
    for i=1:size(L1,2)
        L2=LS(L1{i});
        for j=1:size(L2,2)
            L{size(L,2)+1}=L2{j};
        end
    end
end

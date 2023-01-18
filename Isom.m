function L=Isom(L1,L2)
% Isom  Compute the set of isomorphisms between two Latin squares.
%    L = Isom(L1,L2) returns the set of isomorphism between two Latin 
%    squares L1 and L2.
%
%    Example: 
%    celldisp(Isom([2 4 3 1; 3 1 2 4; 1 3 4 2; 4 2 1 3],[2 3 1 4; 4 1 3 2; 3 2 4 1; 1 4 2 3]))
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

% Construct an ideal associated to the Latin squares L1 and L2, whose 
% algebraic set can be identified with the set of isomorphism between L1
% and L2.
    n=size(L1,1);
    syms x [n];
    s=0;
    for i=1:n
        s=s+1;
        I(s)=F_n(x(i),n);
        for j=1:n
            if i~=j
                s=s+1;
                I(s)=G_n(x(i),x(j),n);
            end
            for k=1:n
                for l=1:n
                    s=s+1;
                    I(s)=F_n(x(i),n)*F_n(x(j),n)*(x(L1(i,j))-L2(k,l))/((x(i)-k)*(x(j)-l));
                end
            end
        end
    end
% Obtain the Gröbner basis of the ideal. If it is 1, then there is no
% solution. Otherwise, the list of isomorphisms is obtained.
    gb=gbasis(I);
    if gb~=1
         SOL=table2array(struct2table(solve(gb)));
         M=zeros(1,n);
         for l=1:size(SOL,1)
             for i=1:n
                 M(1,i)=SOL(l,i);
             end
             L{l}=M;
         end
    end
end

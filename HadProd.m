function H = HadProd(P,Q,L)
% HadProd  Compute the Hadamard quasigroup product of two arrays.
%    H = HadProd(P,Q,L) returns the Hadamard L-product of two nxn arrays P
%    and Q, with entries in {1,...,n} union {·}, with respect to a Latin 
%    square L of the same order n, with entries in {1,...,n}.  Empty cells 
%    those ones containing the symbol ·) within P and Q are represented by 
%    0.
%
%    Example:
%    HadProd([1 2 3; 3 1 2; 2 3 1], [1 2 3; 2 3 1; 3 1 2], [1 3 2; 3 2 1; 2 1 3])
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
    n=size(P,1);
    for i=1:n
        for j=1:n
            if P(i,j)*Q(i,j)==0
                H(i,j)=0;
            else
                H(i,j)=L(P(i,j),Q(i,j));
        end
    end
end

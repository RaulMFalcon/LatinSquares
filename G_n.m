function p=G_n(x,y,n)
% G_n  Auxiliar function to generate Latin squares.
%    p = G_n(x,y,n) returns a polynomial in x whose zeros constitute the 
%    symbols that can take two cells with distinct symbols within the Latin
%    square under consideration.
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

% Construct the polynomial.
    [r,p]=polynomialReduce(F_n(x,n)-F_n(y,n),x-y);
end

function [Q,R] = householder_QR(A)
% Calculates the QR factorisation of a given matrix using 
% Householder reflections.

[n,m] = size(A); % Dimensies van A
I = eye(n); Q = I; % Eenheidsmatrix (wordt orthogonale matrix)
R = A; % Wordt rechterbovendriehoeksmatrix

for col = 1:min(n-1,m)
        
    % Berekening van Householdermatrix
    v = R(col:end,col);
    v(1) = v(1) + (sign(v(1))+(v(1)==0)) * norm(v); % sign(0) geeft anders 0
    H = I(col:end,col:end) - (2*(v*v'))/(v'*v);
    
    % Nullen creëren
    R(col:n,col:m) = H * R(col:n,col:m);
    % R(col+1:end,col) = 0;
    Q(1:n,col:n) = Q(1:n,col:n) * H';

end
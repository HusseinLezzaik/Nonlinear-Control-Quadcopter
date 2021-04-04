function [X,Y,Z] = rotate_and_shift(X1, Y1, Z1, rot_mat, x, y, z)

    X_t = X1 + x;
    Y_t = Y1 + y;
    Z_t = Z1 + z;
    
    X = zeros(size(X_t));
    Y = zeros(size(Y_t));
    Z = zeros(size(Z_t));
    
    for i=1:size(X,1)
        for j=1:size(X,2)
            temp = rot_mat*[X_t(i,j);Y_t(i,j);Z_t(i,j)];
            X(i,j) = temp(1);
            Y(i,j) = temp(2);
            Z(i,j) = temp(3);
        end
    end
    
end
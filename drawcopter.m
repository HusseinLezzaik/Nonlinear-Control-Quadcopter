function ss_ = drawcopter(ss, x, y, z, psi, phi, theta,pausing)

    psi_rot = [cos(psi) sin(psi) 0
                   -sin(psi)  cos(psi) 0
                       0        0     1];

    theta_rot = [ cos(theta)  0  -sin(theta)
                        0       1      0
                   sin(theta)  0  cos(theta)];

    phi_rot = [1      0          0
                   0   cos(phi)  sin(phi)
                   0   -sin(phi)   cos(phi)];
    
    total_rot = phi_rot*theta_rot*psi_rot;

    [X1,Y1,Z1] = sphere(10);
    [X2,Y2,Z2] = sphere(10);
    [X3,Y3,Z3] = sphere(10);
    [X4,Y4,Z4] = sphere(10);
    [X5,Y5,Z5] = sphere(10);

    X1 = X1/10;
    Y1 = Y1/10;
    Z1 = Z1/10;

    X2 = (X2 - 6)/20;
    Y2 = Y2/20;
    Z2 = Z2/20;

    X3 = (X3 + 6)/20;
    Y3 = Y3/20;
    Z3 = Z3/20;

    X4 = X4/20;
    Y4 = (Y4 - 6)/20;
    Z4 = Z4/20;

    X5 = X5/20;
    Y5 = (Y5 + 6)/20;
    Z5 = Z5/20;


    [XX1,YY1,ZZ1] = cylinder(0.01,3);
    [XX2,YY2,ZZ2] = cylinder(0.01,3);
    [XX3,YY3,ZZ3] = cylinder(0.01,3);
    [XX4,YY4,ZZ4] = cylinder(0.01,3);
    
    temp = XX1;
    XX1 = 0.3*ZZ1;
    ZZ1 = temp;
    
    temp = XX2;
    XX2 = -0.3*ZZ2;
    ZZ2 = temp;
    
    temp = YY3;
    YY3 = 0.3*ZZ3;
    ZZ3 = temp;
    
    temp = YY4;
    YY4 = -0.3*ZZ4;
    ZZ4 = temp;

    [X1,Y1,Z1] = rotate_and_shift(X1, Y1, Z1, total_rot, x, y, z);
    [X2,Y2,Z2] = rotate_and_shift(X2, Y2, Z2, total_rot, x, y, z);
    [X3,Y3,Z3] = rotate_and_shift(X3, Y3, Z3, total_rot, x, y, z);
    [X4,Y4,Z4] = rotate_and_shift(X4, Y4, Z4, total_rot, x, y, z);
    [X5,Y5,Z5] = rotate_and_shift(X5, Y5, Z5, total_rot, x, y, z);
    
    [XX1,YY1,ZZ1] = rotate_and_shift(XX1, YY1, ZZ1, total_rot, x, y, z);
    [XX2,YY2,ZZ2] = rotate_and_shift(XX2, YY2, ZZ2, total_rot, x, y, z);
    [XX3,YY3,ZZ3] = rotate_and_shift(XX3, YY3, ZZ3, total_rot, x, y, z);
    [XX4,YY4,ZZ4] = rotate_and_shift(XX4, YY4, ZZ4, total_rot, x, y, z);

    
    s1 = surf(ss, X1,Y1,Z1);
    hold(ss,"on")
    s2 = surf(ss, X2,Y2,Z2);
    s3 = surf(ss, X3,Y3,Z3);
    s4 = surf(ss, X4,Y4,Z4);
    s5 = surf(ss, X5,Y5,Z5);

    s6 = surf(ss, XX1,YY1,ZZ1);
    s7 = surf(ss, XX2,YY2,ZZ2);
    s8 = surf(ss, XX3,YY3,ZZ3);
    s9 = surf(ss, XX4,YY4,ZZ4);
    
    
    if pausing
        pause(0.001)
    end
    
    set(s1,"visible", "off");
    set(s2,"visible", "off");
    set(s3,"visible", "off");
    set(s4,"visible", "off");
    set(s5,"visible", "off");
    set(s6,"visible", "off");
    set(s7,"visible", "off");
    set(s8,"visible", "off");
    set(s9,"visible", "off");
    
    ss_ = ss;
    
end


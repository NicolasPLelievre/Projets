n_inputs = 10;
n_outputs = 1;

function out=G(X)
    // Variables du problème
    E1 = X(:,1);
    E2 = X(:,2);
    A1 = X(:,3);
    A2 = X(:,4);
    P1 = X(:,5);
    P2 = X(:,6);
    P3 = X(:,7);
    P4 = X(:,8);
    P5 = X(:,9);
    P6 = X(:,10);
    
    
    //flèche maximale au point 7 en mètre
    Vmax = 0.12;
    // Longueur des barres
    L1 = 4; // barres horizontales
    L2= sqrt(8); // barres obliques
    
    // Matrice des barres de la structure
    T = zeros(23,2);
    i = 1;
    j = 1;
    while 1
        if j < 23
            T(j:j+1,:) = [i,i+1; i,i+2];
            j = j+2;
            i = i+1;
        else
            T(j,:) = [i,i+1];
            break
        end
    end
    
    // Matrice de rigidité élémentaire pour chacune des barres
    ca1 = 2/L2;
    sa1 = 2/L2;
    ca2 = ca1;
    sa2 = sa1;
    for k=1:length(E1)
        Kepair = E1(k)*A1(k)/L1*[1, 0, -1, 0; 0, 0, 0, 0; -1, 0, 1, 0; 0,0,0,0];
        Keimpair1 = E2(k)*A2(k)/L2*[ca1^2, ca1*sa1, -ca1^2, -ca1*sa1
                              ca1*sa1, sa1^2, -ca1*sa1, -sa1^2
                              -ca1^2, -ca1*sa1, ca1^2, ca1*sa1
                              -ca1*sa1,-sa1^2,ca1*sa1,sa1^2];
        Keimpair2 = E2(k)*A2(k)/L2*[ca2^2, -ca2*sa2, -ca2^2, ca2*sa2
                              -ca2*sa2, sa2^2, ca2*sa2, -sa2^2
                              -ca2^2, ca2*sa2, ca2^2, -ca2*sa2
                              ca2*sa2,-sa2^2,-ca2*sa2,sa2^2];
                              
        // Montage de la matrice de rigidité globale
        K = zeros(26,26);
        for i = 1:size(T,1)
            if modulo(T(i,1)+T(i,2), 2)==1
                if modulo(T(i,1),2)==1
                    K(2*T(i,1)-1:2*T(i,1),2*T(i,1)-1:2*T(i,1)) = K(2*T(i,1)-1:2*T(i,1),2*T(i,1)-1:2*T(i,1)) + Keimpair1(1:2,1:2);
                    K(2*T(i,1)-1:2*T(i,1),2*T(i,2)-1:2*T(i,2)) = K(2*T(i,1)-1:2*T(i,1),2*T(i,2)-1:2*T(i,2)) + Keimpair1(1:2,3:4);
                    K(2*T(i,2)-1:2*T(i,2),2*T(i,1)-1:2*T(i,1)) = K(2*T(i,2)-1:2*T(i,2),2*T(i,1)-1:2*T(i,1)) + Keimpair1(3:4,1:2);
                    K(2*T(i,2)-1:2*T(i,2),2*T(i,2)-1:2*T(i,2)) = K(2*T(i,2)-1:2*T(i,2),2*T(i,2)-1:2*T(i,2)) + Keimpair1(3:4,3:4);
                else
                    K(2*T(i,1)-1:2*T(i,1),2*T(i,1)-1:2*T(i,1)) = K(2*T(i,1)-1:2*T(i,1),2*T(i,1)-1:2*T(i,1)) + Keimpair2(1:2,1:2);
                    K(2*T(i,1)-1:2*T(i,1),2*T(i,2)-1:2*T(i,2)) = K(2*T(i,1)-1:2*T(i,1),2*T(i,2)-1:2*T(i,2)) + Keimpair2(1:2,3:4);
                    K(2*T(i,2)-1:2*T(i,2),2*T(i,1)-1:2*T(i,1)) = K(2*T(i,2)-1:2*T(i,2),2*T(i,1)-1:2*T(i,1)) + Keimpair2(3:4,1:2);
                    K(2*T(i,2)-1:2*T(i,2),2*T(i,2)-1:2*T(i,2)) = K(2*T(i,2)-1:2*T(i,2),2*T(i,2)-1:2*T(i,2)) + Keimpair2(3:4,3:4);
                end
            else
                K(2*T(i,1)-1:2*T(i,1),2*T(i,1)-1:2*T(i,1)) = K(2*T(i,1)-1:2*T(i,1),2*T(i,1)-1:2*T(i,1)) + Kepair(1:2,1:2);
                K(2*T(i,1)-1:2*T(i,1),2*T(i,2)-1:2*T(i,2)) = K(2*T(i,1)-1:2*T(i,1),2*T(i,2)-1:2*T(i,2)) + Kepair(1:2,3:4);
                K(2*T(i,2)-1:2*T(i,2),2*T(i,1)-1:2*T(i,1)) = K(2*T(i,2)-1:2*T(i,2),2*T(i,1)-1:2*T(i,1)) + Kepair(3:4,1:2);
                K(2*T(i,2)-1:2*T(i,2),2*T(i,2)-1:2*T(i,2)) = K(2*T(i,2)-1:2*T(i,2),2*T(i,2)-1:2*T(i,2)) + Kepair(3:4,3:4);
            end
        end
        
        // Vecteur des efforts
        F = [0;(P1(k)+P2(k)+P3(k)+P4(k)+P5(k)+P6(k))/2;0; P1(k);0; 0;0; P2(k);0; 0;0; P3(k);0; 0;0; P4(k);0; 0;0; P5(k);0; 0;0; P6(k);0; (P1(k)+P2(k)+P3(k)+P4(k)+P5(k)+P6(k))/2];
        // Conditions aux limites
        // Point 1 sur x
        K(1,:) = 0;
        K(:,1) = 0;
        K(1,1) = 1;
        // Point 1 sur y
        K(2,:) = 0;
        K(:,2) = 0;
        K(2,2) = 1;
        // Point 13 sur y
        K(26,:) = 0;
        K(:,26) = 0;
        K(26,26) = 1;
        
        // Résolution du problème
        Utot = K\F
        U14 = Utot(14);
        U(k) = U14;
    end
    // Comparaison du déplacement du point central par rapport à Vmax
    out = Vmax - abs(U);
endfunction

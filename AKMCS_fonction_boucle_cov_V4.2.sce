function [pf, cpf, time, nappel, nMC_used, nboucle, najout] = loopingnMC(ncluster)


    timer();
    stacksize('max');
    exec('G4variablesSphere.sce');
//    exec('G7variablesFlambage.sce');
//    exec('G10variablesTruss.sce');
    exec('ChargeFonctions.sce');
    
    ////////////////////////////////////////////////////////////////////////////////
    // Data
    seed = 21;
    nMC = 3e6;
    nMCstep = 1e4;
    Umin = 0.01;
    cmin = 0.02;
    nini = 11;
    ////////////////////////////////////////////////////////////////////////////////
    // 4 variables
    Seuil = 0;
    Loi = [2; 2; 2; 2];
    Moy = [300; 130; 100; 50];
    Stdev = [20; 8; 5; 2.5];
    PtTirage = [0; 0; 0;0];
    
    ////////////////////////////////////////////////////////////////////////////////
    // 7 variables
//    Seuil = -1;
//    Loi = [1; 1; 1; 1; 1; 1; 1];
//    Moy = [24; 600; 156; 10; 120; 24; 200000];
//    Stdev = [0.72; 2; 4.68; 0.3; 3.6; 0.72; 400];
//    PtTirage = [0; 0; 0; 0; 0; 0; 0];
    
    ////////////////////////////////////////////////////////////////////////////////
    // 10 variables
//    Seuil = 0;
//    Loi = [2; 2; 2; 2; 5; 5; 5; 5; 5; 5];
//    Moy = [2.1e11; 2.1e11; 2e-3; 1e-3; 5e4; 5e4; 5e4; 5e4; 5e4; 5e4];
//    Stdev = [2.1e10; 2.1e10; 2e-4; 1e-4; 7.5e3; 7.5e3; 7.5e3; 7.5e3; 7.5e3; 7.5e3];
//    PtTirage = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
    
    nva = length(Loi);
    Ustar = PtTirage;
    
    ////////////////////////////////////////////////////////////////////////////////
    // AKMCS
    
    // Calcul du premier plan d'expérience
    fdoe = grand(nini, nva, 'unf', -3.5, 3.5);
    first_doe = fdoe
    v = variance(first_doe);
    for k = 1:50
        fdoe = grand(nini, nva, 'unf', -3.5, 3.5);
        v2 = variance(fdoe);
        if v2 > v
            first_doe= fdoe
            v = v2;
        end
    end
    
    // Definition de l'état du générateur aléatoire
    grand("setsd", seed);
    
    //options AK-MCS
    if exists('theta0') == 0 then
        theta0 = 0.01 * ones(1, nva);
    end
    lower_b = 1e-10 * ones(1, nva);
    upper_b = 1e2 * ones(1, nva);
    
    p = [1:nini]';
    najout = 0;// Nombre de calculs ajouté à chaque itération, n'est pas forcément égal à ncluster notament à la dernière itération, surtout utile  l'affichage des courbes de CV
    niterj = [];// Nombre de calcul par step de nMC
    minPINF = 1;
    maxPSUP = 0;
    np = 1e4; // Nombre limite de paquet à prédire pour la fonction MultiPredictpour
    i = 0;
    j=0;
    first_base = 0;
    ii = 0;
    
    Base = grand(nva, nMCstep, 'nor', 0, 1)';
    Base = Base + ones(nMCstep,1) * Ustar';
    BBASE =[];
    while 1
        i = i + 1;
        critOK = 0; // Booléen permetant de savoir si le critère sur U est respécté.
        if i == 1 // Si c'est le premier plan d'expériences
            DOE_u = first_doe; // Premier plan d'exp en U par choix aléatoire
            DOE_y = Perf(DOE_u, Seuil); // Calcul de la fonction de performance
        else //si ce n'est pas le premier plan d'experiences
            DOE_u = [DOE_u; new_u]; // On ajoute le point
            DOE_y = [DOE_y; new_y];
        end
        if exists('Gmodel') == 1 then
            theta0 = Gmodel.theta;
        end
        
        // gestion pour la mauvaise optimisation des paramètres theta
        Gmodel = BuildKrigingModel(DOE_u, DOE_y, theta0, lower_b, upper_b);
        theta0 = Gmodel.theta;
        [mu_G, var_G] = MultiPredict(Base, np, Gmodel);//Predict the points of the population x with DACE Kriging toolbox
        
        
        // Calcul de la probabilite de defaillance 
        n_negatif(j+1) = length(mu_G(mu_G<=0))
        Pf_AK(j+1)= sum(n_negatif)/(nMCstep*(j+1)); //Calculate the failure probability prediction
        
        // Si la probabilité AK est nulle
        if Pf_AK($) == 0 | Pf_AK($) == 1
            cov_AK(i) = %inf;
            if j == (nMC/nMCstep)-1
                disp("ATTENTION : Le nombre total de points à classer est trop faible")
                disp("ATTENTION : Aucun point trouvé dans le domaine de défaillance")
                break
            end
            if first_base == 0
                ii = ii + 1;
                new_u = Base(nini+ii,:);
                new_y = Perf(new_u, Seuil);
                
                if ii == 5
                    first_base = 1;
                end
            else
                j = j+1;
                new_u = [];
                new_y = [];
                Base = grand(nva, nMCstep, 'nor', 0, 1)';
                Base = Base + ones(nMCstep,1) * Ustar';
            end
            najout(i) = size(new_u, 1);

        else
            first_base = 1;
            // Coefficient de variation de la probabilité
            cov_AK(i)=sqrt((1-Pf_AK($))/((nMCstep*(j+1))*Pf_AK($))); //Coefficient de variation

            NCALCG_i0 = size(DOE_u, 1);
            NCALCG = size(DOE_u, 1);
            var_G(var_G <= %eps) = %eps;

            // Calcul de U
            U = abs(mu_G)./ sqrt(var_G); // Evaluation du critere d'apprentissage
            [mU pos] = min(U); //Recherche du meilleur point
            // Intervalles de confiance
            Pinf(i)=length(mu_G(mu_G<=0&U>2))/(nMCstep);
            Psup(i)=1-length(mu_G(mu_G>=0&U>2))/(nMCstep);
            if Pinf(i)==0 | Pinf(i)==1
                cov_inf(i)=1;
                PINF(i)=0;
            else
                cov_inf(i)=sqrt((1-Pinf(i))/(nMCstep*Pinf(i)));
                PINF(i)=Pinf(i)*(1-1.96*cov_inf(i));
            end
            if Psup(i)==0 | Psup(i)==1
                cov_sup(i)=1;
                PSUP(i)=1;
            else
                cov_sup(i)=sqrt((1-Psup(i))/(nMCstep*Psup(i)));
                PSUP(i)=Psup(i)*(1+1.96*cov_sup(i));
            end
        
            //Affichage de la courbe de CV
            if i > 1
                minPINF = min(minPINF, PINF(i));
                maxPSUP = max(maxPSUP, PSUP(i));
            end
        
            minU(i)= mU; // Evolution de min de U
            
            // test critère d'arrêt
    
            // critere d'arret
            if Pinf($)==0
                crit=1;
            else
                crit=abs(Pinf($)-Psup($))/Pinf($);
            end
            if crit < Umin
                critOK = 1;
                if cov_AK($) < cmin
                    niterj(j+1) = [sum(najout)-sum(niterj(1:j))];
                    disp("FIN DU CALCUL : CONVERGENCE")
                    break;  // convergence
                elseif j == (nMC/nMCstep)-1
                    niterj(j+1) = [sum(najout)-sum(niterj(1:j))];
                    disp("ATTENTION : Le CoV de Pf souhaité n''a pas été atteint. Le nombre de point total est trop faible.");
                    break;  // convergence
                else
                    if j == 0
                        niterj(j+1) = [sum(najout)];
                        j = j+1;
                        Base = grand(nva, nMCstep, 'nor', 0, 1)';
                        Base = Base + ones(nMCstep,1) * Ustar';
                    else
                        niterj(j+1) = [sum(najout)-sum(niterj(1:j))];
                        if modulo(j,500) == 0 then
                            disp(nMCstep*(j+1)*100/nMC);
                            disp( "% des points ont été classés");
                        end
                        j = j+1;
                        Base = grand(nva, nMCstep, 'nor', 0, 1)';
                        Base = Base + ones(nMCstep,1) * Ustar';
                    end
                end
            elseif crit == 1 & size(U(U>2)) == size(U)
                if j == (nMC/nMCstep)-1
                    disp("ATTENTION : Le nombre total de points à classer est trop faible")
                    break;
                end
                critOK = 1;
                niterj(j+1) = [sum(najout)-sum(niterj(1:j))];
                j = j+1;
                Base = grand(nva, nMCstep, 'nor', 0, 1)';
                Base = Base + ones(nMCstep,1) * Ustar';
            end
            // Si le critère d'arret n'est pas atteint, nouveaux points à ajouter
            if critOK == 0 then 
                if ncluster == 1
                    new_u = Base(pos, :);
                else
                    if size(Base(U < 2, :), 1) > ncluster
                        [new_u, place] = FindCluster(Base(U < 2, :), ncluster);
                        disp("hey hey je suis là");
                    else
                        new_u = Base(U < 2, :);
                    end
                end
                
                
                // Calcul de la fonction de performance
                new_y = Perf(new_u, Seuil);
            else
                new_u = [];
                new_y = [];
            end
            najout(i) = size(new_u, 1);
        end
        
    end
    pf = Pf_AK($);
    cpf = cov_AK($);
    time = timer();
    nappel = nini + sum(najout);
    nMC_used = (j+1) * nMCstep;
    nboucle = i;
endfunction

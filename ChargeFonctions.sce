// Charge toutes les fonctions de l'interface graphique

////////////////////////////////////////////////////////////////////////////
// Fonctions de l'IHM
////////////////////////////////////////////////////////////////////////////

function charge_pdf_aide(dir, version)
    tmp = strsplit(version, '.')
    i1 = tmp(1,1)
    i2 = tmp(2,1)
    winopen(dir + '\Aide\Aide IHM AKMCSIS V' + i1 + '_' + i2 + '.pdf')
endfunction

function verification_version_module()
    if getversion()<>'scilab-5.4.1' then
       disp('ATTENTION: le programme est écrit pour Scilab 5.4.1. Certaines fonctions peuvent ne pas fonctionner correctement.', 'warning', ["OK"], 'modal')
    end
     if ~atomsIsInstalled('dace_scilab') then
        messagebox('ATTENTION: le module dace_scilab n''est pas installé. Veuillez installer le module.', 'error', ["OK"], 'modal')
    end
    if ~atomsIsLoaded('dace_scilab') then
        atomsLoad('dace_scilab')
    end
endfunction

function init_graphe_cv()
    xsetech([0.5,0,0.5,0.5])
    xlabel('Nombre de calculs')
    ylabel('Probabilité')
    set(gca(), 'data_bounds', matrix([0, 1, 0, 1], 1, -1))
    xgrid()
endfunction

function sauveEtude()
    savefic=uiputfile(["*.dat"],pwd(),'Sauvegarde de l''étude au format .dat');
    if isempty(savefic) then
        abort;
    end
    a1 = ficf.string
    a2 = Seuils.string
    a3 = sprintf("%d", n_inputs)
    a4 = LoiVar.string
    a5 = MoyVar.string
    a6 = StdevVar.string
    a7 = nMCs.string
    a8 = seeds.string
    a9 = Umins.string
    a10 = ninis.string
    a11 = Pf.string
    a12 = cv.string
    a13 = NCALCGSTR.string
    a14=Pf_AK
    a15=DOE_u
    a16=Pinf
    a17=Psup
    a18=PINF
    a19=PSUP
    a20=DOE_x;
    a21=DOE_y;
    a22=nclusters.string;
    a30=PtTirage.string;
    a31=BetaFORM.string;
    a32=PfFORM.string;
    a33=niterFORM.string;
    if ~exists('theta0') then
        a34 = 1e-2*ones(1,n_inputs)
    else
        a34=theta0;
    end
    if strcmp(part(savefic,length(savefic)-3:length(savefic)),'.dat') == 0 then
        save(savefic, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a30, a31, a32, a33, a34)
    else
        save(savefic + '.dat', a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a30, a31, a32, a33, a34)
    end
    
endfunction

function exportEtude()
    savefic=uiputfile(["*.sco"],pwd(),'Export de l''étude au format .sco');
    if isempty(savefic) then
        abort;
    end

    a10 = ninis.string
    if ninis.string == '' | ninis.string == '--' then
        savefic = ''
        texte = 'Vous devez spécifier une taille pour le premier plan d''expérience !!'
    else
        if evstr(ninis.string) < n_inputs then
            savefic = ''
            texte = 'ATTENTION : la taille du premier plan d''experiences doit être supérieure au nombre de variables'
        end
    end

    a9 = Umins.string
    if Umins.string == '' then
        savefic = ''
        texte = 'Vous devez spécifier un critère d''arrêt !!'
    end

    a8 = seeds.string
    if seeds.string == '' then
        savefic = ''
        texte = 'Vous devez spécifier un état pour le générateur aléatoire !!'
    end

    a7 = nMCs.string
    if nMCs.string == '' then
        savefic = '';
        texte = 'Vous devez spécifier une taille de population à classer !!';
    end

    a6 = StdevVar.string
    for i=1:n_inputs
        if evstr(StdevVar.string(1,i)) <= 0. then
            savefic = ''
            texte = sprintf("L''écarts type de la variables aléatoire %d est négatif !", i)
        end
    end

    a4 = LoiVar.string
    for i=1:n_inputs
        if LoiVar.string(1,i) <> 'D' & LoiVar.string(1,i) <> 'N' & LoiVar.string(1,i) <> 'LN' & LoiVar.string(1,i) <> 'U' & LoiVar.string(1,i) <> 'W' then
            savefic = ''
            texte = 'Les lois des variables aléatoires sont mal renseignées !!'
        end
    end

    a3 = sprintf("%d", n_inputs)

    a2 = Seuils.string
    if Seuils.string == '' then
        savefic = ''
        texte = 'Vous devez spécifier une valeur de seuil !!'
    end

    a1 = ficf.string
    if ficf.string == '' then
        savefic = ''
        texte = 'Vous devez spécifier une fonction de performance !!'
    end
    
    a22 = nclusters.string
    if nclusters.string == '' then
        savefic = ''
        texte = 'Vous devez spécifier un nb de points à enrichir à chaque itération !!'
    end

    if (savefic <> '') & strcmp(part(savefic,length(savefic)-3:length(savefic)),'.sco')==0 then
        save(savefic, a1, a2, a3, a4, a6, a7, a8, a9, a10, a22)
    elseif (savefic <> '') & strcmp(part(savefic,length(savefic)-3:length(savefic)),'.sco')<>0
        save(savefic + '.sco', a1, a2, a3, a4, a6, a7, a8, a9, a10, a22)
    else
        messagebox(texte)
    end
endfunction

function [Pf_AK,DOE_u,DOE_x,DOE_y,Pinf,Psup,PINF,PSUP,n_inputs,Ustar,G,Loi,Moy,Stdev,Seuil,theta0]=chargeEtude()
    fi=uigetfile(['*.dat'],pwd(),'Ouverture d''un fichier d''étude au format .dat');
    if fi == "" then
        abort;
    end
    try
        load(fi, "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10", "a11", "a12", "a13", "a14", "a15", "a16", "a17", "a18", "a19", "a20", "a21", "a22", "a30", "a31", "a32", "a33", "a34");
    catch
        messagebox("Fichier corrompu!", "Erreur: chargeEtude", "error");
        abort;
    end
    reset_Var()
    if a1 <> '' then
        fici = fileinfo(a1)
        sfici = size(fici)
        if sfici(1,1) <> 0 then
            set(ficf,'string', a1)
            exec(ficf.string, -1)
        else
            messagebox("Le fichier définissant la fonction de performance n''existe pas. Vous devez spécifier un nouveau fichier avant de lancer un calcul.", ...
                       "Avertissement: chargeEtude", "warning", ["OK"], "modal");
            [n_inputs, G] = ChargeG();
        end
    end
    set(nvariables,'string', sprintf("%d", n_inputs))

    set(Seuils,'string', a2)
    adapt_nb_var(n_inputs)
    for i=1:n_inputs
        set(LoiVar(i,1),'string', a4(1,i))
    end
    check_var_det(n_inputs)
    for i=1:n_inputs
        set(MoyVar(i,1),'string', a5(1,i))
        set(StdevVar(i,1),'string', a6(1,i))
        set(PtTirage(i,1),'string', a30(1,i))
    end
    [Loi,Moy,Stdev]=LoiVariables()
    Seuil=evstr(Seuils.string);
    
    set(nMCs,'string', a7)
    set(seeds,'string', a8)
    set(Umins,'string', a9)
    set(ninis,'string', a10)
    set(Pf,'string', a11)
    set(cv,'string', a12)
    set(NCALCGSTR,'string', a13)
    
    Ustar = evstr(PtTirage.string)
    set(BetaFORM,'string', a31)
    set(PfFORM,'string', a32)
    set(niterFORM,'string', a33)

    Pf_AK=a14;
    DOE_u=a15;
    DOE_x=a20;
    DOE_y=a21;
    Pinf=a16;
    Psup=a17;
    PINF=a18;
    PSUP=a19;
    theta0=a34;
    
    set(nclusters,'string', a22)
    ncluster=evstr(a22);

    if DOE_u <> [] then
        nini = ninis.string
        a=get("current_axes")
        delete(a)
        xsetech([0.5,0,0.5,0.5])
        xlabel('Nombre de calculs')
        ylabel('Probabilité')
        minPINF = min(PINF)
        maxPSUP = max(PSUP)
        taillePf = size(Pf_AK)
        tailleDOE = size(DOE_u)
        tailleplot = tailleDOE(1) + (1 - taillePf(1)) * ncluster
        X = [tailleplot:ncluster:tailleDOE(1)]
        plot(X,Pf_AK, 'ko-',...
             X,Pinf, 'k--',...
             X,PINF, 'r--',...
             X,Psup, 'k--',...
             X,PSUP, 'r--')
        xgrid()
        Umin=evstr(Umins.string);
        legend(['Probabilité'; sprintf('IC Krigeage (Umin = 2)'); 'IC Monte-Carlo (95%)'], ...
               pos=1)
    end
endfunction

function [DOE_x, DOE_y]=ChargePLEX()
    if nvariables.string == '--' then
        messagebox('Veuillez commencer par sélectionner un modèle!', "Erreur: ChargePLEX", "error");
        abort;
    end
    fi=uigetfile(['*.csv'],pwd(),'Import d''un plan d''expériences au format .csv');
    if fi == "" then
        abort;
    end
    table = fscanfMat(fi);
    n_inputs = evstr(nvariables.string);
    if size(table, 2) <> n_inputs+1 then
        messagebox('La taille de la table contenue dans le PLEX sélectionné n''est pas compatible avec le modèle sélectionné!', "Erreur: ChargePLEX", "error");
        abort;
    end
    DOE_x = table(:, 1:$-1);
    DOE_y = table(:, $);
endfunction

function SauvePLEX(DOE_x, DOE_y)
    if isempty(DOE_x) then
        messagebox("Il n''y a aucun PLEX à sauver!", "Erreur: SauvePLEX", "error");
        abort;
    end
    savefic=uiputfile(["*.csv"],pwd(),'Export du plan d''expériences au format .csv');
    if isempty(savefic) then
        abort;
    end
    if strcmp(part(savefic,length(savefic)-3:length(savefic)),'.csv')==0 then
        fprintfMat(savefic, [DOE_x, DOE_y], "%.12g,")
    else
        fprintfMat(savefic + '.csv', [DOE_x, DOE_y], "%.12g,")
    end
endfunction

function apropo()
    texte=[titre_programme,version_scilab,'','Programme développé par IFMA/2MATECH/PHIMECA pour SNECMA ',...
          '2015', '','Version ' + version];
    messagebox(texte,'A propos de AK-MCSIS V' + version,'info','modal');
endfunction

function [n_inputs, G]=ChargeG()
    fi=uigetfile(['*.sce'],pwd(),'Choisir un fichier .sce');
    if fi == "" then
        abort;
    end
    set(ficf,'string', fi)
    exec(ficf.string, -1)
    if ~exists("n_inputs") | ~exists("n_outputs") | ~exists("G") then
        messagebox("Le fichier choisi n''est pas complet. Il manque au moins une variable parmi ''n_inputs'', ''n_outputs'', ou bien il n''implémente aucune fonction ''G''.", "Erreur: ChargeG", "error")
        abort;
    elseif typeof(G) <> "function" then
        messagebox("''G'' n''est pas une fonction dans le fichier choisi.", "Erreur: ChargeG", "error");
        abort;
    end

    if n_inputs > 12 | n_inputs < 1 then
        messagebox(sprintf('AK-MCS ne fonctionne que pour les fonctions à moins de 12 variables. Ici n_inputs = %d.', n_inputs), "Error: ChargeG", "error")
        abort;
    end
    if n_outputs <> 1 then
        messagebox(sprintf('AK-MCS ne fonctionne que pour les fonctions à valeurs scalaires (n_outputs = 1). Ici n_outputs = %d.', n_outputs), "Error: ChargeG", "error")
        abort;
    end
    reset_Var()
    adapt_nb_var(n_inputs)
    set(nvariables,'string',sprintf("%d", n_inputs))
endfunction

function openscinotes(fi)
    if fi == '' then
       messagebox('Vous devez choisir la fonction de performance avant de l''éditer.')
       abort;
    else
      scinotes(fi)
    end
endfunction

function adapt_nb_var(nva)
    if nva>12
        disp('ATTENTION : Le nombre de variables doit etre inferieur à 12.')
        abort;
    end
    for i=1:nva
       set(NomVar(i),'visible','on');
       set(LoiVar(i),'visible','on');
       set(MoyVar(i),'visible','on');
       set(StdevVar(i),'visible','on');
       set(PtTirage(i),'visible','on');
       set(PtTirage(i),'string','0');
   end
   if nva<12
      for i=nva+1:12
        set(NomVar(i),'visible','off');
        set(LoiVar(i),'visible','off');
        set(MoyVar(i),'visible','off');
        set(StdevVar(i),'visible','off')
        set(PtTirage(i),'visible','off');
     end
  end
endfunction

function check_var_det(nva)
    for i=1:nva
        if LoiVar(i).string == 'D' then
            set(StdevVar(i),'visible','off')
            set(StdevVar(i,1),'string', '%inf')
        else
            set(StdevVar(i),'visible','on')
        end
    end
endfunction

function [k] = indiceDistribution(loi, i)
    if loi == 'D' then
        k = 0
    elseif loi == 'N' then
        k = 1
    elseif loi == 'LN' then
        k = 2
    elseif loi == 'U' then
        k = 3
    elseif loi == 'W' then
        k = 4
    else
        k = -1
    end
endfunction

function [mu_G,var_G]=MultiPredict(Base,np,Gmodel)
        BaseTmp=Base;
        mu_G=[];
        var_G=[];
        while 1
            if size(BaseTmp,1)>=np
               [m,v]=predictor(BaseTmp([1:np],:),Gmodel)
                BaseTmp(1:np,:)=[];
                mu_G=[mu_G;m];
                var_G=[var_G;v] 
            elseif size(BaseTmp,1)>0 & size(BaseTmp,1)<np
                [m,v]=predictor(BaseTmp(1:size(BaseTmp,1),:),Gmodel)
                mu_G=[mu_G;m];
                var_G=[var_G;v]
                break;
            else
                break
            end 
        end
endfunction

function [nva,Seuil,Loi,Moy,Stdev,nMC,Umin,nini,ncluster,Pf_AK,DOE_u,DOE_x,DOE_y,Pinf,Psup,PINF,PSUP,theta0]=correspondance()
    
    if ninis.string == '--' then
        messagebox('La taille du premier plan d''experiences n''est pas définie !');
        abort
    end

//    if eval(ninis.string) < 5 * n_inputs then
//        answer = messagebox("La taille du premier plan d''experiences est faible : un minimum de 5 fois le nombre de variables est recommandé. Souhaitez-vous poursuivre quand même ?", "modal", "warning", ["Oui" "Non"]);
//        if answer == 2 then
//            abort;
//        end
//    end

    //br = 0
    set(boutonstart,'fontSize',10);
    set(boutonstart,'FontWeight','light');
    set(boutonarret,'fontSize',16);
    set(boutonarret,'FontWeight','bold');

    // Nom du fichier 
    fic = ficf.string;
    exec(fic,-1);

    // Nombre de variables
    nva = n_inputs;

    // Valeur de seuil
    Seuil=evstr(Seuils.string);

    // Lois des variables
    [Loi,Moy,Stdev]=LoiVariables();
    i = find(Loi==-1);
    if i<>[] then
        messagebox('La loi de la variable ' + string(i) + ' est mal definie !');
        abort;
    end
    
    Ustar=evstr(PtTirage.string)';
    

    // Plan d'expériences initial
    if ~isempty(DOE_x) then
        if size(DOE_x, 2) <> n_inputs then
            messagebox(sprintf('Le plan d''expériences initial chargé en mémoire ne correspond pas à la fonction étudiée (n_inputs de ''DOE_x'' = %d).', size(DOE_x, 2)), "Error: n_inputs", "error");
            set(boutonstart,'fontSize',16);
            set(boutonstart,'FontWeight','bold');
            set(boutonarret,'fontSize',10);
            set(boutonarret,'FontWeight','light');
            xinfo('Prêt');
            abort;
        end
        DOE_u = tiso(DOE_x);
        new_u = DOE_u;
        new_y = DOE_y;
        DOE_u = [];
        DOE_x = [];
        DOE_y = [];
        i = 1;
    else
        i = 0;
    end

    // Etat du générateur élatoire
    seed = evstr(seeds.string);
    if seed < 0 | seed >= 2^32 then
        messagebox('Le germe du générateur aléatoire doit être un entier dans [0, 2^32[.', "Erreur: nMC", "error");
        abort;
    end

    // Taille de la population de MC
    nMC=evstr(nMCs.string);

    // Taille du plan d'expériences initial
    Umin=evstr(Umins.string);
    nini=evstr(ninis.string);
    if (nini < nva | nini > nMC) & i == 0 then
       messagebox('La taille du premier plan d''expériences doit être supérieure au nombre de variables et inférieure à la taille de la population de Monte Carlo.', "Erreur: nini", "error");
       set(boutonstart,'fontSize',16);
       set(boutonstart,'FontWeight','bold');
       set(boutonarret,'fontSize',10);
       set(boutonarret,'FontWeight','light');
       xinfo('Prêt');
       abort;
    end
    
    // Nombre de points enrichis
    ncluster=evstr(nclusters.string);

    // Execution
    interface=1;
    exec('Fonctions\AKMCSexec.sce',-1);

    DOE_x = tisoinv(DOE_u);

    set(boutonstart,'fontSize',16);
    set(boutonstart,'FontWeight','bold');
    set(boutonarret,'fontSize',10);
    set(boutonarret,'FontWeight','light');
    xinfo('Prêt');

endfunction

function [BetaFORMvalue,PfFORMvalue,Ustarvalue]=correspondanceFORM()
    
    set(boutonstartFORM,'fontSize',10);
    set(boutonstartFORM,'FontWeight','light');
    set(boutonarretFORM,'fontSize',16);
    set(boutonarretFORM,'FontWeight','bold');

    // Nom du fichier 
    fic = ficf.string
    exec(fic,-1)

    // Nombre de variables
    nva = n_inputs;

    // Valeur de seuil
    Seuil=evstr(Seuils.string);

    // Lois des variables
    [Loi,Moy,Stdev]=LoiVariables()
    i = find(Loi==-1)
    if i<>[] then
        messagebox('La loi de la variable ' + string(i) + ' est mal definie !')
        abort;
    end
    
    //exécution
    interface=1;
    exec('Fonctions\FORMexec.sce',-1)
    
    set(boutonstartFORM,'fontSize',16);
    set(boutonstartFORM,'FontWeight','bold');
    set(boutonarretFORM,'fontSize',10);
    set(boutonarretFORM,'FontWeight','light');
    xinfo('Prêt');
    
endfunction


function recuperer_coordonnee_ustar()
    if Ustar<>[] then
        for i = 1:n_inputs
            set(PtTirage(i,1),'string', sprintf('%3.3f',Ustar(i)))
        end
    end
endfunction

function reset_PtTirage()
    for i = 1:n_inputs
        set(PtTirage(i,1),'string', '0')
    end
endfunction


function stop()
//    xsetech([0.5,0,0.5,0.5])
//    a=get("current_axes")
//    delete(a)
//    init_graphe_cv()
//    set(Pf,'string', '--')
//    set(cv,'string', '--')
//    set(NCALCGSTR,'string', '--')
    set(boutonstart,'fontSize',16);
    set(boutonstart,'FontWeight','bold');
    set(boutonarret,'fontSize',10);
    set(boutonarret,'FontWeight','light');
    xinfo('Prêt');
    abort;
endfunction

function stopFORM()
    set(boutonstartFORM,'fontSize',16);
    set(boutonstartFORM,'FontWeight','bold');
    set(boutonarretFORM,'fontSize',10);
    set(boutonarretFORM,'FontWeight','light');
    xinfo('Prêt');
    abort;
endfunction


function [Loi,Moy,Stdev]=LoiVariables()
    for i=1:n_inputs
        Loi(i)=indiceDistribution(LoiVar(i).string, i)
        Moy(i)=evstr(MoyVar(i).string);
        Stdev(i)=evstr(StdevVar(i).string);
    end
endfunction

function reset_PtTirage()
    for i = 1:n_inputs
        set(PtTirage(i,1),'string', '0')
    end
endfunction

function reset_Var()
    for i=1:12  //nb max de var
        set(LoiVar(i),'string','');
        set(MoyVar(i),'string','');
        set(StdevVar(i),'string','');
        set(PtTirage(i),'string','');
    end
endfunction

function [Loi,Moy,Stdev,Seuil,DOE_x,DOE_y,DOE_u]=reset()
    if ~exists('Loi') then
        Loi = zeros(n_inputs,1);
    end
    if ~exists('Moy') then
        Moy = zeros(n_inputs,1);
    end
    if ~exists('Stdev') then
        Stdev = zeros(n_inputs,1);
    end
    
    [LoiNew,MoyNew,StdevNew]=LoiVariables();
    SeuilNew = evstr(Seuils.string);
    
    if ~isequal(LoiNew,Loi) | ~isequal(MoyNew,Moy) | ~isequal(StdevNew,Stdev) | ~isequal(SeuilNew,Seuil) then
        set(Pf,'string','--');
        set(cv,'string','--');
        set(NCALCGSTR,'string','--');
        DOE_x=[];
        DOE_y=[];
        DOE_u=[];
        set(ninis,'string','--');
        set(BetaFORM,'string','--');
        set(PfFORM,'string','--');
        set(niterFORM,'string','--');
        init_graphe_cv()
        reset_PtTirage()
    end
    Loi = LoiNew;
    Moy = MoyNew;
    Stdev = StdevNew;
    Seuil = SeuilNew;
endfunction

////////////////////////////////////////////////////////////////////////////
// Fonctions pour le clustering
////////////////////////////////////////////////////////////////////////////
function C=Correspondance(P,X)
    // Regroupe les points par subsets en donnant le vecteur C des numeros des 
    // points P le plus proche de chaque point X
    dist=[];
    for i=1:size(P,1)
        dist=[dist sum((X-ones(size(X,1),1)*P(i,:)).^2,2)];
    end
    [a,pos]=gsort(dist,'c');    // pos est une matrice
    C=pos(:,size(pos,2));
endfunction

function [Xcluster,place]=FindCluster(X,n)
    // X : points à synthétsier
    // n : nombre de points clsuter souhaités
    // U est le vecteur des valeurs de U associées à la population X
    // Xclsuter : points issus de l'opération de clustering
    // Les points initiaux sont tirés aléatoirement dans la population à synthéstiser.
    j=1;
    nva=size(X,2);
    // points de depart
    t=grand(n,1,"uin",1,size(X,1));
    P=X(t,:);
    
    while 1
        C=Correspondance(P,X);
        for i=1:size(P,1)
            Subsets=X(C==i,:);
            if Subsets==[] then
                P(i,:)=-100*ones(1,nva);    // On enleve un point de depart
                Pnew(i,:)=-100*ones(1,nva);
            else
            Pnew(i,:)=mean(Subsets,1);
            end
        end
        if  (sum((P-Pnew).^2)<0.0000001) //| j==10
            Pini=Pnew;
            break;
        else
            P=Pnew;
            j=j+1
        end
    end
    Xcluster=Pini
    // On prend les points de la population les plus proches des points issue de l'opération de clsutering
    for i=1:size(Xcluster,1)
        dist=sum((X-ones(size(X,1),1)*Xcluster(i,:)).^2,2);
        [a,b]=min(dist);
        place(i)=b;
        Xcluster(i,:)=X(b,:);
        X(b,:)=[];  // On enlève le point pour être sur qu'il ne soit pas utiliser deux fois
    end
endfunction


////////////////////////////////////////////////////////////////////////////
// Fonctions de la methode AKMCS
////////////////////////////////////////////////////////////////////////////


// Tansformation isoprobabiliste
function X=tisoinv(U, Loi, Moy, Stdev)
    X = zeros(U)
    for i=1:nva
        // Loi déterministe X = moyenne
        if Loi(i) == 0 then
            X(:,i) = Moy(i)
        // Loi normale X = sigma * U + Mu
        elseif Loi(i) == 1 then
            X(:,i) = Stdev(i) * U(:,i) + Moy(i)
        // Loi Lognormale ln(X) = sigmaLog * U + MuLog
        elseif Loi(i) == 2 then
            mulog = log((Moy(i) ** 2)/(sqrt(Moy(i) ** 2 + Stdev(i) ** 2)))
            siglog = sqrt(log(1. + Stdev(i) ** 2 / Moy(i) ** 2))
            X(:,i) = exp(siglog * U(:,i) + mulog)
        // Loi Uniforme X = (b-a) * Phi(U) + a [avec a et b bornes inf et sup calculées à partir de Mu et Sigma]
        elseif Loi(i) == 3 then
            a = Moy(i) - sqrt(3) * Stdev(i)
            b = Moy(i) + sqrt(3) * Stdev(i)
            X(:,i) = (b-a) * cdfnor("PQ", U(:,i), zeros(U(:,i)), ones(U(:,i))) + a
        // Loi de Weibull X = f(U, lambdaW, kW) [avec lambdaW et kW paramètres internes calculés à partir de Mu et Sigma]
        elseif Loi(i) == 4 then
            ratio = 1.0 + (Stdev(i) / Moy(i)) ** 2.
            betaMin = 1.
            betaMax = 1.
            step = 0.5
            epsilon = 1.e-12

            if ratio > 2 then
                betaMin = betaMin - step
                step = step * 0.5
                t = exp(log(gamma(1.0 + 2.0 / betaMin)) - 2.0 * log(gamma(1.0 + 1.0 / betaMin)))
                while t < ratio
                    betaMin = betaMin - step
                    step = step * 0.5
                    t = exp(log(gamma(1.0 + 2.0 / betaMin)) - 2.0 * log(gamma(1.0 + 1.0 / betaMin)))
                end
                betaMax = betaMin + 2.0 * step

            else
                betaMax = betaMax + step
                step = step * 2.0
                t = exp(log(gamma(1.0 + 2.0 / betaMax)) - 2.0 * log(gamma(1.0 + 1.0 / betaMax)))
                while t >= ratio
                    betaMax = betaMax + step
                    step = step * 2.0
                    t = exp(log(gamma(1.0 + 2.0 / betaMax)) - 2.0 * log(gamma(1.0 + 1.0 / betaMax)))
                end
                betaMin = betaMax - 0.5 * step
            end

            while 1
                kW = 0.5 * (betaMin + betaMax)
                if (betaMax - betaMin <= epsilon * (1.0 + abs(betaMax + betaMin))) then
                    lambdaW = Moy(i) / gamma(1.0 + 1.0 / kW)
                    break
                end
                t = exp(log(gamma(1.0 + 2.0 / kW)) - 2.0 * log(gamma(1.0 + 1.0 / kW)))
                if (t < ratio) then
                    betaMax = kW
                else
                    betaMin = kW
                end
            end

            X(:,i) = lambdaW * (-1.0 * log(1.0 - cdfnor("PQ",U(:,i),zeros(U(:,i)), ones(U(:,i))))) ** (1./kW)
        // Loi de Gumbel
        elseif Loi(i) == 5 then
            cstgamma = 0.5772;
            Beta = sqrt(6)*Stdev(i)/%pi;
            Mu = Moy(i) - Beta*cstgamma;
            X(:,i) = Mu - Beta*log(-log(cdfnor("PQ",U(:,i), zeros(U(:,i)), ones(U(:,i)))));
            
        end
    end
endfunction

function U=tiso(X, Loi, Moy, Stdev)
    U = zeros(X)
    for i=1:nva
        // Loi déterministe U = 0
        if Loi(i) == 0 then
            U(:,i) = 0
        // Loi normale U = (X - Mu) / Sigma
        elseif Loi(i) == 1 then
            U(:,i) = (X(:,i) - Moy(i)) / Stdev(i)
        // Loi Lognormale U = (ln(X) - MuLog) / sigmaLog [avec MuLog et SigmaLog paramètres internes calculés à partir de Mu et Sigma]
        elseif Loi(i) == 2 then
            mulog = log((Moy(i) ** 2)/(sqrt(Moy(i) ** 2 + Stdev(i) ** 2)))
            siglog = sqrt(log(1. + Stdev(i) ** 2 / Moy(i) ** 2))
            U(:,i) = (log(X(:,i)) - mulog) / siglog
        // Loi Uniforme U = Phi^{-1}((X-a)/(b-a)) [avec a et b bornes inf et sup calculées à partir de Mu et Sigma]
        elseif Loi(i) == 3 then
            a = Moy(i) - sqrt(3) * Stdev(i)
            b = Moy(i) + sqrt(3) * Stdev(i)
            U(:, i) = cdfnor("X", zeros(X(:,i)), ones(X(:,i)), (X(:,i) - a) / (b - a), 1. - (X(:,i) - a) / (b - a))
        // Loi de Weibull U = f(X, lambdaW, kW) [avec lambdaW et kW paramètres internes calculés à partir de Mu et Sigma]
        elseif Loi(i) == 4 then
            ratio = 1.0 + (Stdev(i) / Moy(i)) ** 2.
            betaMin = 1.
            betaMax = 1.
            step = 0.5
            epsilon = 1.e-12

            if ratio > 2 then
                betaMin = betaMin - step
                step = step * 0.5
                t = exp(log(gamma(1.0 + 2.0 / betaMin)) - 2.0 * log(gamma(1.0 + 1.0 / betaMin)))
                while t < ratio
                    betaMin = betaMin - step
                    step = step * 0.5
                    t = exp(log(gamma(1.0 + 2.0 / betaMin)) - 2.0 * log(gamma(1.0 + 1.0 / betaMin)))
                end
                betaMax = betaMin + 2.0 * step

            else
                betaMax = betaMax + step
                step = step * 2.0
                t = exp(log(gamma(1.0 + 2.0 / betaMax)) - 2.0 * log(gamma(1.0 + 1.0 / betaMax)))
                while t >= ratio
                    betaMax = betaMax + step
                    step = step * 2.0
                    t = exp(log(gamma(1.0 + 2.0 / betaMax)) - 2.0 * log(gamma(1.0 + 1.0 / betaMax)))
                end
                betaMin = betaMax - 0.5 * step
            end

            while 1
                kW = 0.5 * (betaMin + betaMax)
                if (betaMax - betaMin <= epsilon * (1.0 + abs(betaMax + betaMin))) then

                    lambdaW = Moy(i) / gamma(1.0 + 1.0 / kW)
                    break
                end
                t = exp(log(gamma(1.0 + 2.0 / kW)) - 2.0 * log(gamma(1.0 + 1.0 / kW)))
                if (t < ratio) then
                    betaMax = kW
                else
                    betaMin = kW
                end
            end

            U(:,i) = cdfnor("X", zeros(X(:,i)), ones(X(:,i)), 1. - exp(- (X(:,i) / lambdaW) ** kW), exp(- (X(:,i) / lambdaW) ** kW))
        // Loi de Gumbel
        elseif Loi(i) == 5 then
            cstgamma = 0.5772;
            Beta = sqrt(6)*Stdev/%pi;
            Mu = Moy - Beta*cstgamma;
            
            U(:,i) = cdfnpr("X", zeros(X(:,i)), ones(X(:,i)), exp(-exp((Mu-X(:,i))/Beta)), 1-exp(-exp((Mu-X(:,i))/Beta)));
            
        end
    end
endfunction

// Fonction de performance
function P=Perf(U,Seuil)
    X=tisoinv(U,Loi,Moy,Stdev);
    fonctionG = G(X);
//    non vectorisé
//    sX = size(X)
//    for i=1:sX(1)
//        fonctionG(i,1)=G(X(i,:));
//    end
    P=fonctionG-Seuil;
endfunction


// Construction du modèle de krigeage avec verification de la bonne interpolation
// sur les points du modèle
function Hmodel = BuildKrigingModel(S, H_S, theta0, lower_b, upper_b)
    predictok = %f;
    flag_sortie = %f;
    relative_error = 1e-2
    while ~predictok
        [Hmodel,Hperf] = dacefit(S, H_S, regpoly0, corrgauss, theta0,lower_b,upper_b);
        [mu_S, var_S] = predictor(S, Hmodel);
        if or(abs((mu_S - H_S)./H_S) > relative_error) then // si la prédiction est mauvaise sur les points du DOE
            theta0 = theta0 / 10  // changement du point de départ  initial
            if or(theta0 <= lower_b) then // si theta0 a atteint la borne inférieure
                theta0 = upper_b /10
                if flag_sortie then // Sécurité si theta0 a déjà fait un tour
                    disp('Attention, mauvaise construction du metamodèle')
//                    set(boutonstart,'fontSize',16);
//                    set(boutonstart,'FontWeight','bold');
//                    set(boutonarret,'fontSize',10);
//                    set(boutonarret,'FontWeight','light');
//                    xinfo('Prêt');
                    abort;
                end
                flag_sortie = %t;
            end
            disp('Searching for theta0...');
        else
            predictok = %t;
        end
    end
endfunction


///////////////////////////////////////////////////////////////////////////////////////
// Fonctions FORM
//////////////////////////////////////////////////////////////////////////////////////

function [G0,grad]=gradiantGFORM(U)
  G0=Perf(U,Seuil)
  h=0.01;
  for i=1:nva
    UplusH=U;
    UplusH(:,i)=UplusH(:,i)+h;
    GplusH=Perf(UplusH,Seuil);
    grad(:,i)=(GplusH-G0)/h;
  end
endfunction

// Algorithme de RF 
function [BetaFORM,Ustar,iter]=RAF(Udepart)
  [G0,grad0]=gradiantGFORM(Udepart);
  beta0=0;
  iter=1;
  bet=[];
  U=[];
  while 1
     norme_gradiant0=sqrt(sum(grad0.^2,2))
     alpha0=grad0.*((norme_gradiant0*ones(1,nva)).^(-1));
     beta1=-sum(Udepart.*alpha0)+G0.*(norme_gradiant0.^(-1));
     u1=-beta1*ones(1,nva).*alpha0
     if abs(beta1-beta0)<0.01 | iter=100 then
       BetaFORM=beta1;
       alpha=alpha0;
       Ustar=u1;
       if iter==100 then
           disp('ATTENTION : nombre limites d''itération atteint, problème de convergence')
       end
       break;
     else
       bet=[bet beta1];
       beta0=beta1;
       iter=iter+1;
       Udepart=u1;
       [G0,grad0]=gradiantGFORM(Udepart);
     end
  end
  
endfunction



%%% Critère de Routh %%%

function[stable, gauche, droite] = routh(coefficients)
    nbr_coeff = length(coefficients);                       % Longueur des coéfficients
    lignes = nbr_coeff;
    colonnes = round(nbr_coeff/2)+1;
    if mod(nbr_coeff, 2) == 1
        coefficients = cat(2, coefficients, [0 0]);
    elseif mod(nbr_coeff, 2) == 0
        coefficients = cat(2, coefficients, 0);
    end 

    % construction des vecteurs de routh initiaux
    routh_array = zeros(lignes, colonnes);
    routh_array(1,1:colonnes-1) = coefficients(1:2:nbr_coeff);
    routh_array(2,1:colonnes-1) = coefficients(2:2:nbr_coeff);

    fprintf("Le matrice de routh",routh_array);
    

    % construction des autres lignes
    %reste_lignes = 3:lignes;
    iterCol = 1;
    iterLig = 3;
    suivant = true;
    while iterLig < lignes+1 && suivant
        
        for i = iterCol:colonnes-1
            ligneSuivante = -det([routh_array(iterLig-2,1) routh_array(iterLig-1, 1) ;
            routh_array(iterLig-2,i + 1) routh_array(iterLig-1, i+1)])/(routh_array(iterLig-1,1));
            routh_array(iterLig, i) = ligneSuivante;
            
        end
        verification = routh_array(iterLig) == zeros(1, colonnes)
        if sum(verification) == colonnes
            ligneSuivante2  = ligneNulle(routh_array(iterLig-1, 1:colonnes), iterLig, lignes);
            routh_array(iterLig, 1:colonnes-1) = ligneSuivante2;
        end
        if (routh_array(iterLig,1) == 0) && sum(verification)  ~= colonnes
            [dpg, dpd] = routhNulle(routh_array, iterLig, iterCol, lignes, colonnes);
            suivant = false;
            break
        end
        iterLig = iterLig+1;
    end
    if suivant == true
        dpg = routh_array(1:lignes, 1:colonnes-1);
        dpd = [];
    end

    signChdpg = 0;
    signdpg = 1;
    signChdpd = 0;
    signdpd = 0;
    stable = true;
    droite = 0;
    gauche = lignes - 1;
    if isempty(dpd)
        col_dpg = dpg(1:lignes,1)';
        for i = 1:lignes
            if col_dpg(i) < 0 && signdpg == 1
                signdpg = 0;
                signChdpg = signChdpg + 1;
            end
            if col_dpg(i) > 0 && signdpg == 0
                signdpg = 1;
                signChdpg = signChdpg + 1;
            end
        end
        if signChdpg == 0
            stable = true;
            droite = 0;
        end
        if signChdpg > 0
            stable = false;
            droite = signChdpg;
            gauche = gauche - droite;
        end
    end
    if  ~isempty(dpd)
        col_dpg = dpg(1:lignes,1)';
        col_dpd = dpd(1:lignes,1)';
        for i = 1:lignes
            if col_dpg(i) < 0 && signdpg == 1
                signdpg = 0;
                signChdpg = signChdpg + 1;
            end
            if col_dpg(i) > 0 && signdpg == 0
                signdpg = 1;
                signChdpg = signChdpg + 1;
            end

            if col_dpd(i) < 0 && signdpd == 1
                signdpd = 0;
                signChdpd = signChdpd + 1;
            end
            if col_dpd(i) > 0 && signdpd == 0
                signdpd = 1;
                signChdpd = signChdpd + 1;
            end
        end  
    end    
end

% fonction de gérer le cas où un zéro dans la serie de routh
function [routh_array1, routh_array2] = routhNulle(routh, iterLigf, iterColf, lignesf, colonnesf)
        r1 = routh;
        r2 = routh;
        r1(iterLigf, 1) = 0.001;
        r2(iterLigf, 1) = -0.001;
        iterLig1 = iterLigf+1;
        while iterLig1 < lignesf+1 
            for i = iterColf:colonnesf-1
                ligneSuivante1 = -det([r1(iterLig1-2,1) r1(iterLig1-1, 1) ;
                r1(iterLig1-2,i + 1) r1(iterLig1-1, i+1)])/(r1(iterLig1-1,1));
                r1(iterLig1, i) = ligneSuivante1;

                ligneSuivante2 = -det([r2(iterLig1-2,1) r2(iterLig1-1, 1) ;
                r2(iterLig1-2,i + 1) r2(iterLig1-1, i+1)])/(r2(iterLig1-1,1));
                r2(iterLig1, i) = ligneSuivante2;
            end
            verification = sum(routh(iterLigf, 1:colonnesf) == zeros(colonnesf));
        
        if (sum(verification) == 4)
            break
        end
            iterLig1 = iterLig1+1;
        end
        routh_array1 = r1;
        routh_array2 = r2;
end

% fonction de gérer le cas où une ligne est nulle
function [ligneSuivante] = ligneNulle(routh, iterLigf, lignesf)
    exposant = lignesf-1:-1:0 ;
    if mod(iterLigf, 2 ) == 1
        exposant = cat(2, exposant, [0 0]);
    elseif mod(iterLigf, 2 ) == 0
        exposant = cat(2, exposant, 0);
    end
    exposant = exposant(iterLigf-1:2:length(exposant));
    routh = routh(1:(length(routh)-1));
    ligneSuivante = exposant .* routh;
end
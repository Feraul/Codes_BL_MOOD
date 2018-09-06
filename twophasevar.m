%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 03/05/2012
%Modify data:  / /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION calculate the 

%--------------------------------------------------------------------------

function [fw,fo,gama,krw,kro,nw,no,kwmax,komax] = twophasevar(Sw,...
    numcase)
%Define global parameters:
global satlimit visc;

%Initialize the parameters:
krw = zeros(length(Sw),1);

%"gama" is equal to "fw*kro/visc(2)"
gama = krw;

%Calculate the parameters to each element
    %Special case (Hurtado et al., 2007)
    %Cases 34.3 and 34.4 ==> 1 CV sourse; 46.3 ==> 4 CV source 
    %(boundary free)
    if numcase == 34.3 || numcase == 34.4 || numcase == 46.3
        %Define Mobility ratio:
        M = visc(2)/visc(1);
        %Definition of fractional flow (WATER)
        fw = Sw.^2;
        %Definition of fractional flow (OIL)
        fo = 1 - fw;
        %Definition of relative permeability (WATER)
        krw = fw./(M.*(1 - fw) + fw); 
        %Definition of relative permeability (OIL)
        kro = 1 - krw;
    
    %Adapted from Bastian (2002) for the benchmark 31.1, with lambda = 2; 
    elseif numcase == 31.1 || numcase == 31.6 
        %Normalizes the saturation:
        Swn = ((Sw - satlimit(1))./(1 - satlimit(1) - satlimit(2))); 
        
        %Definition of relative permeability (WATER)
        krw = Swn.^4; 
        %Definition of relative permeability (OIL)
        kro = ((1 - Swn).^2).*(1 - (Swn.^2));

        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
        fw = (krw./visc(1))./((krw./visc(1)) + (kro./visc(2)));
        %Definition of fractional flow (OIL)
        fo = (kro./visc(2))/((krw./visc(1)) + (kro./visc(2)));
    
        %------------------------------------------------------------------
        %Define "gama". It is used when gravity effects are account
    
        gama = fw.*kro./visc(2);
        
    %Another examples
    else
        %It Chose the properties according to "numcase"
        %Example 36: Two-Phase Flow case. Lamine and Edwards, 2010.
        booleanlin = (numcase == 36 || numcase == 36.1);
        %Example 45: Two-Phase Flow case. Kozdon et al., 2011. 
            %"Multidimensional upstream weighting for multiphase transport 
            %in porous media". Section 5.2 the benchmark 45 to 46.
            boolean4th = ((numcase >= 45 && numcase < 46) || ...
                numcase == 43.1);
            %Example 34.8: Two-Phase Flow case. Adapted from Nikitin and 
            %Vassilevisky, 2012. 
            %It evaluates three different quadrangular meshes.
            boolean5th = (numcase == 34.8);
            %Any another case
            boolean2nd = (booleanlin + boolean4th + boolean5th) == 0; 
            %Get the exponent:
            nw = booleanlin + 4*boolean4th + 5*boolean5th + 2*boolean2nd; 
            no = booleanlin + 2*boolean4th + boolean5th + 2*boolean2nd;
            %Fit parameter (water and oil)
            kwmax = 1;
            komax = 1;
            
        %------------------------------------------------------------------
        %Normalized Saturation
    
        Swn = ((Sw - satlimit(1))./(1 - satlimit(1) - satlimit(2))); 
    
        %------------------------------------------------------------------
        %Relative Permeability:
    
        %Definition of relative permeability (WATER)
        krw = kwmax*(Swn).^nw; 
        %Definition of relative permeability (OIL)
        kro = komax*(1 - Swn).^no; 
    
        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
        fw = (krw./visc(1))./((krw./visc(1)) + (kro./visc(2)));
        %Definition of fractional flow (OIL)
        fo = (kro./visc(2))./((krw./visc(1)) + (kro./visc(2)));
    
        %------------------------------------------------------------------
        %Define "gama". It is used when gravity effects are account
    
        gama = fw.*kro./visc(2);
    end  %End of IF (special case)

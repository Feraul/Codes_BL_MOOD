%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: FUNCTION
%Criate date: 13/03/2013
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Calculate the water saturation field using an explicit formulation. 

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

function [d2fdS,dgamadS] = calcder2dS(fw,gama,Sw,dertype)
%Define global parameters
global satlimit visc numcase keygravity;

%Define "d2fdS" and "dgamadS" according to "dertype"
switch dertype
    %Discrete derivative. It needs of two saturation value 
    %(used in timestep definition, for example)
    case 0
        %There saturation difference bigger than zero
        if (Sw(1) - Sw(2)) ~= 0
            d2fdS = (fw(1) - fw(2))/(Sw(1) - Sw(2));
            dgamadS = (gama(1) - gama(2))/(Sw(1) - Sw(2));
        %The saturation difference is zeros
        else
            d2fdS = 0;
            dgamadS = 0;
        end  %End of IF
    %Analitical derivative. It uses only one value into argument.
    case 1
        %Initialize some propoerties (two-phase flow)
        Swi = satlimit(1);
        Sor = satlimit(2);
        miw = visc(1);
        mio = visc(2);
        %Initialize "d2fdS" and "dgamadS"
        d2fdS = zeros(length(Sw),1);
        dgamadS = d2fdS; 
        
        %------------------------------------------------------------------
        %Chose analitical derivative according to "numcase"
        
        %Adapted from Bastian (2002) for the benchmark 31.1 
        %(van Genushten model); 
        if numcase == 31.1 || numcase == 31.6 
            %for i = 1:length(Sw)
                %Calculate "d2fdS"
                d2fdS = ...
               -(4*(Sw.^2).*miw.*mio.*((Sw.^4).*(5 - 6.*Sw + (Sw.^3))*mio ...
               - ((Sw - 1).^4).*(3 + 4*Sw + 4.*(Sw.^2) + (Sw.^3))*miw))./...
               (((Sw.^4)*mio - ((Sw - 1).^3).*(Sw + 1)*miw).^3) ;
            %end  %End of FOR
        %Hurtado et al. (2007)
        elseif numcase == 34.3
            d2fdS = 2;
            dgamadS = 0;
        %Another examples. (Brooks-Corey model). It can vary the 
        %coefficients "nw", "no", ...
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

            %Calculate the derivate
            for i = 1:length(Sw)
                %Define some terms:
                term1 = (Sw(i) + Sor - 1);
                term2 = (Swi + Sor - 1);
                term3 = (Swi - Sw(i));
                term4 = (term1/term2)^no;
                term5 = (term3/term2)^nw;
                term6 = (komax*miw*term4 + kwmax*mio*term5);
                denom = ((term1^2)*((-term3)^2)*(term6^3)) + 1e-16; 
                %Calculate "dfwdSw"
                d2fdS(i) = (komax*kwmax*miw*mio*term4*...
                    term5*(-(nw - no)*term1*term3*term6 - ...
                    term1*(nw*term1 + no*term3)*term6 + nw*term1*...
                    (nw*term1 + no*term3)*term6 + term3*(nw*term1 + ...
                    no*term3)*term6 - no*term3*(nw*term1 + no*term3)*...
                    term6 + 2*term1*term3*(nw*term1 + ...
                    no*term3)*(((komax*miw*no*term4)/term1) - ...
                    ((kwmax*mio*nw*term5)/term3))))/denom;

                %Calculate the gravity contribution (when necessary)
                if strcmp(keygravity,'y')
                    %Define another terms:
                    %Calculate "dgamadSw"
                    dgamadS(i) = (komax*kwmax*term4*...
                        term5*(komax*miw*nw*term1*...
                        term4 - kwmax*mio*no*term3*...
                        term5))/denom;
                end  %End of IF
            end  %End of FOR
        end  %End of IF (is the Bastian example?)
end  %End of SWITCH


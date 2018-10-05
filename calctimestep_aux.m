%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 08/05/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Determinate the saturation and presure fields (2D) in a eithe homogen or 
%heterogen domain such as isotropic and anisotropic media for each time 
%step or in the steady state when will be important.  

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

function [dt] = calctimestep_aux(flowrate,Fg,Sw,satinbound,injecelem,klb)
%Define global parameters:
global pormap elemarea courant order inedge bedge normals numcase ...
    smethod;

%Define tolrance
tol = 1e-12;
%Define the degree of the reconstruction polynomium "n"
n = order - 1;

%Get the flow rate in whole edge when MULTIDIMENSIONAL Schemes are used
%Multidimens. Applic. (Lamine and Edwards, 2010; Kozdon et al.,2011)
if strcmp(smethod,'mwic') || strcmp(smethod,'mwec') || ...
        strcmp(smethod,'rtmd')
    %Join the flowrate calculated for each half-edge in a unic flowrate
    %over the whole edge.
    [flowrate] = joinflowrate(flowrate);
end  %End of IF

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Swept all internal edges
i = 1:inedgesize;
%Get the fractional flux on the left:
[fw_left,~,gama_left,] = twophasevar(Sw(inedge(:,3)),numcase);
%Get the fractional flux on the right:
[fw_right,~,gama_right,] = twophasevar(Sw(inedge(:,4)),numcase);
    
%Obtain the apropriated deltax:
%Calculate middle volume: the mean between volume shared by edge 
vol = (elemarea(inedge(:,3)) + elemarea(inedge(:,4)))/2;
    
%Define delta t:
%Chose according physical effects (gravity existence etc)

%Calculate the derivative of functions "fw" and "gama"
%Initialize "dfwdS" and "dgamadS" (if necessary)
dfwdS = zeros(inedgesize,1);
%Calculate "dfwdS"
diff_fw = fw_left - fw_right;
diffsat = Sw(inedge(:,3)) - Sw(inedge(:,4));
%Get the accuracy
diffsat = diffsat.*(abs(diffsat) > tol);
pointnozero = diffsat ~= 0;
%Fill "dfwdS"
dfwdS(pointnozero) = diff_fw(pointnozero)./diffsat(pointnozero); 

%There is gravity effects
if size(Fg,2) > 1
    %Initialize "dgamadS"
    dgamadS = zeros(inedgesize,1); 
    %Calculate "dgamadS"
    diff_gama = gama_left - gama_right;
    %Fill "dgamadS"
    dgamadS(pointnozero) = diff_gama(pointnozero)./diffsat(pointnozero); 
    
    %Calculate "dt" by edge (inedge)
    dtbyedge = abs((courant/((2*n) + 1))*pormap.*vol./...
        (dfwdS.*flowrate(bedgesize + i) + ...
        dgamadS.*dot(Fg(inedge(:,3),:)',normals(bedgesize + i,1:2)')' + ...
        1e-16));
%There is no gravity effects
else
    %Calculate "dt" by edge (inedge)
    dtbyedge = abs((courant/((2*n) + 1))*pormap.*vol./...
        (dfwdS.*flowrate(bedgesize + i) + 1e-16));
end  %End of IF

%--------------------------------------------------------------------------
%Boundary Tratment (Buckley-Leverett Applications)

if any(klb)
    %Swept edges in "bedge" associated with boundary (injection)
    i = 1:length(klb);
    %Get the fractional flux on the left:
    [fw_bound,~,gama_bound,] = twophasevar(satinbound(i),numcase);
    %Get the fractional flux on the right:
    [fw_left,~,gama_left,] = twophasevar(Sw(injecelem),numcase);

    %Calculate middle volume: the mean between volume shared by edge 
    vol = elemarea(injecelem);

    %Initialize "dfwdS" and "dgamadS" (if necessary)
    dfwdS = zeros(length(injecelem),1);
    %Calculate "dfwdS"
    diff_fw = fw_bound - fw_left;
    diffsat = satinbound(i) - Sw(injecelem);
    diffsat = diffsat.*(abs(diffsat) > tol);
    pointnozero = diffsat ~= 0;
    %Fill "dfwdS"
    dfwdS(pointnozero) = diff_fw(pointnozero)./diffsat(pointnozero);
    
    %Define delta t:
    %Chose according physical effects (gravity existence etc)
    %There is gravity effects
    if size(Fg,2) > 1
        %Initialize "dgamadS"
        dgamadS = zeros(length(injecelem),1); 
    
        %Calculate "dgamadS"
        diff_gama = gama_bound - gama_left;
        %Fill "dgamadS"
        dgamadS(pointnozero) = ...
            diff_gama(pointnozero)./diffsat(pointnozero); 
        
        %Calculate "dt" by edge (inedge)
        dtbyboundedge = abs((courant/((2*n) + 1))*pormap.*vol./...
            (dfwdS.*flowrate(klb(i)) + dgamadS.*dot(Fg(injecelem(i),:)',...
            normals(klb(i),1:2)')' + 1e-16));
    %There is no gravity effects
    else
        %Calculate "dt" by edge (inedge)
        dtbyboundedge = abs((courant/((2*n) + 1))*pormap.*vol./...
            (dfwdS.*flowrate(klb(i)) + 1e-16));
    end  %End of IF

    %Do the union between "dtbyedge" and "dtbyboundedge".
    dtbyedge = [dtbyedge; dtbyboundedge];
end  %End of IF (boundary contribution)

%Finally, define the minor "dt"
dt = min(dtbyedge(dtbyedge ~= 0));

        
%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with
%high-order resolution
%Type of file: FUNCTION
%Criate date: 13/01/2013
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%.

%--------------------------------------------------------------------------
%Additional comments:
%

%--------------------------------------------------------------------------

function [advecterm,entrineqterm,earlysw] = calcnumflux(Sw,Fg,flowrate,...
    taylorterms,limiterflag,flagknownvert,satonvertices,flagknownedge,...
    satonboundedges,pointbndedg,pointinedg,orderbedgdist,orderinedgdist,...
    constraint,mlplimiter,earlysw)
%Define global parameters:
global coord elem bedge inedge normals dens numcase centelem order courant;

%Initialize a tolerance. It is a computational zero
tol = 1e-12;
%Initialize "advecterm" and "entrineqterm"
advecterm = zeros(size(elem,1),1);
entrineqterm = advecterm;
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);

paramk = 1;
lenginjface = 1.166822;

%--------------------------------------------------------------------------
%Boundary edges (when it exists):

%In cases where the producer well edges are evaluated, it is necessary
%verify if the producer well shares some boundary edge.
if any(pointbndedg)
    %Swept "bedge"
    for i = 1:length(pointbndedg)
        %Initialize some parameters:
        ibedg = pointbndedg(i);
        %Define the order for this edge.
        faceorder = orderbedgdist(i);
        
        %Define the "vertices"
        vertices = bedge(ibedg,1:2);
        %Get the coordinate for the vertices:
        verticescoord = coord(vertices,:);
        %Define left elements
        leftelem = bedge(ibedg,3);
        
        %------------------------------------------------------------------
        %Define velocity due gravity
        
        %There is gravity
        if size(Fg,2) > 1
            dotvg = dot(Fg(leftelem,:),normals(ibedg,1:2))*(dens(1) - ...
                dens(2))/lenginjface;
            %There is NO gravity
        else
            dotvg = 0;
        end  %End of IF
        
        %------------------------------------------------------------------
        
        %Define the elements that share the edge evaluated
        elemeval = leftelem;
        
        %Verify if there is saturation prescribed on boundary:
        %There is a prescribed saturation
        if flagknownedge(ibedg) == 1
            %Attribute the saturation on boundary
            Sleft = satonboundedges(ibedg);
            if Sleft>1
                Sleft=1;
                
            end
            Sleft = Sleft*(Sleft >= 0);
            %There is no prescribed saturation. It is necessary calculate.
        else
            %Get the statment regarding to mlp limiter:
            boolmlp = (length(mlplimiter) == 1);
            mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
            
            %Get the saturation value recovered
            Sleft = getsatonedge(elemeval,vertices,verticescoord,...
                taylorterms,Sw,limiterflag,faceorder,constraint,flagknownvert,...
                satonvertices,mlpbyelem,centelem(leftelem,1:2));
            
            %======================================================================
            limiterflagaux=limiterflag;
            limiterflagaux{1}='on';
            limiterflagaux{4}='wf';
            limiterflagaux{7}=0.33333;
            %======================================================================
            Sleftaux = getsatonedge(elemeval,vertices,verticescoord,...
                taylorterms,Sw,limiterflagaux,faceorder,constraint,flagknownvert,...
                satonvertices,mlpbyelem,centelem(leftelem,1:2));
            if Sleft>1 || Sleft<0
                Sleft=Sleftaux;
                
            end
            Sleft = Sleft*(Sleft >= 0);
            
        end  %End of IF
        
        Sleft = Sleft*(Sleft >= 0);
        %Fill "earlysw"
        earlysw(ibedg) = Sleft;
        
        %Calculate the fractional flow in boundary ("fwbound")
        [fw,~,gama,] = twophasevar(Sleft,numcase);
        
        %Define the normal velocity into face
        dotvn = flowrate(ibedg);
        %Get accuracy for "dotvn"
        dotvn = dotvn*(abs(dotvn) > tol);
        
        %     [dfwdSleft,] = calcdfunctiondS(0,0,Sleft,1);
        %     wcmax = dfwdSleft*dotvn;
        %     difusion = wcmax*sign(dotvn)*abs(Sleft - Sw(leftelem))/2.4;
        
        %Calculate the numerical flux through interface.
        numflux = dotvn*fw + dotvg*gama;
        %Obtain the contribution of interface over element to LEFT
        advecterm(leftelem) = advecterm(leftelem) + numflux;
        
        %Evaluate the ENTROPY condition (flux)
        %Get the entropy flux (left value)
        % entrvar = getineqentropy(Sleft,5,2,1);
        %Calculate the entropy flux (right value)
        % entrflux = entrvar*dotvn;
        %Obtain the contribution of interface over element to LEFT
        %  entrineqterm(leftelem) = entrineqterm(leftelem) + entrflux;
    end  %End of FOR (Swept "bedge")
end  %End of IF (Does evaluate the boundary edges?)

%--------------------------------------------------------------------------
%Internal edges:

%Swept "inedge" evaluating left and right elements by edge. Apply
%approximated Riemann Solver through edge.

for i = 1:length(pointinedg)
    %Initialize some parameters:
    inedg = pointinedg(i);
    
    %Define "vertices"
    vertices = inedge(inedg,1:2);
    %Get the coordinate for the vertices:
    verticescoord = coord(vertices,:);
    %Define left and right elements
    leftelem = inedge(inedg,3);
    rightelem = inedge(inedg,4);
    
    %Left Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,1);
    %Define the elements that share the edge evaluated
    elemeval = [leftelem rightelem];
    %----------------------------------------------------------------------
    %Calculate the velocity due to GRAVITY effect
    
    %There is gravity
    if size(Fg,2) > 1
        dotvg = dot(Fg(leftelem,:),normals(bedgesize + inedg,1:2))*...
            (dens(1) - dens(2))/lenginjface;
        %There is NO gravity
    else
        dotvg = 0;
    end  %End of IF
    
    %----------------------------------------------------------------------
    
    %Get the saturation value recovered on each quadrature point ("on_q")
    %Get the statment regarding to mlp limiter:
    boolmlp = (length(mlplimiter) == 1);
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    %Left Contribution:
    Sleft = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
        limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
    %======================================================================
    limiterflagaux=limiterflag;
    limiterflagaux{1}='on';
    limiterflagaux{4}='wf';
    limiterflagaux{7}=0.33333;
    %======================================================================
    Sleftaux = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
        limiterflagaux,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
    %Right Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,2);
    %Define the elements that share the edge evaluated
    elemeval = [rightelem leftelem];
    %Get the statment regarding to mlp limiter:
    mlpbyelem = mlplimiter(boolmlp + (1 - boolmlp)*elemeval(1),:);
    %Get the saturation value recovered on each quadrature point ("on_q")
    Sright = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
        limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
    Srightaux = getsatonedge(elemeval,vertices,verticescoord,taylorterms,Sw,...
        limiterflagaux,faceorder,constraint,flagknownvert,satonvertices,...
        mlpbyelem,centelem(elemeval,1:2));
    %Get an accurated value:
    if Sleft>1 || Sleft<0
        %Sleft=Sw(leftelem);
        Sleft = Sleftaux;
    end
    
    if Sright>1 || Sright<0
        %Sright=Sw(rightelem);
        Sright = Srightaux;
    end
    
    Sleft = Sleft*(abs(Sleft) > tol);
    Sright = Sright*(abs(Sright) > tol);
    
    %PAD
    
    Sleft = Sleft*(Sleft >= 0);
    Sright = Sright*(Sright >= 0);
    
    %     if i==1
    %         if Sright<Sleft
    %
    %             Sleft=Sw(leftelem);
    %             %Sleft = Sleftaux;
    %             if Sright<Sleft
    %                 Sright=Sw(rightelem);
    %                 %Sright = Srightaux;
    %             end
    %
    %         else
    %             Sleft=Sleft;
    %         end
    %     else
    %         if Sleft<Sright
    %             Sright=Sw(rightelem);
    %             %Sright = Srightaux;
    %
    %             if Sleft<Sright
    %
    %                 Sleft=Sw(leftelem);
    %                 %Sleft = Sleftaux;
    %             end
    %         else
    %             Sright=Sright;
    %         end
    %     end
    
    
    %Calculate the derivative dfw/dSw:
    %Analytical:
    [dfwdS_analeft,dgamadS_analeft] = calcdfunctiondS(0,0,Sleft,1);
    [dfwdS_analright,dgamadS_analright] = calcdfunctiondS(0,0,Sright,1);
    
    %Get accuracy (second derivative):
    dfwdS_analeft = dfwdS_analeft*(abs(dfwdS_analeft) > tol);
    dfwdS_analright = dfwdS_analright*(abs(dfwdS_analright) > tol);
    
    %Get the analytical for a point between "Sleft" and "Sright"
    Smid = (Sleft + Sright)/2;
    %Smid = ((1-courant)*Sleft + courant*Sright);
    
    %Discrete:
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    [fw,~,gama,] = twophasevar([Sleft Smid Sright],numcase);
    % fw(1) ---> fluxo fracional no elemento fw(Sleft)
    % fw(2) ---> fluxo fracional no elemento fw(Smid)
    % fw(3) ---> fluxo fracional no elemento fw(Sright)
    
    %Calculate the Rankine-Hugoniot ratio:
    % veja o tese Darlan pag. 165
    [dfwdS_rh,dgamadS_rh] = ...
        calcdfunctiondS([fw(1) fw(3)],[gama(1) gama(3)],[Sleft Sright],0);
    %%
    %There saturation difference bigger than zero
    %     if (Sleft - Smid) ~= 0
    %         dfdSlef = (fw(2) - fw(3))/( Smid-Sleft);
    %
    %         %The saturation difference is zeros
    %     else
    %         dfdSlef = 0;
    %
    %     end  %End of IF
    %     if (Smid-Sright) ~= 0
    %         dfdSrel= (fw(2) - fw(1))/(Smid-Sright);
    %
    %         %The saturation difference is zeros
    %     else
    %         dfdSrel = 0;
    %
    %     end  %End of IF
    %               if ( dfdSrel< dfwdS_rh || dfdSrel== dfwdS_rh)  && (dfwdS_rh< dfdSlef|| dfwdS_rh==dfdSlef )
    %              else
    %                  i
    %      %            disp('Não satisfaz condição de Rankine Hugoniot')
    %                  pause
    %              end
    
    %%
    %Get accuracy (RH ratio):
    dfwdS_rh = dfwdS_rh*(abs(dfwdS_rh) > tol);
    
    %Calculate the second derivative d2fw/dSw2:
    [d2fwdS2_discleft,d2gamadS2_discleft] = calcder2dS(0,0,Sleft,1);
    [d2fwdS2_discright,d2gamadS2_discright] = calcder2dS(0,0,Sright,1);
    %Get accuracy (second derivative):
    d2fwdS2_discleft = d2fwdS2_discleft*(abs(d2fwdS2_discleft) > tol);
    d2fwdS2_discright = d2fwdS2_discright*(abs(d2fwdS2_discright) > tol);
    d2gamadS2_discleft = d2gamadS2_discleft*(abs(d2gamadS2_discleft) > tol);
    d2gamadS2_discright = d2gamadS2_discright*(abs(d2gamadS2_discright) > tol);
    
    %---------------------------------------
    %Define the normal velocity in each face
    dotvn = flowrate(bedgesize + inedg);
    
    %Get accuracy:
    dotvn = dotvn*(abs(dotvn) > tol);
    dotvg = dotvg*(abs(dotvg) > tol);
    %---------------------------------------
    
    %Get the sign of the first derivative:
    signder_left = sign(dfwdS_analeft*dotvn + dgamadS_analeft*dotvg);
    signder_right = sign(dfwdS_analright*dotvn + dgamadS_analright*dotvg);
    %Get the sign of second derivative:
    sign2der_left = sign(d2fwdS2_discleft*dotvn + d2gamadS2_discleft*dotvg);
    sign2der_right = sign(d2fwdS2_discright*dotvn + d2gamadS2_discright*dotvg);
    
    %Define "charvel" (characteristic velocity):
    %Define "charvel" on the left
    %   charvel_left = dotvn*dfwdS_analeft;
    %Define "charvel" on the right
    %   charvel_right = dotvn*dfwdS_analright;
    
    %Discrete:
    %Discrete:
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    %    [fwaux,fo,gama,] = twophasevar([Sleft Smid Sright],numcase);
    %"fw" has three values: fw(Sleft) is fw(1), fw(Sright) is fw(3)
    %    [dfwdS_discleft,~] = calcdfunctiondS([fwaux(1) fwaux(2)],[gama(1) gama(2)],[Sleft Smid],0);
    %    [dfwdS_discright,~] = calcdfunctiondS([fwaux(2) fwaux(3)],[gama(2) gama(3)],[Smid Sright],0);
    %     %Define "numcharvel" (numerical characteristic velocity):
    %     %Define "charvel" on the left
    %    numcharvel_left = dotvn*dfwdS_discleft;
    %     %Define "charvel" on the right
    %     numcharvel_right = dotvn*dfwdS_discright;
    %Define the Rankine-Hugoniout velocity
    charvel_rh = dotvn*dfwdS_rh + dotvg*dgamadS_rh;
    
    %MOOD is turned ON
    %     if strcmp(limiterflag{8},'on')
    %          %Verify the entropy condition for the left and right states:
    %          %"entropycond" == "0", MOOD acts (get the order down);
    %          %"entropycond" == "1", MOOD remains with higher order;
    %          entropycond = (numcharvel_left*(sign2der_left) >= charvel_rh && ...
    %              charvel_rh >= numcharvel_right*(sign2der_right));
    %          %Verify NON physical condition:
    %          booleanpad = (Sleft < 0 || Sright < 0 || Sleft > 1 || Sright > 1);
    %          entropycond = (1 - booleanpad)*entropycond*(entropycond == 1);
    %      end  %End of IF (is the MOOD turned ON?)
    
    %----------------------------------------------------------------------
    %Choise according second derivative sign (see Serma, 2009)
    
    %It use Upwind (Roe)
    if signder_left*signder_right >= 0 && sign2der_left*sign2der_right >= 0
        %Verify the sign of the characteristic velocity:
        %It uses the saturation on the left
        if charvel_rh > 0 || charvel_rh==0
            %Calculate the numerical flux through interface
            numflux = fw(1)*dotvn + gama(1)*dotvg;
            %Fill "earlysw"
            earlysw(bedgesize + inedg) = Sleft;
            
            %Entropy:
            %Get "entrvar":
            %   entrvar = getineqentropy(Sleft,5,2,(1 - entropycond)*paramk);
            %Calculate the entropy flux (left value)
            %   entrflux = entrvar*dotvn;
            %It uses the saturation on the right
        else
            
            %Calculate the numerical flux through interface
            numflux = fw(3)*dotvn + gama(3)*dotvg;
            %Fill "earlysw"
            earlysw(bedgesize + inedg) = Sright;
            
            %Entropy:
            %Get "entrvar":
            %      entrvar = getineqentropy(Sright,5,2,(1 - entropycond)*paramk);
            %Calculate the entropy flux (right value)
            %      entrflux = entrvar*dotvn;
        end  %End of IF (Upwind flux)
        
        %It uses the LLF to define the saturation through edge.
    else
        %Get the max value of characteristic velocity
        %Define a range for the saturtion
        
        Sranglr = [Sleft Sright];
        
        %Get the analitical derivative:
        [dfwdS,dgamadS] = calcdfunctiondS(0,0,Sranglr,1);
        
        %Finally, get the maximun value of saturation:
        alfamax = max(abs(dfwdS*dotvn + dgamadS*dotvg));
        %Denine the numerical flux
        Fleft = fw(1)*dotvn + gama(1)*dotvg;
        Fright = fw(3)*dotvn + gama(3)*dotvg;
        %Define Local Lax-Friedrichs Flux
        
        LLFlux = 0.5*((Fleft + Fright) - alfamax*(Sright - Sleft));
        
        %Calculate the numerical flux through interface using LLF.
        numflux = LLFlux;
        
        earlysw(bedgesize + inedg) = 0.5*(Sleft + Sright);
        
        %Entropy:
        %Get "entrvar":
        %      entrvar = getineqentropy(0.5*(Sleft + Sright),5,2,(1 - entropycond)*paramk);
        %Calculate the entropy flux (right value)
        %      entrflux = entrvar*dotvn;
    end  %End of IF (type of flux)
    
    
    %Obtain the contribution of interface over element to LEFT
    advecterm(leftelem) = advecterm(leftelem) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(rightelem) = advecterm(rightelem) - numflux;
    
    
    
    %     if i==1
    %         satmedio(i,1:2)=[Sright Sleft];
    %     else
    %         satmedio(i,1:2)=[Sleft Sright];
    %     end
    %     axesxx(i,1)=i;
    %Evaluate ENTROPY condition (flux terms)
    %Obtain the contribution of interface over element to LEFT
    %    entrineqterm(leftelem) = ...
    %        entrineqterm(leftelem) + entrflux;   %max(0,sign(dotvn))*(1 - entropycond)*(orderinedgdist(i,1) > 1);  %
    %Obtain the contribution of interface over element to RIGHT
    %    entrineqterm(rightelem) = ...
    %        entrineqterm(rightelem) - entrflux;   %min(0,sign(dotvn))*(1 - entropycond)*(orderinedgdist(i,2) > 1);  %
end  %End of FOR ("inedge")
%plot(axesxx(:,1),satmedio(:,1),'o',axesxx(:,1),satmedio(:,2),'s')
%legend('MAIOR VALOR','MENOR VALOR')
%grid

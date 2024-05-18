% % ** Brittle PFF（AT12）            **
% % ** code by P.M.H @bit.edu.cn (CN) **
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% %  ---------------------------------------
% % Last update: 2024-05-16;
% % Create date: 2024-04-02; 

clear; close all

%%  ***  Reas Ansys Mesh  ***
YourModel = 'LPlate';  % Choose your model

readdir  = ['./ansys_result\',YourModel, '\'];
% read the element and boundary condation 
fprintf(1,'read the mesh\n')
node = load([readdir,'NLIST.DAT']);
sumNode = size(node,1);
elem = load([readdir,'ELIST.DAT']);
fixNode = load([readdir,'fixNode.dat']);

filedir = mkResultsDir('LPlate\');
fdc = [filedir, 'force_displacement.txt']; % as filename tell

%% ***  Material para  *** (Ambati's Paper)
Para.PFModel = 2; % 1-AT2; 2-AT1
Para.ndim = 2; % dim
Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.lambda = 6160; % Lame Constant 
Para.mu = 10950; % Lame Constant
Para.E = Para.mu*(2*Para.mu+3*Para.lambda)/(Para.mu+Para.lambda); % Young's Modulus based on (N/mm2)
Para.nu = Para.lambda/(2*(Para.mu+Para.lambda)); % Poisson's Ratio
Para.Gc = 0.089; % Critical energy release for unstable crack (Gc, N/mm)
Para.Len = 3; % 

Para.NNd = size(node,1); % number of nodes

elem(:,1:2) = [];
node(:,1)   = [];
node = node(:, 1 : Para.ndim);

[GaussInfo] = shapeFunc_valueDeriv(elem, node, Para);
Phi = zeros(Para.NNd,1);

%% Time integration parameters
loadrate = 0.1; % clamped velocity (mm/s)
dt = 1d-2/loadrate; % delta load = loadrate * dt

%% Miehe's staggered scheme  
fdcfid = fopen(fdc,'w');  % force-displacement 
AMtol = 1d-4;
loaddisp = 0;
for inc = 1:1000
    if loaddisp > 0.2
        dt = 1d-3/loadrate;
    end
    loaddisp = loaddisp + dt * loadrate; % quasi-static
    BC = ElasSENT(fixNode, sumNode*2, loaddisp);
    
    AMres = 1; it = 0; 
    while AMres > AMtol
        % compute the disp sub-problem
        [Disp] = assembleElasKK(GaussInfo, elem, Phi, Para, BC);
        
        % compute the phase-field sub-problem
        [Phi] = assembleElasKPhi(GaussInfo, elem, Disp, Para);
        %
        % AMres = norm((Phiold-Phi))/norm(Kphi*Phi); % multi-pass
        AMres = 1d-4; % one-pass
        it = it+1; % iteration counts
    end
    
    % update history & compute internal force
    [GaussInfo, InF] = updateRefEnerg(GaussInfo, elem, Disp, Para);
   
    BDF = sum(InF(BC.BDforce));
    fprintf(fdcfid, ['%6d' repmat('%16.10f ',1,2) '\n'], it, loaddisp, full(BDF) );
    
%% Plot results
if mod(inc-1,10) == 0
    disp(['Exporting, disp: ',num2str(loaddisp),'mm, Load: ', num2str(BDF), 'N'])
    figure(1)
    axis equal;
    PlotContour(node,elem,full(Phi));
    axis off;
end

end

fclose(fdcfid);

LoadForc = textread(fdc);
figure()
plot(LoadForc(:,2),LoadForc(:,3));


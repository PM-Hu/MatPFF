% % ** Brittle PFF（PFCZM）           **
% % ** code by P.M.Hu @bit.edu.cn     ** 
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% %  ---------------------------------------
% % Last update: 2024-12-12
% % Create date: 2024-03-23;

clear; close all

%%  ***  Reas Ansys Mesh  ***
YourModel = 'LPlatePFCZM';  % Choose your model

readdir  = ['./ansys_result\',YourModel, '\'];
% read the element and boundary condation 
node = load([readdir,'NLIST.DAT']);
sumNode = size(node,1);
elem = load([readdir,'ELIST.DAT']);
sumElem = size(elem,1);
fixNode = load([readdir,'fixNode.dat']);
nodeForce = load([readdir,'nodeForce.dat']); %
nodeForce(1,:) = [];

filedir = mkResultsDir('WuLPlate\');
fdc = [filedir, 'force_displacement.txt']; % as filename tell

%% ***  Material para  *** (Wu's Paper, A unified phase-field theory for the mechanics of damage and quasi-brittle failure)
% % Cornelissen soften law is used for L-Plate
Para.PFModel = 3; % 1-AT2; 2-AT1; 3-PFCZM
Para.ndim = 2; % dim
Para.isStress = 1;  % 1 - plane stress, 2 - plane strain
Para.E = 20000; % Young's Modulus based on (N/mm2)
Para.nu = 0.18; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu) * (1-2*Para.nu)); % Lame Constant 
Para.mu = Para.E /(2* (1+Para.nu)); % Lame Constant
Para.Gc = 0.12; % Critical energy release for unstable crack (Gc, N/mm)
Para.Len = 6.25; % 

% for PFCZM
Para.ft = 2.5; % MPa
Para.Lch = Para.E*Para.Gc/Para.ft^2; % Irwin’s internal length

Para.NNd = size(node,1); % number of nodes

elem(:,1:2) = [];
node(:,1)   = [];
node = node(:, 1 : Para.ndim);

[GaussInfo] = shapeFunc_valueDeriv(elem, node, Para);

% initial phi
Phi = zeros(Para.NNd,1);

%% Time integration parameters
loadrate = 0.1; % clamped velocity (mm/s) # quasi-static?
dt = 1d-3/loadrate; % delta load = loadrate * dt

%% Miehe's staggered scheme  
fdcfid = fopen(fdc,'w');  % store force-displacement
AMtol = 1d-4;
loaddisp = 0;
for inc = 1:1000
    % fine dt is required when brutal fracture approaches
    loaddisp = loaddisp + dt * loadrate; % quasi-static
    BC = ElasSENT(fixNode, sumNode*2, loaddisp);
    
    AMres = 1; it = 0; 
    while AMres > AMtol
        % compute the disp sub-problem
        [Disp] = assembleElasKK(GaussInfo, elem, Phi, Para, BC);
        
        % compute the phase-field sub-problem
        Phiold = Phi;
        Phi = NewtonItPhaseField(GaussInfo, elem, Disp, Phi, Para);
        %
%         AMres = norm((Phiold-Phi)); % multi-pass
        AMres = 1d-4; % one-pass
        it = it+1; % iteration counts
    end
    
    % update history & compute internal force
    [GaussInfo, InF] = updateRefEnerg(GaussInfo, elem, Disp, Para);
   
    BDF = sum(InF(BC.BDforce));
    fprintf(fdcfid, ['%6d' repmat('%16.10f ',1,2) '\n'], it, loaddisp, full(BDF) );
    
%% Plot the contour
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


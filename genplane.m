% To compute the best fit plane

close all
%clear
%clc
format long


%% Input data

N = 1000;ntypes = 8;
trialnum = 3;
datafilenum = 2;
lx = 3.5088212900000002e+02;
ly = 3.5088212900000002e+02;
lz = 3.5088212900000002e+02;
masses = zeros(ntypes,1);

%% Analysis

aidarr  = zeros(length(N),3);
rxyzarr = zeros(length(N),3);
fid = fopen(sprintf('./n%d/repmain_singlechain/trial%d_main.data',N,trialnum),'r');

disp('Reading structure')

while(~feof(fid))
    tline = fgetl(fid);
    lvals = strsplit(tline,' ');
    atline = cell2mat(lvals(1));
    if strcmp(atline,'Masses') == 1
        disp('Reading Atom Masses')
        tline = fgetl(fid);
        for i = 1:ntypes
            tline = fgetl(fid);
            lvals = strsplit(tline,' ');
            masses(i,1)  = str2double(cell2mat(lvals(2)));
        end
    end
    if strcmp(atline,'Atoms') == 1
        disp('Reading Atom Data')
        tline = fgetl(fid);
        for i = 1:N
            tline = fgetl(fid);
            lvals = strsplit(tline,' ');
            atomid  = str2double(cell2mat(lvals(1)));
            aidarr(atomid,1) = atomid;
            aidarr(atomid,2) = str2double(cell2mat(lvals(2)));
            aidarr(atomid,3) = str2double(cell2mat(lvals(3)));
            rx = str2double(cell2mat(lvals(4)));
            ry = str2double(cell2mat(lvals(5)));
            rz = str2double(cell2mat(lvals(6)));
            ix = str2double(cell2mat(lvals(7)));
            iy = str2double(cell2mat(lvals(8)));
            iz = str2double(cell2mat(lvals(9)));
            rxyzarr(atomid,1) = rx + ix*lx;
            rxyzarr(atomid,2) = ry + iy*ly;
            rxyzarr(atomid,3) = rz + iz*lz;
        end
        break;
    end
end

%masses(:,1) = 1.0;
disp('Removing COM')
rxcm = 0.0; rycm = 0.0; rzcm = 0.0; totmass = 0.0;
for i = 1:N
    rxcm = rxcm + masses(aidarr(atomid,3),1)*rxyzarr(atomid,1);
    rycm = rycm + masses(aidarr(atomid,3),1)*rxyzarr(atomid,2);
    rzcm = rzcm + masses(aidarr(atomid,3),1)*rxyzarr(atomid,3);
    totmass = totmass + masses(aidarr(atomid,3),1);
end

rxcm = rxcm/totmass; rycm = rycm/totmass; rzcm = rzcm/totmass;


%% Copy arrays before subtracting for plotting

rxyzplot = rxyzarr;

rxyzarr(:,1) = rxyzarr(:,1) - rxcm;
rxyzarr(:,2) = rxyzarr(:,2) - rycm;
rxyzarr(:,3) = rxyzarr(:,3) - rzcm;

disp('Fitting Plane')

% Augment arrays

rxyzarr(:,4) = ones(N,1);

% SVD

[u, s, v] = svd(rxyzarr, 0);
P = v(:,4);


% Normalized Values of P

fid = fopen(sprintf('normvec_N%d_Trial%d.txt',N,trialnum),'w');
fprintf(fid,'%s\t%s\t%s\n','Px','Py','Pz');
fprintf(fid,'%g\t%g\t%g\n',P(1),P(2),P(3));

pnorm = sqrt(P(1)^2+P(2)^2+P(3)^2);
pnx = abs(P(1))/pnorm; pny = abs(P(2))/pnorm; pnz = abs(P(3))/pnorm;



% Generate Points on the Plane

checkpoints(:,1) = rxyzplot(:,1) + pnx*100;
checkpoints(:,2) = rxyzplot(:,2) + pny*100;
checkpoints(:,3) = rxyzplot(:,3) + pnz*100;



disp('Plotting Data')
h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$x$','FontSize',20,'Interpreter','Latex')
ylabel('$y$','FontSize',20,'Interpreter','Latex')
zlabel('$z$','FontSize',20,'Interpreter','Latex')
plot3(rxyzplot(:,1),rxyzplot(:,2),rxyzplot(:,3),'rx','MarkerSize',8,'MarkerFaceColor','r')
plot3(checkpoints(:,1), checkpoints(:,2), checkpoints(:,3),'kx','MarkerSize',8)
view(-35,45)
legendinfo{1} = 'Initial Datapoints';
legendinfo{2} = 'Check Points to Plane Fit';
legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
saveas(h1,sprintf('fitdata_N%d_Trial%d.png',N,trialnum))


fprintf(fid,'%s\t%s\t%s\n','Pxnorm','Pynorm','Pznorm');
fprintf(fid,'%g\t%g\t%g\n',pnx,pny,pnz);
fprintf('%g\t%g\t%g\n',pnx,pny,pnz);





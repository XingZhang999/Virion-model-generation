% Virus sample
% after search 2D, need to fit the search results to a 3D model
% the CCG image will be used to fitting the 3D model !
% asign full particles from 3D model,
% and Z value will be asigned from 3D model useing closest point 

% step 2
% read the fitted virus model and use it to asign particles 
% to Virus, used to genarate a closer discription of NM protein array
% not finished !!!

% Xing Zhang @20220628

clear;clc
addpath(genpath('/ssd/Software/emClarity/emClarity_Aug26_noPTmask_Bfactor_smoothXY_new'))
tic;

apix = 0.87;
layerfix = 82; % pixel unit
% samplingRate = 1 % for mod and star file
handness = -1;
dfoffset = 0; %300 % for lacey carben grid, the estimated defocus maybe not the same as virus, so need to addjust
ifplot = 0; % or 4
ifcheck = 1;
% searchbin = 2

files=struct2cell(dir('./cryosparc_DW_*_pt*_srh.star'));

%
for i=1  %:size(files,2)
    tic;
    name=char(files(1,i));
    fprintf('\n*** Processing on %s ***\n',name);
    idend=strfind(name,'_srh'); %_rem
    namenoext=name(1:idend-1);
    remfile=sprintf('%s_rem.txt',namenoext);
    remfilefit=sprintf('%s_remfit.txt',namenoext);
    fitmrc=sprintf('%s_fit.mrc',namenoext);
    starfile_1img=sprintf('%s_fit1.star',namenoext);
    realfolder=sprintf('%s_regular',namenoext);
    
    if( exist(realfolder,'dir')~=7 )
    %fprintf('Making new scratch directory... [%s]\n',scratchDir);
    mkdir(realfolder);
    end

    [b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12] = textread(starfile_1img,'%s %d %d %f %f %f %f %f %f %f %d %d');
    
    %readimg=sprintf('%s/search_valmax.mrc',dirname)
    %CCGmap=single(getVolume(MRCImage(readimg)));
    %CCGmap = imresize(CCGmap,size(CCGmap).*searchbin );
    num=size(b2,1);
    virusrem = load(remfile);
    virusfit = load(remfilefit);
    
    if exist(remfile,'file')~=2 || exist(remfilefit,'file')~=2 || num <=3
        continue
    end
    
    count=num;
    poslist = zeros(num,8);
    meandf = mean((b4+b5)/2);
    % X Y coords b2 b3
    poslist(:,1) = b2; 
    poslist(:,2) = b3;
    % defocus to hight
    poslist(:,3) = ((b4+b5)/2-meandf)/apix;
    % could be used to 
    % poslist(:,4) = 0;
    % eulor
    poslist(:,5) = b7;
    poslist(:,6) = b8;
    poslist(:,7) = b9;
    
    centerrem = virusrem(8:9)-[600 600+(5760-4092)/2]; % for pad sqare, add 600 
    height = virusrem(15);
    layeridst = virusrem(16)-5; 
    poslist(:,1:2) = poslist(:,1:2) - ones(count,1)*centerrem;
    
        if ifplot >= 1 
        figure; hold on;grid on; box on;axis equal;
        plot3(poslist(:,1), poslist(:,2), poslist(:,2)*0,'bo');
        %plot3(newtable(:,11), newtable(:,12), ptCloudtomo.Location(:,3),'r.');
        end
        
    % average from vector should be better than from eulor angle
        Xaxis = [1,0,0];
        Yaxis = [0,1,0];
        Zaxis = [0,0,1];
        Xvector = zeros(num,3);
        Yvector = zeros(num,3);
        Zvector = zeros(num,3);
    for j=1:num
        Rspa = BH_defineMatrix(-[poslist(j,7) poslist(j,6) poslist(j,5)], 'SPIDER', 'inv');
        Rspa = reshape(Rspa,3,3);
        Xvector(j,1:3)=Rspa*Xaxis';
            if ifplot == 1 
               quiver3(poslist(j,1),poslist(j,2),poslist(j,2)*0,Xvector(j,1)*100,Xvector(j,2)*100,Xvector(j,3)*100,'d');
            end
        Yvector(j,1:3)=Rspa*Yaxis';
            if ifplot == 2 
               quiver3(poslist(j,1),poslist(j,2),poslist(j,2)*0,Yvector(j,1)*100,Yvector(j,2)*100,Yvector(j,3)*100,'d');
            end
        Zvector(j,1:3)=Rspa*Zaxis';
            if ifplot == 3 
               quiver3(poslist(j,1),poslist(j,2),poslist(j,2)*0,Zvector(j,1)*100,Zvector(j,2)*100,Zvector(j,3)*100,'d');
            end
    end
    
    Virusvector=sum(Zvector)./num;
    
    if Virusvector(1)>=0
    theta=acosd(Virusvector(3));
    psi=acosd(Virusvector(2)/sqrt(Virusvector(1)^2+Virusvector(2)^2));
    end
    if Virusvector(1)<0
    theta=180+acosd(-Virusvector(3));
    psi=acosd(-Virusvector(2)/sqrt(Virusvector(1)^2+Virusvector(2)^2));
    end
    Rz=[cosd(-psi),-sind(-psi),0;...
        sind(-psi),cosd(-psi),0;...
        0,0,1];
    Rx=[1,0,0; ...
        0,cosd(-theta),-sind(-theta);...
        0,sind(-theta),cosd(-theta)];
    
    a = virusfit(4);
    b = virusfit(5);
    bestshift = virusfit(6:7);
    
    % generate a cylinder model for filtering particles
    rounds = height/layeridst/2*1.2;
    t = -rounds*360:0.5:rounds*360; % try to change it to discrete number
    x3=a*cosd(t);
    y3=b*sind(t);
    z3=layerfix/*t/360;
    
    rot3D=zeros(size(t,2),3);
    
    % write out mod file to view
    modelpos = sprintf('./%s/%s_regu3D.pos',realfolder,namenoext);
    modelmod = sprintf('./%s/%s_regu3D.mod',realfolder,namenoext);
    
    fileID0 = fopen(modelpos,'w');
    for ii=1:size(t,2)
                rot3D(ii,:)=(Rz*Rx*[x3(ii),y3(ii),z3(ii)]')' +[bestshift 0]; 
                fprintf(fileID0,'%f %f %f\n',rot3D(ii,1:3));
    end
    fclose(fileID0);
    system(sprintf('point2model -number 1 -sphere 5 -scat %s %s', modelpos, modelmod));
    
    % caculate Z and filtering particles
    % write out star file of filtered particles
    % write out mod file to view
    %pos_in = sprintf('./%s/%s_old.pos',realfolder,namenoext);
    %mod_in = sprintf('./%s/%s_old.mod',realfolder,namenoext);

    pos_out = sprintf('./%s/%s_regu2D.pos',realfolder,namenoext);
    mod_out = sprintf('./%s/%s_regu2D.mod',realfolder,namenoext);

    %particlestar0 = sprintf('%s_fit0.star',namenoext);
    %particlestar1 = sprintf('%s_fit1.star',namenoext); % fit1 is better
    %particlestar2 = sprintf('%s_fit2.star',namenoext);
    regulatestar = sprintf('%s_regu.star',namenoext);
    
    %fileID2 = fopen(particlestar1,'w'); % fit2 seems in right handness, so close  
    %fileID3 = fopen(particlestar2,'w'); 
    
    %fileID4 = fopen(pos_in,'w');
    fileID5 = fopen(pos_out,'w');
    fileID6 = fopen(regulatestar,'w');

        if ifplot == 4 
        figure; hold on;grid on; box on;axis equal;
        plot3(rot3D(:,1),rot3D(:,2),rot3D(:,3),'.');
        end
        
        used = 0;
    
    for ii=1:count
        
        deltaZestimate = b*Xvector(ii,3)+poslist(ii,1:2)*Virusvector(1:2)'/norm(Virusvector(1:2))*Virusvector(3); %sin(angletoplan);
        particlexyz = [poslist(ii,1:2) deltaZestimate];
        
        modeldist = sqrt( sum((rot3D-particlexyz).^2,2) );
        %modeldist = sqrt( sum((rot3D(:,1:2)-poslist(ii,1:2)).^2,2) );
        [mindist,minID] = min(modeldist);
        deltaZ = rot3D(minID,3); % not wright !!! one particle could belong to different classes
        
        
        particlexyz = [poslist(ii,1:2) deltaZ];
        
            if ifplot == 4 
            %plot3(poslist(ii,1),poslist(ii,2),deltaZestimate,'bo');
            quiver3(poslist(ii,1),poslist(ii,2),deltaZestimate,Xvector(ii,1)*100,Xvector(ii,2)*100,Xvector(ii,3)*100,'o');
            end
            
        dist2axis = norm(cross([poslist(ii,1:2) 0],[Virusvector(1:2) 0]));
        if mindist <= 50 % && dist2axis < 0.7*a && dist2axis > 0.3*a % A
        fprintf(fileID5,'%f %f %f\n',poslist(ii,1:2),deltaZ);
        
        dist2axis;
        angle = acosd(dist2axis/a)
        b = abs(Xvector(ii,3))/sqrt(Xvector(ii,1)^2+Xvector(ii,2)^2)*a*cosd(angle)./sind(angle)
        
        deltaDF = deltaZ*handness*apix+dfoffset;
        %fprintf(fileID2,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii)+deltaDF,b5(ii)+deltaDF,b6(ii),b7(ii),b8(ii),b9(ii),b10(ii),b11(ii),b12(ii)); % write virus ID as class ID
        %fprintf(fileID3,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii)-deltaDF,b5(ii)-deltaDF,b6(ii),b7(ii),b8(ii),b9(ii),b10(ii),b11(ii),b12(ii)); % write virus ID as class ID
        used = used+1;
        end
    end  
        
        fclose(fileID5);
        
    system(sprintf('point2model -number 1 -sphere 20 -scat ./%s ./%s', pos_out, mod_out));
    
    toc;

end

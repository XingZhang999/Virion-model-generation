% run under MATLAB2023a
% Virus sample
% After template matching 2D
% Try to ganerate Virus model for each virion
% This 3D model is used to exclude false picking and asign full sub-particles 
% The Z values will be asigned from 3D model from closest point 
% By Xing Zhang @Tsinghua University, 20240407

clear;clc

% input parameters
apix = 0.87 % the pixel size of micrographs
layer_inter = 80.5 % the average layer distance in pixel unit
dfoffset = 0 % the estimated defocus maybe not the same as virus, so need to addjust
ifplot = 2 % if do plot the intermedia results for check

% input particles as a format of RELION star file
files=struct2cell(dir('./cryosparc_DW_*_pt*_srh.star'));

% run for all the particle files
%for i= 1 % :size(files,2)
i=1
    
    tic;
    name=char(files(1,i));
    fprintf('\n*** Processing on %s ***\n',name);
    idend=strfind(name,'_srh'); 
    namenoext=name(1:idend-1);
    
    remfile=sprintf('%s_rem.txt',namenoext); % try to skip rem file later
    starfile_1img=sprintf('%s_srh.star',namenoext);
    virusrem = load(remfile);
    [b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12] = textread(starfile_1img,'%s %d %d %f %f %f %f %f %f %f %d %d');
    
    % output folder and files
        realfolder=sprintf('%s_regular',namenoext);
        if( exist(realfolder,'dir')~=7 )
            mkdir(realfolder);
        end
    
        particlepos = sprintf('./%s/%s_particle.pos',realfolder,namenoext);
        particlemod = sprintf('./%s/%s_particle.mod',realfolder,namenoext);
        modelpos = sprintf('./%s/%s_model.pos',realfolder,namenoext);
        modelmod = sprintf('./%s/%s_model.mod',realfolder,namenoext);
    
        particlestar0 = sprintf('%s_fit0.star',namenoext);
        particlestar1 = sprintf('%s_fit1.star',namenoext);
        particlestar2 = sprintf('%s_fit2.star',namenoext);
        remfilefit=sprintf('%s_remfit.txt',namenoext);
        
    % creat a array for save points
        num=size(b2,1);
        poslist = zeros(num,12);
        % X Y coords
        poslist(:,1) = b2; 
        poslist(:,2) = b3;
        % defocus value
        poslist(:,3) = b4;
        poslist(:,4) = b5;
        % Eulor angles: rot, tilt, psi
        poslist(:,5) = b7;
        poslist(:,6) = b8;
        poslist(:,7) = b9;
    
    % estimate virus vector, either use mean value or peak value
        if 0
        % angletilt = mean(poslist(:,6))
        [aa,ss] = hist(poslist(:,6));
        [~,maxid] = max(aa);
        angletilt = ss(maxid);
        % anglepsi = mean(poslist(:,7))
        [aa,ss] = hist(poslist(:,7));
        [~,maxid] = max(aa);
        anglepsi = ss(maxid);
        % vector from high probability eulor angle
    	R = BH_defineMatrix(-[anglepsi angletilt 0], 'SPIDER', 'inv');
    	R = reshape(R,3,3);
        Virusvector=(R*Zaxis')';
        end

    % read virus paramters from rem.txt
        centerrem = virusrem(8:9)-[600 600+(5760-4092)/2]; % for padded micrography, added 600 pixel
        Virusvector = virusrem(11:12);
        width = virusrem(14);
        height = virusrem(15);
        layer_inter2D = virusrem(16)-5; % check why -5 later
        poslist(:,1:2) = poslist(:,1:2) - ones(num,1)*centerrem; 
        radius = width / 2 - 70;
        
    % paramters can be read from rem.txt or caculate from particles.star
        Xaxis = [1,0,0];
        Yaxis = [0,1,0];
        Zaxis = [0,0,1];
    
    %  vector from average of particle Zvector, should be better than from eulor angle
        Xvector = zeros(num,3);
        Yvector = zeros(num,3);
        Zvector = zeros(num,3);
        
        if ifplot >= 4
        figure; hold on;grid on; box on;axis equal;
        plot(poslist(:,1), poslist(:,2), 'b.');
        end
        
        % skip the points have wrong Z direction
        poslist(:,8)=1;
    for j=1:num
        
        Rspa = BH_defineMatrix(-[poslist(j,7) poslist(j,6) poslist(j,5)], 'SPIDER', 'inv');
        Rspa = reshape(Rspa,3,3);
        Zdirction = Rspa*Zaxis';
        Zvector(j,1:3)=Zdirction;
        Xvector(j,1:3)=Rspa*Xaxis';
        Yvector(j,1:3)=Rspa*Yaxis';

        if  acosd(dot(Zdirction(1:2),Virusvector(1:2))) > 12 
            poslist(j,8)=-1;            
        end

        if norm(cross([Virusvector 0],[poslist(j,1:2) 0])) > radius
        	poslist(j,8)=-1;            
        end
        
        if norm(dot([Virusvector 0],[poslist(j,1:2) 0])) > height/2
        	poslist(j,8)=-1;            
        end
    end
    
        Avevector=sum(Zvector(poslist(:,8)>0,:))./sum(poslist(:,8)>0);
        Vectorchange=acosd(dot(Avevector(1:2),Virusvector(1:2)));
        
        if  Vectorchange > 20
            fprintf('Virus vector not match ! skip !\n');
            %continue
        else
            fprintf('Virus vector change %.1f\n',Vectorchange);
        end

        if ifplot >= 2
        figure; hold on;grid on; box on;axis equal;
        plot(poslist(poslist(:,8)>0,1), poslist(poslist(:,8)>0,2),'r.');
        quiver(poslist(poslist(:,8)>0,1), poslist(poslist(:,8)>0,2),Zvector(poslist(:,8)>0,1)*100,Zvector(poslist(:,8)>0,2)*100,'r');
        end
        
       %
    Virusvector = [Virusvector Avevector(3)];
    Virusvector = Virusvector./norm(Virusvector);
    
    if Virusvector(1)>=0
        theta=acosd(Virusvector(3));
        psi=acosd(Virusvector(2)/sqrt(Virusvector(1)^2+Virusvector(2)^2));
    else
        theta=180+acosd(-Virusvector(3));
        psi=acosd(-Virusvector(2)/sqrt(Virusvector(1)^2+Virusvector(2)^2));
    end
    
    Rz=[cosd(-psi),-sind(-psi),0;...
        sind(-psi),cosd(-psi),0;...
        0,0,1];
    Rx=[1,0,0; ...
        0,cosd(-theta),-sind(-theta);...
        0,sind(-theta),cosd(-theta)];
    
    % fitting may genarate acurate a and rough b, but b is still necessary 
    starta = radius;
    [am,bm,anglesearch] = meshgrid(0.86:0.02:1.00,0.6:0.05:1,0:1:0);
    
    ratiolist = [am(:) bm(:) anglesearch(:)];
    rationum = size(ratiolist,1);
    fprintf('Samplling number %d\n',rationum);
    score = zeros(rationum,1);
    shiftlist = zeros(rationum,16,'double');
         
    for r=1:rationum
    
        fprintf('--- run %d, ratio %.2f - %.2f, Xangle %.1f, ',r,ratiolist(r,1),ratiolist(r,2),ratiolist(r,3));
    ratio=ratiolist(r,1:2);
    
    a = starta*ratiolist(r,1);
    b = a*ratiolist(r,2);
    xangle = ratiolist(r,3);
    Rxangle = [1,0,0; ...
        0,cosd(-xangle),-sind(-xangle);...
        0,sind(-xangle),cosd(-xangle)];

    rounds = height/layer_inter2D/2*1;
    t = -rounds*360:0.5:rounds*360;
    modelsize = size(t,1);

    x3=a*cosd(t);
    y3=b*sind(t);
    z3=layer_inter*t/360;
    
    dx3=b*cosd(t);
    dy3=a*sind(t);
    dx2=-a*sind(t);
    dy2=b*cosd(t);

    rot3D=zeros(size(t,2),3);
    Back3Dv=zeros(size(t,2),3);
    Back3Ds=zeros(size(t,2),3);
    
    % rotated 3D model
    for ii=1:size(t,2)
                rot3D(ii,:)=(Rz*Rx*Rxangle*[x3(ii),y3(ii),z3(ii)]')'; % + [bestshift 0]; %+[centerrem 0];
                
                Back3Dv(ii,:)=(Rz*Rx*Rxangle*[dx3(ii),dy3(ii),0]')';
                Back3Dv(ii,:)=Back3Dv(ii,:)./norm(Back3Dv(ii,:));
                
                Back3Ds(ii,:)=(Rz*Rx*Rxangle*[dx2(ii),dy2(ii),0]')';
                Back3Ds(ii,:)=Back3Ds(ii,:)./norm(Back3Ds(ii,:));
    end
                     
        % project the 3D model to 2D, and register them to 2D detections
        % point cloud regist
        ptCloud1 = pointCloud([poslist(poslist(:,8)>0,1:2) poslist(poslist(:,8)>0,3).*0]);
        ptCloud2 = pointCloud([rot3D(:,1:2) rot3D(:,3).*0]);
        tform = pcregistericp(ptCloud2,ptCloud1);
        %tform = pcregisterndt(ptCloud2,ptCloud1,0.5);
        
        %ptCloud1 = pointCloud([poslist(poslist(:,8)>0,1:2)]);
        %ptCloud2 = pointCloud([rot3D(:,1:2) rot3D(:,3)]);
        %tform = pcregistericp_vmg(ptCloud2,ptCloud1);

        %ptCloud1 = pointCloud([poslist(poslist(:,8)>0,1:2)]);
        %ptCloud2 = pointCloud(rot3D(:,1:2));
        %tform = fitgeotrans(rot3D(:,1:2),poslist(poslist(:,8)>0,1:2),'simularity');
        %[tform,movingreg,rmse] = pcregistercpd(ptCloud2,ptCloud1,'Transform','Rigid');
        
        %disp(tform.T)
        ptCloudTformed = pctransform(ptCloud2,tform);
        
        %
        transferm = tform.T';
        shiftlist(r,:) = transferm(:)'; % to save the shifts
        
        % find the clost point
        diff_record1 = [];
        for ii=1:modelsize
            modeldist = sqrt( sum((ptCloudTformed.Location(ii,1:2)-poslist(:,1:2)).^2,2) ); % distance in 2D
            [mindist,minID] = min(modeldist);
            
            diff1 = norm(poslist(minID,1:2)-ptCloudTformed.Location(ii,1:2));
            diffang1 = acosd(dot(Back3Dv(ii,1:3),Xvector(minID,1:3)));

            diff_record1 = [diff_record1; diff1]; %*diffang1];
        end
        
        if ifplot >= 4 && mod(r,20)==1
            figure; hold on; axis equal;
            plot(ptCloudTformed.Location(:,1),ptCloudTformed.Location(:,2),'r.');
            plot(poslist(poslist(:,8)>0,1),poslist(poslist(:,8)>0,2),'go');
        end

        % need to caculate the nearest point in shifts and eulor angles
        diff_record2 = [];
        
        for ii=find(poslist(:,8)>0')'
            modeldist = sqrt( sum((ptCloudTformed.Location(:,1:2)-poslist(ii,1:2)).^2,2) ); % distance in 2D
            [mindist,minID] = min(modeldist);
            poslist(ii,10) = minID;
            diff2 = norm(poslist(ii,1:2)-ptCloudTformed.Location(minID,1:2));
            diffang2 = acosd(dot(Back3Dv(minID,1:3),Xvector(ii,1:3)));

            %if diffang2 < 5
            %    diffang2 = 5;
            %end
            %if diffang2 > 40
            %    diffang2 = 180;
                %continue
            %end

            diff_record2 = [diff_record2; diff2]; %*diffang2]; %diff1+diff2]; 
        
        end

        score(r) = mean(diff_record1)*mean(diff_record2);
        fprintf('mean diff %.2f \n',score(r));

    end
   %
    [~,minID]=min(score(:));
    a=starta*ratiolist(minID,1)
    b=a*ratiolist(minID,2)
    bestscore = score(minID)
    bestshift = shiftlist(minID,:)
    xangle = ratiolist(r,3);
    Rxangle = [1,0,0; ...
        0,cosd(-xangle),-sind(-xangle);...
        0,sind(-xangle),cosd(-xangle)];

    fprintf('the best ratio %.2f - %.2f, Xangle %.1f, diff %.2f \n',ratiolist(minID,1),ratiolist(minID,2),ratiolist(minID,3),score(minID));
%
        if ifplot >= 2
            figure; hold on; 
            plot3(ratiolist(:,1),ratiolist(:,2),score(:),'.');
            plot3(ratiolist(minID,1),ratiolist(minID,2),score(minID),'o')
        end
    %
    % plot the best fit
    %rounds = height/layer_inter2D/2*0.8;
    %t = -rounds*360:1:rounds*360;
    
    x3=a*cosd(t);
    y3=b*sind(t);
    z3=layer_inter*t/360;

        for ii=1:size(t,2)
                rot3D(ii,:)=(Rz*Rx*Rxangle*[x3(ii),y3(ii),z3(ii)]')'; % + [bestshift 0]; %+[centerrem 0];
                
                %Back3Dv(ii,:)=(Rz*Rx*[dx3(ii),dy3(ii),0]')';
                %Back3Dv(ii,:)=Back3Dv(ii,:)./norm(Back3Dv(ii,:));
                
                %Back3Ds(ii,:)=(Rz*Rx*[dx2(ii),dy2(ii),0]')';
                %Back3Ds(ii,:)=Back3Ds(ii,:)./norm(Back3Ds(ii,:));
        end
        ptCloud2 = pointCloud(rot3D);
        besttform = rigid3d([bestshift(1:3);bestshift(5:7);bestshift(9:11)],bestshift(13:15));
        ptCloudTformed = pctransform(ptCloud2,besttform);

        if ifplot >= 1 
            figure; hold on; axis equal;
            plot(ptCloudTformed.Location(:,1),ptCloudTformed.Location(:,2),'r.');
            plot(poslist(poslist(:,8)>0,1),poslist(poslist(:,8)>0,2),'go');
        end

    %%
    % use particle pointed to the axis for futher processing
    
    for ii=1:num
        
        % re check the functions !!!
        deltaZestimate = b*Xvector(ii,3)+poslist(ii,1:2)*Virusvector(1:2)'/norm(Virusvector(1:2))*Virusvector(3); 
        dist2axis2D = norm(cross([poslist(ii,1:2) 0],[Virusvector(1:2) 0]/norm([Virusvector(1:2) 0])));
        dist2axis3D = norm(cross([poslist(ii,1:2) deltaZestimate],Virusvector(1:3)));

        pjaxis = dot([poslist(ii,1:2) deltaZestimate],Virusvector)*Virusvector;
        angledot = dot(([poslist(ii,1:2) deltaZestimate]-pjaxis)/norm([poslist(ii,1:2) deltaZestimate]-pjaxis),Xvector(ii,1:3));
        
        poslist(ii,9)=dist2axis2D;
        poslist(ii,10)=dist2axis3D;
        poslist(ii,11)=angledot;
        
        if dist2axis2D > a*0.75 || dist2axis2D < a*0.25 || angledot < 0.707 
            poslist(ii,8)=-1;
        end
    end

    fitlist=find(poslist(:,8)>=0);
    
    % a aecond round fitting for b, using normals
    bratiolist = [0.8:0.01:1.2]';
    brationum = size(bratiolist,1);
    bscore = zeros(brationum,1);
    
    for jj=1:brationum
        
    b = a*bratiolist(jj);
    rounds = height/layer_inter2D/2*1.2;
    t = -rounds*360:0.5:rounds*360;
    x3=a*cosd(t);
    y3=b*sind(t);
    z3=layer_inter*t/360;
    
    dx3=b*cosd(t);
    dy3=a*sind(t);
    dx2=-a*sind(t);
    dy2=b*cosd(t);

    rot3D=zeros(size(t,2),3);
    Back3Dv=zeros(size(t,2),3);
    Back3Ds=zeros(size(t,2),3);
    
    for ii=1:size(t,2)
                rot3D(ii,:)=(Rz*Rx*[x3(ii),y3(ii),z3(ii)]')' + [bestshift 0]; %+[centerrem 0];
                
                Back3Dv(ii,:)=(Rz*Rx*[dx3(ii),dy3(ii),0]')';
                Back3Dv(ii,:)=Back3Dv(ii,:)./norm(Back3Dv(ii,:));
                
                Back3Ds(ii,:)=(Rz*Rx*[dx2(ii),dy2(ii),0]')';
                Back3Ds(ii,:)=Back3Ds(ii,:)./norm(Back3Ds(ii,:));
    end
     
        diff_record2 = [];
        
    for ii=fitlist'
        
        deltaZestimate = b*Xvector(ii,3)+poslist(ii,1:2)*Virusvector(1:2)'/norm(Virusvector(1:2))*Virusvector(3); 
        %particlexyz = [poslist(ii,1:2) deltaZestimate];
        
        %modeldist = sqrt( sum((rot3D-particlexyz).^2,2) ); % distance in 3D
        modeldist = sqrt( sum((rot3D(:,1:2)-poslist(ii,1:2)).^2,2) )+1000*(-sign(deltaZestimate*rot3D(:,3))+1); % distance in 2D
        
        [mindist,minID] = min(modeldist);
        %deltaZ = rot3D(minID,3); % not wright !!! one particle could belong to different classes
        
        %particlexyz = [poslist(ii,1:2) deltaZ];
        
        if mindist <= 50 % A
        diff = dot(Back3Dv(minID,1:3),Xvector(ii,1:3));
        diff_record2 = [diff_record2 diff]; % sqrt(sum(diff(:).^2));
        end
    end
    
    bscore(jj) = mean(diff_record2);

    end
    
    % find out the best match of Xvector
    [~,minIDb]=min(bscore(:));
    bestratio2=bratiolist(minIDb,1)

    %b=a*0.98 ; % trust a high probability ratio rather than fitted ratio
        
    if ifplot == 3
        figure; 
        plot3(ratiolist(:,1),ratiolist(:,2),score(:),'.');
        figure; 
        plot(bratiolist(:,1),bscore(:),'o');
    end
    
    if ifplot == 4 
    figure; hold on;grid on; box on;axis equal;
    end
    
    % generate a cylinder model for filtering particles
    rounds = height/layer_inter2D/2*1.1;
    t = -rounds*360:5:rounds*360;
    x3=a*cosd(t);
    y3=b*sind(t);
    z3=layer_inter*t/360; 
    
    dx3=b*cosd(t);
    dy3=a*sind(t);
    dx2=-a*sind(t);
    dy2=b*cosd(t);

    rot3D=zeros(size(t,2),3);
    Back3Dv=zeros(size(t,2),3);
    Back3Ds=zeros(size(t,2),3);
    
    fileIDfit = fopen(remfilefit,'w');
    fprintf(fileIDfit,'%f %f %f %f %f %f %f \n',starta,bestratio1,bestratio2,a,b,bestshift);
    fclose(fileIDfit);
    
    fileID0 = fopen(modelpos,'w');
    for ii=1:size(t,2)
                rot3D(ii,:)=(Rz*Rx*[x3(ii),y3(ii),z3(ii)]')'+[bestshift 0]; %+[centerrem 0];
                
                Back3Dv(ii,:)=Rz*Rx*[dx3(ii),dy3(ii),0]';
                Back3Dv(ii,:)=Back3Dv(ii,:)./norm(Back3Dv(ii,:));
               
                Back3Ds(ii,:)=Rz*Rx*[dx2(ii),dy2(ii),0]';
                Back3Ds(ii,:)=Back3Ds(ii,:)./norm(Back3Ds(ii,:));
                fprintf(fileID0,'%f %f %f\n',rot3D(ii,1:3));
    end
    fclose(fileID0);
    system(sprintf('point2model -number 1 -sphere 5 -scat %s %s', modelpos, modelmod));
    
    if ifplot == 4 
        for ii=1:size(t,2)
            %plot3(poslist(ii,1),poslist(ii,2),deltaZestimate,'bo');
            quiver3(rot3D(ii,1),rot3D(ii,2),rot3D(ii,3),Back3Dv(ii,1)*100,Back3Dv(ii,2)*100,Back3Dv(ii,3)*100,'r');
        end
    end
    
    %% caculate Z and filtering particles
	particles_used = 0;

	fileID1 = fopen(particlepos,'w');
    fileID2 = fopen(particlestar1,'w'); % fit1 seems is the right handness
    fileID3 = fopen(particlestar2,'w');
    fileID4 = fopen(particlestar0,'w'); % without defocus correction

    for ii=1:num
        deltaZestimate = b*Xvector(ii,3)+poslist(ii,1:2)*Virusvector(1:2)'/norm(Virusvector(1:2))*Virusvector(3);

        particlexyz = [poslist(ii,1:2) deltaZestimate];
        
        modeldist = sqrt( sum((rot3D-particlexyz).^2,2) );
        [mindist,minID] = min(modeldist);
        deltaZ = rot3D(minID,3); % not wright !!! one particle could belong to different classes
        
        particlexyz = [poslist(ii,1:2) deltaZ];
        
        if mindist <= 50 % A
            fprintf(fileID1,'%f %f %f\n',poslist(ii,1:2),deltaZ);
        
            if ifplot == 4 && sum(ismember(fitlist',ii))
                %plot3(poslist(ii,1),poslist(ii,2),deltaZestimate,'bo');
                quiver3(poslist(ii,1),poslist(ii,2),deltaZestimate,Xvector(ii,1)*100,Xvector(ii,2)*100,Xvector(ii,3)*100,'og');
            end
        
        deltaDF = deltaZ*(-1)*apix+dfoffset;
        fprintf(fileID2,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii)+deltaDF,b5(ii)+deltaDF,b6(ii),b7(ii),b8(ii),b9(ii),b10(ii),b11(ii),b12(ii)); % write virus ID as class ID
        fprintf(fileID3,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii)-deltaDF,b5(ii)-deltaDF,b6(ii),b7(ii),b8(ii),b9(ii),b10(ii),b11(ii),b12(ii)); % write virus ID as class ID
        fprintf(fileID4,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii),b5(ii),b6(ii),b7(ii),b8(ii),b9(ii),b10(ii),b11(ii),b12(ii)); % write virus ID as class ID
        particles_used = particles_used+1;
        end
    end
    
    fclose(fileID1);
    fclose(fileID2);
    fclose(fileID3);
    fclose(fileID4);
    
    % particles of each Virus in 3D
    system(sprintf('point2model -number 1 -sphere 5 -scat %s %s', particlepos, particlemod));
    fprintf('--- %d in, %d fit, %d out, %.1f precenter --- \n ',count,size(fitlist',2),particles_used,particles_used/count*100);
    toc;
    
%end

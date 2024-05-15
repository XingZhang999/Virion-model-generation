% Step 2
% After template matching,
% try to ganerate model for each target,
% writeout a image file model.png for checking the fitting,
% and a text file model.log to save the paramters.

% By Xing Zhang @20240106

clear;clc

% input particles as a format of RELION star file
files=struct2cell(dir('./cryosparc_DW_*_pt*_srh.star'));

% input parameters
apix = 0.87; % the pixel size of micrographs
layer_inter = 81; % the average layer distance in pixel unit

% optional input 
ifplot = 0; % whether do plot the intermedia results for check
            % 0 or 1
            % 2 for check the ratio space 
            % 3 for check the normals 
            % 4 for check the fitting for each sampling ratio 

% run for all the particle files
for i= 1:size(files,2)
    
    tic;
    name=char(files(1,i));
    idend=strfind(name,'_srh'); 
    namenoext=name(1:idend-1);
    fprintf('\n*** Processing on %s ***\n',namenoext);

    remfile=sprintf('%s_rem.txt',namenoext); 
    starfile_1img=sprintf('%s_srh.star',namenoext);

    virusrem = load(remfile);
    [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12] = textread(starfile_1img,'%s %d %d %f %f %f %f %f %f %f %d %d'); %#ok<DTXTRD>
    
    % output files
        modelpng = sprintf('%s_model.png',namenoext);
        modelparam = sprintf('%s_model.log',namenoext);
        
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

    % read virus paramters from rem.txt
        centerrem = virusrem(8:9)-[600 600+(5760-4092)/2]; % for padded micrography, added 600 pixel
        Virusvector = [virusrem(11:12) 0];
        width = virusrem(14);
        height = virusrem(15);
        layer_inter2D = virusrem(16); 
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
        
        if ifplot >= 1
            ff=figure('Position',[680   558   1400   420]); 
            tiledlayout(1,3); nexttile; hold on;grid on; box on;axis equal;
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

        if norm(cross(Virusvector,[poslist(j,1:2) 0])) > radius
        	poslist(j,8)=-1;            
        end
        
        if norm(dot(Virusvector,[poslist(j,1:2) 0])) > height/2
        	poslist(j,8)=-1;            
        end
    end
    
    poslist = poslist(poslist(:,8)>0,:);
    Xvector = Xvector(poslist(:,8)>0,:); % the normal vector
    Yvector = Yvector(poslist(:,8)>0,:); % not used 
    Zvector = Zvector(poslist(:,8)>0,:); % the virus vector
    num = size(poslist,1);

        Avevector=sum(Zvector(:,:))./num;
        Vectorchange=acosd(dot(Avevector(1:2),Virusvector(1:2)));
        
        if  Vectorchange > 20
            fprintf('--- Virus vector not match ! skip !\n');
            continue
        else
            fprintf('--- Virus vector change %.1f\n',Vectorchange);
        end

        if ifplot >= 1
        nexttile; hold on;grid on; box on;axis equal;
        plot(poslist(:,1), poslist(:,2),'r.');
        quiver(poslist(:,1), poslist(:,2),Zvector(:,1)*100,Zvector(:,2)*100,'r');
        ss = norm(Virusvector);
        quiver(-Virusvector(1)*height/2/ss,-Virusvector(2)*height/2/ss,Virusvector(1)*height/ss,Virusvector(2)*height/ss,'b-');
        ss = norm(Avevector(1:2));
        quiver(-Avevector(1)*height/2/ss,-Avevector(2)*height/2/ss,Avevector(1)*height/ss,Avevector(2)*height/ss,'g-');
        end
        
       %
    Virusvector(3) = Avevector(3);
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
    
    % fitting for a, b and outplan rotate angle
    starta = radius;
    [am,bm,anglesearch] = meshgrid(0.86:0.02:1.00,0.6:0.05:1,-6:1.5:6);
    
    ratiolist = [am(:) bm(:) anglesearch(:)];
    rationum = size(ratiolist,1);
    fprintf('--- Total samplling number %d\n',rationum);
    score = zeros(rationum,3);
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
    t = -rounds*360:1:rounds*360;
    modelsize = size(t,2);

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
        ptCloud1 = pointCloud([poslist(:,1:2) poslist(:,3).*0]);
        ptCloud2 = pointCloud([rot3D(:,1:2) rot3D(:,3).*0]);
        tform = pcregistericp(ptCloud2,ptCloud1);
        %tform = pcregisterndt(ptCloud2,ptCloud1,0.5);
        
        %disp(tform.T)
        ptCloudTformed = pctransform(ptCloud2,tform);
        transferm = tform.T';

        % rotate the normals
        ptCloud3 = pointCloud(Back3Dv);
        bestshift = transferm(:)';
        besttform = rigid3d([bestshift(1:3);bestshift(5:7);bestshift(9:11)],[0 0 0]);
        tmpBack3Dv = pctransform(ptCloud3,besttform);
        Back3Dv = tmpBack3Dv.Location;

        % to save the transfer
        shiftlist(r,:) = transferm(:)'; % to save the shifts
        
        % to check the transfer
        if ifplot == 4 && mod(r,round(rationum/5))==1
            figure; hold on; axis equal;
            plot(ptCloudTformed.Location(:,1),ptCloudTformed.Location(:,2),'r.');
            plot(poslist(:,1),poslist(:,2),'go');
        end

        % find the nearest data points for model points
        diff_record1 = zeros(modelsize,2);

        
        % caculate the nearest model point for data points 
        diff_record2 = zeros(num,2);

        score(r,1:2) = [mean(diff_record1(:,1)) mean(diff_record2(:,1))];
        score(r,3) = mean(diff_record1(:,1))+mean(diff_record2(:,1));
        
        fprintf('diff %.1f - %.1f, %.1f ---\n',score(r,:));

    end

    %
    [~,minID]=min(score(:,3));
    a=starta*ratiolist(minID,1);
    b=a*ratiolist(minID,2);
    bestscore = score(minID,3);
    bestshift = shiftlist(minID,:);
    xangle = ratiolist(minID,3);

    Rxangle = [1,0,0; ...
        0,cosd(-xangle),-sind(-xangle);...
        0,sind(-xangle),cosd(-xangle)];

        if ifplot == 2
            figure; hold on; 
            plot3(ratiolist(:,1),ratiolist(:,2),score(:,3),'.');
            plot3(ratiolist(minID,1),ratiolist(minID,2),score(minID,3),'o')
        end
    
    % plot the best fit
        rounds = height/layer_inter2D/2*1;
        t = -rounds*360:1:rounds*360;
        modelsize = size(t,2);

        x3=a*cosd(t);
        y3=b*sind(t);
        z3=layer_inter*t/360;

        dx3=b*cosd(t);
        dy3=a*sind(t);
        %dx2=-a*sind(t);
        %dy2=b*cosd(t);

        rot3D=zeros(size(t,2),3);
        Back3Dv=zeros(size(t,2),3);

        for ii=1:size(t,2)
                rot3D(ii,:)=(Rz*Rx*Rxangle*[x3(ii),y3(ii),z3(ii)]')'; 
                
                Back3Dv(ii,:)=(Rz*Rx*Rxangle*[dx3(ii),dy3(ii),0]')';
                Back3Dv(ii,:)=Back3Dv(ii,:)./norm(Back3Dv(ii,:));
        end

        ptCloud2 = pointCloud(rot3D);
        besttform = rigid3d([bestshift(1:3);bestshift(5:7);bestshift(9:11)],bestshift(13:15));
        ptCloudTformed = pctransform(ptCloud2,besttform);
        rot3D = ptCloudTformed.Location;

        % rotate the normals
        ptCloud3 = pointCloud(Back3Dv);
        besttform = rigid3d([bestshift(1:3);bestshift(5:7);bestshift(9:11)],[0 0 0]);
        tmpBack3Dv = pctransform(ptCloud3,besttform);
        Back3Dv = tmpBack3Dv.Location;

        if ifplot >= 1
            figure(ff.Number);nexttile; hold on; axis equal;
            plot(ptCloudTformed.Location(:,1),ptCloudTformed.Location(:,2),'r.');
            plot(poslist(:,1),poslist(:,2),'go');
        end

        % check the normals
        if ifplot == 3
            figure; hold on; axis equal;

            for ii=round(modelsize*0.4):10:round(modelsize*0.6)
                plot3(rot3D(ii,1),rot3D(ii,2),rot3D(ii,3),'r.');
                quiver3(rot3D(ii,1),rot3D(ii,2),rot3D(ii,3),Back3Dv(ii,1)*100,Back3Dv(ii,2)*100,Back3Dv(ii,3)*100,'r');
            end
            
            view(45,45)
        end

            figure('Visible','off'); hold on; axis equal;
            plot(ptCloudTformed.Location(:,1),ptCloudTformed.Location(:,2),'r.');
            plot(poslist(:,1),poslist(:,2),'go');
            saveas(gcf,modelpng);

        % writeout a paramter file which can be used to generate the model
        fileIDfit = fopen(modelparam,'w');
        fprintf(fileIDfit,['%.1f %.1f %.1f %.1f %.1f %.4f %.4f %.4f %.4f '...
        '%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f '...
        '%.4f %.4f %.4f %.4f %.1f %.1f %.1f %.1f %.1f %.1f\n']...
        ,a,b,psi,theta,xangle,bestshift,centerrem,score(minID,1:3),rounds);
        fclose(fileIDfit);
        toc;
        
    fprintf('*** the best ratio %.2f - %.2f, Xangle %.1f, diff %.1f - %.1f, %.1f ***\n',ratiolist(minID,1),ratiolist(minID,2),ratiolist(minID,3),score(minID,1:3));

end % loop on models

% Step 3
% After genaration of 3D models,
% exclude the bad fitting manually,
% asign all sub-particles from 3D model,
% the Z value will also be asigned.
% Three output .star file are in RELION format,
% mod0 for cleaning particles only,
% mod1 for asign the Z high for cleaned particles,
% mod2 for asign Z high and Euler angles from model,
% mod3 for asign full particles from model.
% By Xing Zhang @20240106

%clear;clc

files=struct2cell(dir('./cryosparc_DW_*_pt*_model.png'));

% input parameters
apix = 0.87; % the pixel size of micrographs
layerinter = 80.5; % the average layer distance in pixel unit
posinter = 26.4; % the distance of two nearby positions in pixel unit
handness = -1; % whether to do a flip on handness 
dstep = 0.2; % the angle step (in degree) for caculation, 0.2 is good enough
discut = 9;  % within distance of X pixels
anglecut = 12; % within angle of X degree

% optional input 
dfoffset = 0; % the estimated defocus maybe not the same as virus, so need to addjust
ifplot = 1; % whether do plot the intermedia results for check

% run for all the particle files
for i=1:size(files,2)

    name=char(files(1,i));
    idend=strfind(name,'_model'); 
    namenoext=name(1:idend-1);
    fprintf('--- Model %s: ',namenoext);

    paramfile=sprintf('%s_model.log',namenoext);
    starfile_1img=sprintf('%s_srh.star',namenoext);
        
    % file for output

        asignpng = sprintf('%s_asign.png',namenoext);
        asignparam = sprintf('%s_asign.log',namenoext);

        particlestar0 = sprintf('%s_mod0.star',namenoext);
        particlestar1 = sprintf('%s_mod1.star',namenoext);
        particlestar2 = sprintf('%s_mod2.star',namenoext);
        particlestar3 = sprintf('%s_mod3.star',namenoext);

    [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12] = textread(starfile_1img,'%s %d %d %f %f %f %f %f %f %f %d %d'); 
    num=size(b2,1);
    poslist=[b2 b3 b3.*0];

        Xaxis = [1,0,0];
        Xvector = zeros(num,3);
        for j=1:num
            Rspa = BH_defineMatrix(-[b9(j) b8(j) b7(j)], 'SPIDER', 'inv');
            Rspa = reshape(Rspa,3,3);
            Xvector(j,1:3)=Rspa*Xaxis';
        end

    % read model paramters from file
    parameters = load(paramfile);
    a = parameters(1);
    b = parameters(2);
    psi = parameters(3);
    theta = parameters(4);
    xangle = parameters(5);
    bestshift = single(parameters(6:21));
    centerrem = parameters(22:23);
    %height = parameters(24);
    rounds = parameters(end)*1.2;

    % prepare global matrix
            Rz=[cosd(-psi),-sind(-psi),0;...
                sind(-psi),cosd(-psi),0;...
                0,0,1];
            Rx=[1,0,0; ...
                0,cosd(-theta),-sind(-theta);...
                0,sind(-theta),cosd(-theta)];
            Rxangle = [1,0,0; ...
                0,cosd(-xangle),-sind(-xangle);...
                0,sind(-xangle),cosd(-xangle)];
            besttform = rigid3d([bestshift(1:3);bestshift(5:7);bestshift(9:11)],bestshift(13:15));

    % prepare the model
    t = -rounds*360:dstep:rounds*360;
    modelsize = size(t,2);

    x3=a*cosd(t);
    y3=b*sind(t);
    z3=layerinter*t/360;
    
    dx3=b*cosd(t);
    dy3=a*sind(t);
    dx2=-a*sind(t);
    dy2=b*cosd(t);

    rot3D=zeros(size(t,2),3);
    Back3Dv=zeros(size(t,2),3);
    Back3Ds=zeros(size(t,2),3);
    
    % rotated 3D model
    for ii=1:modelsize
                rot3D(ii,:)=(Rz*Rx*Rxangle*[x3(ii),y3(ii),z3(ii)]')'; 
                
                Back3Dv(ii,:)=(Rz*Rx*Rxangle*[dx3(ii),dy3(ii),0]')';
                Back3Dv(ii,:)=Back3Dv(ii,:)./norm(Back3Dv(ii,:));
                
                Back3Ds(ii,:)=(Rz*Rx*Rxangle*[dx2(ii),dy2(ii),0]')';
                Back3Ds(ii,:)=Back3Ds(ii,:)./norm(Back3Ds(ii,:));
    end

        % rotate and shift the positions
        ptCloud2 = pointCloud(rot3D);
        ptCloudTformed = pctransform(ptCloud2,besttform);
        rot3D = ptCloudTformed.Location + ones(modelsize,1)*[centerrem 0];

        % rotate the normals
        ptCloud3 = pointCloud(Back3Dv);
        besttform = rigid3d([bestshift(1:3);bestshift(5:7);bestshift(9:11)],[0 0 0]);
        tmpBack3Dv = pctransform(ptCloud3,besttform);
        Back3Dv = tmpBack3Dv.Location;

        if ifplot == 1
            figure('Position',[680   558   1400   420]); 
            tiledlayout(1,3); nexttile; hold on; axis equal;
            plot(rot3D(:,1),rot3D(:,2),'r.');
            plot(b2(:),b3(:),'go');
        end

    % caculate Z and filtering particles

    model_ind = zeros(modelsize,3);
    model_ind(:,1) = (1:modelsize)';
    model_ind(:,4) = t';
    model_ind(1,2) = 0;

                for ii=2:modelsize
                    ss=norm(rot3D(ii,:)-rot3D(ii-1,:));
                    model_ind(ii,2)=model_ind(ii-1,2)+ss;
                end

    fileID1 = fopen(particlestar0,'w'); 
    fileID2 = fopen(particlestar1,'w'); 

    for ii=1:num
       
        % find the closest particles in model for each data point
            iffind = ( abs(poslist(ii,1)-rot3D(:,1))<=discut ) .* ( abs(poslist(ii,2)-rot3D(:,2))<=discut );
            iffindnum = sum(iffind);
            if iffindnum == 0
                continue
            else
                globaID = find(iffind);
            end

            modeltheta = acosd(dot( (Back3Dv(logical(iffind),1:3))' , (ones(iffindnum,1)*Xvector(ii,1:3))' ))';
            modeldist = sqrt( sum((rot3D(logical(iffind),1:2)-poslist(ii,1:2)).^2,2) );

            by2way = (modeldist<=discut) .* (modeltheta<=anglecut);

            if sum(by2way) == 0
                continue
            else
                localID = find(by2way);
            end
            [mindist,minID] = min(modeldist(logical(by2way)).*modeltheta(logical(by2way)));

            deltaZ = rot3D(globaID(localID(minID)),3);
            deltaDF = deltaZ*(-1)*apix+dfoffset;
            fprintf(fileID2,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii)+deltaDF,b5(ii)+deltaDF,b6(ii),b7(ii),b8(ii),b9(ii),b10(ii),b11(ii),b12(ii)); 
            fprintf(fileID1,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii),b5(ii),b6(ii),b7(ii),b8(ii),b9(ii),b10(ii),b11(ii),b12(ii)); 
            model_ind(globaID(localID(minID)),3) = 1;
            poslist(ii,3) = 1;
        
    end
    
        particles_keep = nnz(poslist(:,3)==1);
        fclose(fileID1);
        fclose(fileID2);

        if ifplot == 1
            nexttile; hold on; axis equal;
            plot(rot3D(model_ind(:,3)==1,1),rot3D(model_ind(:,3)==1,2),'r+'); 
            plot(poslist(poslist(:,3)==1,1),poslist(poslist(:,3)==1,2),'o','color',[0.25,0.25,1]);
        end

    %% check the filted particles

        % write the closest particles with asigned angles 
        fileID3 = fopen(particlestar2,'w');
        matched = find(model_ind(:,3)==1)';
        particles_matched = size(matched,2);

        for ii=1:particles_matched

            eularphi=atan2(dy3(matched(ii)),dx3(matched(ii)));
            angle=[-eularphi/pi*180,0,0];
            r = BH_defineMatrix(angle, 'Bah', 'forward');
            r = Rz*Rx*Rxangle*r; % for points, invert the order for map space?
            [angles] = rot_M2eZYZ(inv(reshape(angles,3,3)));

            fprintf(fileID3,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii)-deltaDF,b5(ii)-deltaDF,b6(ii),angles(1:3),b10(ii),b11(ii),b12(ii)); 

            if ii == 1
                continue
            end

                curvedis = model_ind(matched(ii),2)-model_ind(matched(ii-1),2);

                if curvedis > posinter*1.5
                    fills = round(curvedis/posinter)-1;
                    IDdis = round((matched(ii)-matched(ii-1))/(fills+1));
                    for k=1:fills
                        model_ind(matched(ii-1)+IDdis*k,3) = 2;
                    end
                end

        end
        fclose(fileID3);

        if ifplot == 1
            nexttile; hold on; axis equal;
            plot(rot3D(model_ind(:,3)>=1,1),rot3D(model_ind(:,3)>=1,2),'r+','color',[1,0.8,0]);
            plot(poslist(poslist(:,3)==1,1),poslist(poslist(:,3)==1,2),'bo'); 
        end
            figure('Visible','off'); hold on; axis equal;
            plot(rot3D(model_ind(:,3)>=1,1),rot3D(model_ind(:,3)>=1,2),'r+','color',[1,0.7,0]);
            plot(poslist(poslist(:,3)==1,1),poslist(poslist(:,3)==1,2),'o','color',[0.25,0.25,1]);
            saveas(gcf,asignpng);

        %% write all particles with asigned angles 
        fileID4 = fopen(particlestar3,'w');
        matched = find(model_ind(:,3)>=1)';
        particles_asign = size(matched,2);

        for ii=1:particles_asign



            % fprintf(fileID4,'%s %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.1f %d %d\n', char(b1(ii)),b2(ii),b3(ii),b4(ii)-deltaDF,b5(ii)-deltaDF,b6(ii),angles(1:3),b10(ii),b11(ii),b12(ii)); 

        end
        fclose(fileID4);
    
    % report particles number of each Virus 
    fprintf('%d read, %d kept, %d matched, %d asigned ---\n',num,particles_keep,particles_matched,particles_asign);

    fileIDfit = fopen(asignparam,'w');
    fprintf(fileIDfit,'%s %d %d %d %d\n',namenoext,num,particles_keep,particles_matched,particles_asign);
    fclose(fileIDfit);

end



function eu = rot_M2eZYZ(R)
% copy from https://github.com/KudryashevLab/dyn2rel) 
% citation:
% NATURE COMMUNICATIONS | (2020) 11:3709 | https://doi.org/10.1038/s41467-020-17466-0 | www.nature.com/naturecommunications

            eu = zeros(1,3);
            tol = 5e-5;

            if( R(3,3) < (1-tol))
                if( R(3,3) > (tol-1) )
                    % GENERAL CASE
                    eu(2) = acos ( R(3,3) )*180/pi;
                    eu(1) = atan2( R(2,3), R(1,3) )*180/pi;
                    eu(3) = atan2( R(3,2),-R(3,1) )*180/pi;

                else
                    % r22 <= -1
                    eu(2) = 180;
                    eu(1) = -atan2( R(2,1), R(2,2) )*180/pi;
                    eu(3) = 0;

                end
            else

                % r22 <= -1
                eu(2) = 0;
                eu(1) = atan2( R(2,1), R(2,2) )*180/pi;
                eu(3) = 0;

            end

end
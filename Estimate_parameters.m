% Step 1
% Crop out a virion from labeled local area, 
% estimate defocus from this local position, 
% estimate the angle of in plane rotatation, 
% and then ganerate a list of Eulor angles that cover the rotational range.

% read micrograph .mrc files and manual-pick .mod files 
% writeout two files .rem and .mtx,
% .rem keeps the particle area, defocus, center and vector
% .mtx keeps the list of rotation matrix for using in template matching.

% This script using BH_defineMatrix.m and BH_multi_gridCoordinates.m from emClarity
% So makesure to download source code of emClarity and copy the two functions into working folder

% By Xing Zhang @20240423

clear;clc

% input manual-pick as a format of IMOD .mod file
files=struct2cell(dir('./cryosparc_DW_*.mod'));

% input parameters
apix = 0.87; % the pixel size of micrographs
layer_inter = 81; % the average inter-layer distance
doctf = 1; % whether do a ctf estimation on local area
padadd = 600; % pad the micrographs to aviod edge problem

% optional input, whether to plot a figure for check 
ifplot = 0; 

% keep this to 1
interbin = 1; 

padadd = [padadd padadd+(5760-4092)/2 0];

samplingRate = 1; % for mod and star file
radius = round(800/apix/samplingRate/interbin);
boxsize = radius*2+1;

% prepare the masks
mask1 = creatmask2D(boxsize,20);
mask2 = creatmask2D(boxsize,50);
maskforvector = mask2 - mask1;

modellog = sprintf('all_Virus_models.rem');
fileID = fopen(modellog,'w');

%% run for all the particle files

for i=1:size(files,2)
    
    name=char(files(1,i));
    fprintf('\n*** Processing on No.%d %s ***\n',i,name);
    idend=strfind(name,'.mod');
    name=name(1:idend-1);
    
    imgname = sprintf('%s-img-pad.mrc',name);
    pickmod = sprintf('%s.mod',name);

    % read the mod
    system(sprintf('cp %s tmp.mod',pickmod));
    system(sprintf('model2point tmp.mod Virusforread.txt'));
    Viruspick = load('Virusforread.txt')./interbin;

    Viruspick(:,1:3) = (Viruspick(:,1:3)+padadd)./interbin;

    if mod(size(Viruspick,1),4)~=0 
        fprintf('picking maybe wrong, double check!');
        fprintf(fileID,'%d %s picking maybe wrong \n',i,name);
        continue
    end
    
    modnum = size(Viruspick,1)/4;

    % read the image
    m = MRCImage(sprintf('%s', imgname),0);
    imageorg = single(getVolume(m));
    [size1,size2]=size(imageorg);
    image = single(imresize(imageorg,[round(size1/interbin) round(size2/interbin)],'lanczos3'))*(-1);

for imodel=1:modnum
    
    % caculate center position, direction, height and width of each target 
    Viruscenter = (Viruspick((imodel-1)*4+1,:)+Viruspick((imodel-1)*4+2,:))./2;
    Viruscenter = round(Viruscenter);
    Virusvector = (Viruspick((imodel-1)*4+2,:)-Viruspick((imodel-1)*4+1,:));
    Virusvectormanual = Virusvector;
    Virusheight = norm(Virusvector)*interbin;
    
    try
        localimg = image(Viruscenter(1)-radius:Viruscenter(1)+radius,Viruscenter(2)-radius:Viruscenter(2)+radius);
    catch
        fprintf(fileID,'%d %s cropping out of images size\n',i,name);
        continue
    end
    
    localfft = abs(fftshift(fftn(localimg))) ;
    
    [curve,imgpsi] = find_Virus_dirction(localfft,maskforvector);

    %SAVE_IMG( MRCImage(gather(localfft( ((end+1)/2-100):((end+1)/2+100),((end+1)/2-100):((end+1)/2+100) ))),sprintf('%s_pt%d_fft.mrc',name,imodel) ); 
    
    Viruspsi = imgpsi+90;
    if dot([cosd(Viruspsi) sind(Viruspsi)],Virusvector(1:2) ) < 0
        Viruspsi = mod(Viruspsi+180,360);
    end
    
    Virusvector = [cosd(Viruspsi),sind(Viruspsi),0];    
    
    widthmanual = norm(dot((Viruspick((imodel-1)*4+3,:)-Viruspick((imodel-1)*4+4,:)),[cosd(imgpsi) sind(imgpsi) 0] ))*interbin;
    Viruside = [cosd(imgpsi),sind(imgpsi),0].*sign(cosd(imgpsi))*sign(Virusvector(1));
    side1 = dot(Viruspick((imodel-1)*4+3,:)-Viruscenter,Viruside)*Viruside;
    side2 = dot(Viruspick((imodel-1)*4+4,:)-Viruscenter,Viruside)*Viruside;
    newcentermanul = (Viruscenter+(side1+side2)/2)*interbin;
        
    angledis=acosd(dot(Virusvector,Virusvectormanual)/norm(Virusvectormanual)/norm(Virusvector));
    if angledis >=15 % in degree
        fprintf(' --- vector too different from manualpick! %.1f \n',angledis)
        fprintf(fileID,'%d %s pt %d angle change %.1f \n',i,name,imodel,angledis);
        continue
    end
    
    Varea(1,1:3)=newcentermanul+Virusvector*Virusheight/2+cross(Virusvector,[0 0 1])* widthmanual/2; %nrom(Virusvector)
    Varea(2,1:3)=newcentermanul+Virusvector*Virusheight/2+cross(Virusvector,[0 0 -1])* widthmanual/2; 
    Varea(3,1:3)=newcentermanul-Virusvector*Virusheight/2+cross(Virusvector,[0 0 1])* widthmanual/2; 
    Varea(4,1:3)=newcentermanul-Virusvector*Virusheight/2+cross(Virusvector,[0 0 -1])* widthmanual/2; 
    Varea = round(Varea);
    
    % may be should crop out virus after whiting
    VirusX = ceil(Virusheight+2); % y direction
    VirusY = ceil(widthmanual+200); % x direction
    
    centerX = VirusX/2;
    centerY = VirusY/2;
    
    [meshX,meshY]= meshgrid(1:VirusX,1:VirusY);
    newmeshX = cosd(Viruspsi)*(meshX-centerX)-sind(Viruspsi)*(meshY-centerY)+newcentermanul(1);
    newmeshY = sind(Viruspsi)*(meshX-centerX)+cosd(Viruspsi)*(meshY-centerY)+newcentermanul(2);
    Virusrot = interpn(1:(size1),1:(size2), imageorg, newmeshX,newmeshY,'cubic');
    SAVE_IMG(MRCImage(gather(Virusrot)),sprintf('%s_pt%d_rot.mrc',name,imodel) ); 
    
    if doctf == 1
        fileID4 = fopen('run_ctffind4_tmp.sh','w');
        fprintf(fileID4,'#!/bin/csh -f\n');
        fprintf(fileID4,'ctffind > ctffind_tmp.txt << EOSc\n');
        fprintf(fileID4,'$1\ndiag_temp.mrc\n0.87\n300\n0.01\n0.07\n256\n35\n5\n5000\n');
        fprintf(fileID4,'20000\n100\nno\nno\nyes\n1000\nno\nno\nEOSc\n\n');
        fclose(fileID4);
        
    system(sprintf('chmod +x run_ctffind4_tmp.sh'));
    system(sprintf('./run_ctffind4_tmp.sh %s_pt%d_rot.mrc',name,imodel) );
    system(sprintf('awk ''!/^#/{ print $0 }'' diag_temp.txt > diag_temp_simple.txt'));
    ctffind = load(sprintf('diag_temp_simple.txt') );
    defFit1 = ctffind(2);
    defFit2 = ctffind(3);
    angleFit = ctffind(4)+Viruspsi; 
    ctfres = ctffind(7);
    system(sprintf('rm diag_temp*'));
    end
    
    % caculate layer distance from Fourier space
    [freqlayer,~] = find_Virus_layer(localfft,maskforvector,imgpsi,ifplot);
     layerdist_f = interbin./freqlayer;
    
    Virusrot([1:130 end-129:end],:) = 0;
    Virusrot(round(end/2-widthmanual/2+230):round(end/2+widthmanual/2-230),:) = 0;   

    % caculate layer distance from real space
    [layerdist_r,layermap] = find_Virus_layer_CC(Virusrot,2,ifplot);

    %SAVE_IMG(MRCImage(layermap(((end+1)/2-100):((end+1)/2+100),((end+1)/2-100):((end+1)/2+100))),sprintf('%s_pt%d_layer.mrc',name,imodel) ); 
    
    % writeout matrix file
    matrixname = sprintf('%s_pt%d_mtx.txt',name,imodel);
    samplingfile =  'relion_sampling.star' ;
    usepsi = (-Viruspsi)-180;
    range_psi = 5; % 7.2 in run4, 0 in run5, 5 in run5.6
    step_psi = 1.8; % 1.8 in run4
    range_theta = 12; % 15 in run4
    usetheta = acosd(min(layerdist_r,layer_inter)/layer_inter);
    
    % writeout rem file
    remname = sprintf('%s_pt%d_rem.txt',name,imodel);
    fileID0 = fopen(remname,'w');
    fprintf(fileID0,'%d %d %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.4f %.4f %.4f %.2f %.2f %.1f %.1f\n',min(Varea(:,1)),max(Varea(:,1)),min(Varea(:,2)),max(Varea(:,2)),defFit1,defFit2,ctfres,newcentermanul,Virusvector,widthmanual,Virusheight,layerdist_r,usetheta);
    fclose(fileID0);
    
    fprintf(fileID,'%d %s pt %d %.1f %.1f %.4f %.4f %.1f %.1f %.1f %.1f\n',i,name,imodel,newcentermanul(1:2),Virusvector(1:2),widthmanual,Virusheight,layerdist_r,usetheta);
    fprintf('pt %d: center %.1f %.1f, vector %.4f %.4f, width %.1f, height %.1f, layer_dist %.1f / %.1f, theta_use %.1f\n',imodel,newcentermanul(1:2),Virusvector(1:2),widthmanual,Virusheight,layerdist_r,layerdist_f,usetheta);

    fileID1 = fopen(matrixname,'w');
    
    [angleRot,angleTilt] = textread(samplingfile,'%f %f ','headerlines',25);

    anglenum = size(angleRot,1);
    psilist = usepsi-range_psi:step_psi:usepsi+range_psi;
    psinum = size(psilist,2);
    anglecount = 1;
    
    for num = 1:anglenum 
        for psi = psilist 
        rot = [ angleRot(num),angleTilt(num),psi ];
            %if angleRot(num)>=0
            %if abs(angleTilt(num)-90)<=max(range_theta,usetheta)
            if abs(abs(angleTilt(num)-90) - usetheta) < range_theta
                R=BH_defineMatrix(rot, 'SPIDER', 'inv');
                for Ri=1:3
                fprintf(fileID1,'  %d %.4f %.4f %.4f \n',anglecount,R(Ri,:));
                end
                anglecount = anglecount+1;
            end
        end
    end
    fclose(fileID1);
    
end

end

fclose(fileID);

% find membrane dirction
function [curve,imgpsi] = find_Virus_dirction(imgorg,mask)

imgorg = imgorg.*mask;
imgmean=mean(mean(imgorg(mask>0.5)));
imgorg(mask>0.5) = imgorg(mask>0.5)-imgmean;

imgpsi = 0:0.1:179.5;
[R,~] = radon(imgorg,imgpsi);

[~,maxind] = max( R( (size(R,1)+1)/2 ,:) );

curve = R(:,maxind);
imgpsi = imgpsi(maxind);

end

function [freqlayer,fpixel] = find_Virus_layer(imgfft,mask,imgpsi,ifplot)

    imgfft = imgfft.*mask;
    [radialForCTF,~,~,~,~,~] = BH_multi_gridCoordinates(size(imgfft),'Cylindrical','cpu',{'none'},1,1,0);

% increase the accuracy !
center = (size(imgfft,1)+1)/2;
imgfftexpand = imresize(imgfft(center-35:center+35,center-35:center+35),[710 710]);
radialForCTFexp = imresize(radialForCTF(center-35:center+35,center-35:center+35),[710 710]);

[R,~] = radon(imgfftexpand,imgpsi);
Rsize = (size(R,1)+1)/2;
Rplot = R;
R(round(end*2/5):round(end*3/5))=0;

[peakss,peaklocal] = findpeaks(R);
peaklist = [peakss,peaklocal];

peaklist = sortrows(peaklist,1,'descend');

peak1 = peaklist(1,2);
peak2 = peaklist(2,2);

fpixel = peak1-Rsize;

peak1c = round(peak1-Rsize+356.0);

peak2c = round(peak2-Rsize+356.0);

freqlayer1 = radialForCTFexp(peak1c,356);
freqlayer2 = radialForCTFexp(peak2c,356);
freqlayer = mean([freqlayer1,freqlayer2]);

if ifplot == 1
figure; 
plot(Rplot);
hold on;
plot([peak1 peak1],[min(R) max(R)],'--g');
plot([peak2 peak2],[min(R) max(R)],'--r');
xlim([1,size(R,1)]);
end

end

function [layerdist,Xcorrimg] = find_Virus_layer_CC(imgrot,bin,ifplot)

size1=size(imgrot,1)/bin;
size2=size(imgrot,2)/bin;
imgrot = single(imresize(imgrot,[round(size1/bin) round(size2/bin)],'lanczos3'));

    Xcenter = size(imgrot(:,:,1),1);
    Ycenter = size(imgrot(:,:,1),2);
    
    % real space   
        Xcorrimg = xcorr2(imgrot,imgrot); 
        R = Xcorrimg([(Xcenter-4):(Xcenter-1) (Xcenter+1):(Xcenter+4)],:);
        R = sum(R,1);   
        R(round(Ycenter-60/bin):round(Ycenter+60/bin)) = 0;

        Xcorrimg([1:round(Xcenter-20/bin) round(Xcenter+20/bin):end],:) = 0; 
        Xcorrimg(:, round(Ycenter-70/bin):round(Ycenter+70/bin)) = 0; 
        Xcorrimg(:, [1:round(Ycenter-90/bin) round(Ycenter+90/bin):end]) = 0; 
        
    [~,index] = max(Xcorrimg(:));
    [x,y] = ind2sub(size(Xcorrimg),index);
    
    % subpixel acuracy
    peakCoord = [x,y];
    peakCOM = [1,1].*2; % 
    peakLOW = peakCoord - peakCOM;
    peakTOP = peakCoord + peakCOM;
    
    [cmX, cmY] = ndgrid(-1*peakCOM(1):peakCOM(1), -1*peakCOM(2):peakCOM(2));
    boX = Xcorrimg(peakLOW(1):peakTOP(1), peakLOW(2):peakTOP(2));
    boX = boX - min(min(boX));
    cMass = [ sum(sum(boX.*cmX)) sum(sum(boX.*cmY)) ] ./ sum(boX(:));
            
    layerdist = abs((Ycenter-y-cMass(2))*bin);
        
if ifplot == 1
figure; 
plot(R);
hold on;
plot([y-cMass(2) y-cMass(2)],[min(R) max(R)],'--g');
plot([Ycenter*2-y+cMass(2) Ycenter*2-y+cMass(2)],[min(R) max(R)],'--r');
xlim([1,size(R,2)]);
end

end

function mask2d = creatmask2D(imgsize,rad)
    
    Xcenter = (imgsize+1)/2;
    mask2d = zeros(imgsize,imgsize);
    for x=1:imgsize
        for y=1:imgsize
            dis = sqrt((x-Xcenter)^2+(y-Xcenter)^2);
            if dis <= rad
                mask2d(x,y)=1;
            end
        end
    end
end

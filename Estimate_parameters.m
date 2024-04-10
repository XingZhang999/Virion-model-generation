% Virus sample, perpare for search 2D
% crop out a virus labeled 
% estimate defocus from this local position 
% estimate the angle of rotatation 
% ganerate a list of Eulor angle that cover the range
% writeout two files .rem and .mtx
% .rem keeps the particle area, defocus, center and vector
% .matx keeps the posible matrix of eulor angles for this particle
% Xing Zhang @20210106

clear;clc
% addpath(genpath('/ssd/Software/emClarity/emClarity_Aug26_noPTmask_Bfactor_smoothXY_new'))

% need to estimate defocus from each virus:
% 1, use virus related defocus change compare to Micrograph
% 2, search defocus handness by template matching ? not accurate !

padadd=600
apix = 0.87
doctf = 1
ifplot = 0 
ifcheck = 0 
interbin= 5 

padadd=[padadd padadd+(5760-4092)/2 0]
padaddbin=padadd./interbin;

samplingRate = 1 % for mod and star file
radius = round(800/apix/samplingRate/interbin);
boxsize = radius*2+1;

%system('rm cryosparc_DW*_pt?.*');
files=struct2cell(dir('./cryosparc_DW_*.mod'));

modellog = sprintf('all_Virus_models.rem');
fileID = fopen(modellog,'w');

for i=1:size(files,2)
    
    name=char(files(1,i));
    fprintf('\n*** Processing on No.%d %s ***\n',i,name);
    idend=strfind(name,'.mod');
    name=name(1:idend-1);
    
    imgname = sprintf('%s-ice-pad.mrc',name);
    pickmod = sprintf('%s.mod',name);

system(sprintf('cp %s tmp.mod',pickmod));
system(sprintf('model2point tmp.mod Virusforread.txt'));
Viruspick = load('Virusforread.txt')./interbin;

Viruspick(:,1:3) = Viruspick(:,1:3)+padaddbin;

if ifcheck==1
    if mod(size(Viruspick,1),4)~=0 
        fprintf('picking maybe wrong, double check!');
        return
    else
        continue
    end
end
    
modnum = size(Viruspick,1)/4;

% read image
m = MRCImage(sprintf('%s', imgname),0);
imageorg = single(getVolume(m));
[size1,size2]=size(imageorg);
image = single(imresize(imageorg,[round(size1/interbin) round(size2/interbin)],'lanczos3'))*(-1);

mask1 = creatmask2D(boxsize,10);
mask2 = creatmask2D(boxsize,50);
maskforvector = mask2 - mask1;
SAVE_IMG(MRCImage(gather(maskforvector)),'check-test-maskforvector.mrc' ); 

mask1 = creatmask2D(boxsize,17);
mask2 = creatmask2D(boxsize,35);
maskforwidth = mask2 - mask1;
SAVE_IMG(MRCImage(gather(maskforwidth)),'check-test-maskforwidth.mrc' ); 

for imodel=1:modnum
    
    % caculate center posotion, direction and width of each labeled virus 
    Viruscenter = (Viruspick((imodel-1)*4+1,:)+Viruspick((imodel-1)*4+2,:))./2;
    Viruscenter = round(Viruscenter);
    Virusvector = (Viruspick((imodel-1)*4+2,:)-Viruspick((imodel-1)*4+1,:));
    Virusvectormanual = Virusvector;
    Virushight = norm(Virusvector)*interbin;
    
    if (Viruscenter(1)-radius>=1)*(Viruscenter(1)+radius<=size(image,1))*(Viruscenter(2)-radius>=1)*(Viruscenter(2)+radius<=size(image,2)) == 0
        fprintf(' --- model %d too close to edge ! skip ! \n',imodel);
        continue
    end
    
    localimg = image(Viruscenter(1)-radius:Viruscenter(1)+radius,Viruscenter(2)-radius:Viruscenter(2)+radius);
    SAVE_IMG(MRCImage(gather(localimg)),'check-test-localimg.mrc' ); 
    
    localfft = abs(fftshift(fftn(localimg))).*maskforvector;
    
    [curve,imgpsi] = find_Virus_dirction(localfft,maskforvector);

    SAVE_IMG(MRCImage(gather(localfft)),'check-test-localfft.mrc'); 
    SAVE_IMG(MRCImage(gather(localfft)),sprintf('%s_pt%d_fft.mrc',name,imodel) ); 
    
    [freqlayer,fpixel] = find_Virus_layer(localfft.*maskforwidth,imgpsi,ifplot);
    layerdist = interbin./freqlayer;
    
    masklayer = creatmask2D_theta(boxsize,imgpsi,round(fpixel));
    
    localimg = real(ifftn(ifftshift(fftshift(fftn(localimg)).*masklayer))); %maskforwidth
    SAVE_IMG(MRCImage(gather(localimg)),'check-test-localimg-processed.mrc' ); 
    
    anglevector = imgpsi+90;
    if dot([cosd(anglevector) sind(anglevector)],Virusvector(1:2) ) < 0
        anglevector = mod(anglevector+180,360);
    end
    
    %if sign(cosd(anglevector)) ~=0
    %Virusvector = [cosd(anglevector),sind(anglevector),0].*sign(cosd(anglevector))*Virusvector(1)/abs(Virusvector(1));
    %elseif sign(Virusvector(1)) ~=0
    %Virusvector = [cosd(anglevector),sind(anglevector),0].*sign(Virusvector(1));
    %else
    Virusvector = [cosd(anglevector),sind(anglevector),0];    
    %end
    
    [width,newcenter] = find_Virus_width(localimg,anglevector,Viruscenter,ifplot); 
    newcenter=newcenter*interbin;
    width = width*interbin;
    
    widthmanual = norm(cross((Viruspick((imodel-1)*4+3,:)-Viruspick((imodel-1)*4+4,:)),Virusvector ))*interbin;
    Viruside = [cosd(imgpsi),sind(imgpsi),0].*cosd(imgpsi)/abs(cosd(imgpsi))*Virusvector(1)/abs(Virusvector(1));
    side1 = dot(Viruspick((imodel-1)*4+3,:)-Viruscenter,Viruside)*Viruside;
    side2 = dot(Viruspick((imodel-1)*4+4,:)-Viruscenter,Viruside)*Viruside;
    newcentermanul = (Viruscenter+(side1+side2)/2)*interbin;
    %newcentermanul = newcentermanul(1:2);
    
    fprintf('center %f %f | %f %f --- vector %f %f --- width %f | %f --- layerdist %f \n',newcenter(1:2),newcentermanul(1:2),Virusvector(1:2),width,widthmanual,layerdist);

    width = widthmanual;
    shift2center = 106;
    
    angledis=acosd(dot(Virusvector,Virusvectormanual)/norm(Virusvectormanual)/norm(Virusvector));
    if angledis >=15 %degree
        fprintf(' --- vector too different from manualpick! skip !\n')
        continue
    end
    
    if layerdist>=(79/apix) || layerdist<=(65/apix)
        fprintf(' --- layerdist too different ! skip !\n')
        continue
    end
    
    if width>=(1050/apix) || width<=(620/apix)
        fprintf(' --- width too different ! skip !\n')
        continue
    end

    % now record as a good virus
    fprintf(fileID,'%d %d %s %d %s %f %f %f %f %f %f %f %f\n',1,i,name,imodel,sprintf('%s_pt%d',name),newcentermanul,Virusvector,widthmanual,layerdist);
    
    % ganerate 3D model for each Virus
    %if Virusvector(1)>=0
    %Viruspsi=acosd(Virusvector(1)/sqrt(Virusvector(1)^2+Virusvector(2)^2));
    %end
    %if Virusvector(1)<0
    %Viruspsi=acosd(-Virusvector(1)/sqrt(Virusvector(1)^2+Virusvector(2)^2));
    %end
    
    Viruspsi=atan2d(Virusvector(2),Virusvector(1))
    
    Varea(1,1:3)=newcentermanul+Virusvector*Virushight/2+cross(Virusvector,[0 0 1])* widthmanual/2; %nrom(Virusvector)
    Varea(2,1:3)=newcentermanul+Virusvector*Virushight/2+cross(Virusvector,[0 0 -1])* widthmanual/2; 
    Varea(3,1:3)=newcentermanul-Virusvector*Virushight/2+cross(Virusvector,[0 0 1])* widthmanual/2; 
    Varea(4,1:3)=newcentermanul-Virusvector*Virushight/2+cross(Virusvector,[0 0 -1])* widthmanual/2; 
    Varea = round(Varea);
    
    % may be should crop out virus after whiting
    VirusX = ceil(Virushight+2); % y direction
    VirusY = ceil(widthmanual+2); % x direction
    
    centerX = VirusX/2;
    centerY = VirusY/2;
    
    [meshX,meshY]= meshgrid(1:VirusX,1:VirusY);
    newmeshX = cosd(Viruspsi)*(meshX-centerX)-sind(Viruspsi)*(meshY-centerY)+newcentermanul(1);
    newmeshY = sind(Viruspsi)*(meshX-centerX)+cosd(Viruspsi)*(meshY-centerY)+newcentermanul(2);
    Virusrot = interpn(1:(size1),1:(size2), imageorg, newmeshX,newmeshY,'cubic');
    SAVE_IMG(MRCImage(gather(Virusrot)),sprintf('%s_pt%d_rot.mrc',name,imodel) ); 
    
    if doctf == 1
    system(sprintf('./PRE-ctffind4.sh %s_pt%d_rot.mrc',name,imodel) );
    system(sprintf('awk ''!/^#/{ print $0 }'' diag_temp.txt > diag_simple.txt'));
    ctffind = load(sprintf('diag_simple.txt') );
    defFit1 = ctffind(2);
    defFit2 = ctffind(3);
    angleFit = ctffind(4)+Viruspsi % will be replaced with image ctf angle
    ctfres = ctffind(7);
    %Vdef = [defFit1 defFit2 angleFit];
    end
    
    % writeout rem file
    remname = sprintf('%s_pt%d_rem.txt',name,imodel);
    fileID0 = fopen(remname,'w');
    fprintf(fileID0,'%d %d %d %d %.1f %.1f %.1f %.1f %.1f %.1f %.2f %.2f %.2f\n',min(Varea(:,1)),max(Varea(:,1)),min(Varea(:,2)),max(Varea(:,2)),defFit1,defFit2,ctfres,newcentermanul,Virusvector);
    fclose(fileID0);
    
    % writeout matrix file
    matrixname = sprintf('%s_pt%d_mtx.txt',name,imodel);
    samplingfile =  'run_it017_sampling.star' %'run_it000_sampling.star''run_it007_sampling.star'
    usepsi = (-Viruspsi)-180
    range_psi = 7.2
    step_psi = 1.8 % 6 
    range_theta = 15
    
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
            if abs(angleTilt(num)-90)<=range_theta
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

imgmean=mean(mean(imgorg(mask>0.5)));
imgorg(mask>0.5) = imgorg(mask>0.5)-imgmean;

imgpsi = 0:0.1:179.5;
[R,~] = radon(imgorg,imgpsi);
%[weight,~] = radon(mask,0);
SAVE_IMG(MRCImage(gather(R)),'check-test-localradon.mrc' ); 
%[~,maxind] = max(max(R,[],1));
%size(R)
[~,maxind] = max( R( (size(R,1)+1)/2 ,:) );
%maxind
curve = R(:,maxind);
imgpsi = imgpsi(maxind);

end

function [freqlayer,peak] = find_Virus_layer(imgfft,imgpsi,ifplot)

    [radialForCTF,~,~,~,~,~] = ...
     BH_multi_gridCoordinates(size(imgfft),'Cylindrical','cpu',{'none'},1,1,0);
    fposition = (size(imgfft,1)-1)/2+1;
    freqx = radialForCTF(fposition:size(imgfft),fposition);

[R,~] = radon(imgfft,imgpsi);

Rsize = (size(R,1)+1)/2;
SAVE_IMG(MRCImage(gather(R)),'check-test-fftlayer.mrc' ); 

% increase the accuracy !
center = (size(imgfft,1)+1)/2;
imgfftexpand = imresize(imgfft(center-35:center+35,center-35:center+35),[701 701]);
SAVE_IMG(MRCImage(gather(imgfftexpand)),'check-test-localfft-expand.mrc' ); 
[R,~] = radon(imgfftexpand,imgpsi);
Rsize = (size(R,1)+1)/2;
SAVE_IMG(MRCImage(gather(R)),'check-test-fftlayer-expand.mrc' ); 

[~,peak1] = max(R(1:Rsize-10));
peak1c = Rsize-peak1;
[~,peak2] = max(R(Rsize+10:end));
peak2c = peak2+9;
peak = (peak1c+peak2c)/2/10;
%freqlayer = freqx(round(peak)-1)
%size(freqx)
freqlayer = interpn(1:length(freqx),freqx(:)',peak);

if ifplot == 1
figure(1)
plot(R);
hold on;
plot([peak1 peak1],[min(R) max(R)],'--g');
plot([peak2+Rsize+9 peak2+Rsize+9],[min(R) max(R)],'--r');
end

end

function [XSH] = line_max(Z)
    Z
	C2 = (0 + 5.*Z(1) - 8.*Z(2) + 3.*Z(3) +0)/(-6.);
	C3 = (0 + Z(1) -2.*Z(2) + Z(3) + 0)/6.;
	C4 = (0 -8.*Z(1) -8.*Z(2) - 8.*Z(3) +0)/(-6.);
	C5 = (0)/4.;
	C6 = (0- 2.*Z(1) - 2.*Z(2) - 2.*Z(3) +0)/6.;
	
	XSH   = 0.0;
    DENOM = 4. * C3 * C6 - C5 * C5;
	if DENOM == 0 
        return
	end
	XSH = (C2*C5 - 2.*C4*C3) / DENOM - 2.

	if (XSH < -1.) ; XSH = -1 ; end
	if (XSH > +1.) ; XSH = 1 ; end

end


function [width,center] = find_Virus_width(imgorg,theta,center,ifplot)

imgmean=mean(mean(imgorg));
imgorg= imgorg-imgmean;
imgorg=imgorg./std(std(imgorg));
imgorg(imgorg<0)=imgorg(imgorg<0)*(-1);

%theta = 0:1:179;
[R,~] = radon(imgorg,theta);
%[weight,~] = radon(mask,0);
SAVE_IMG(MRCImage(gather(R)),'check-test-radontheta.mrc' ); 
%[~,maxind] = max(max(R,[],1));

th = mean(R)+std(R)*1.5;

if ifplot == 1
figure(2)
plot(R);
hold on;
plot(1:size(R,1),ones(size(R,1))*(th),'--r');
end

startp = find(R>th,1,'first');
endp = find(R>th,1,'last');
Rsize = (size(R,1)+1)/2;
centershift  = (startp+endp)/2 - Rsize;
width = endp-startp;
center = center + ([centershift*cosd(theta-90) centershift*sind(theta-90) 0]);
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
    %SAVE_IMG(MRCImage(mask2d),sprintf('check-mask2D-sharp.mrc'));
end

function mask2d = creatmask2D_theta(imgsize,theta,fpixel)
    
    Xcenter = (imgsize+1)/2;
    mask2d = zeros(imgsize,imgsize);
    for x=1:imgsize
        for y=Xcenter-fpixel-2:Xcenter-fpixel+2
                mask2d(x,y)=1;
        end
        for y=Xcenter+fpixel-2:Xcenter+fpixel+2
                mask2d(x,y)=1;
        end
    end
    mask2d = imrotate(mask2d,theta,'crop');
    %mask2d = 1-mask2d;
    SAVE_IMG(MRCImage(mask2d),sprintf('check-test-masktheta.mrc'));
end


function eu = rot_M2eZYZ(R)
% added by Xing Zhang @20210928
% copy from https://github.com/KudryashevLab/dyn2rel) 
% which also used in
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

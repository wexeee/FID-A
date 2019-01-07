%io_loadspec_twix.m
%Jamie Near, McGill University 2014.
%Edits from Franck Lamberton, 2017.
%Major edits from Will Clarke in 2018 to add in imaLoopAssignment argument
%
% USAGE:
% out=io_loadspec_twix(filename,imaLoopAssignment);
% 
% DESCRIPTION:
% Reads in siemens twix raw data (.dat file) using the mapVBVD.m and 
% twix_map_obj.m functions from Philipp Ehses (philipp.ehses@tuebingen.mpg.de).
% 
% op_loadspec_twix outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% filename   = filename of Siemens twix data to load.
% imaLoopAssignment = 
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out=io_loadspec_twix_wtc(filename,imaLoopAssignment)
%read in the data using the new mapVBVD.  This code has been adapted to 
%handle both single RAID files and multi-RAID files.  The vast majority of
%Siemens twix data comes as a single RAID file, but I've encoundered a few 
%multi-RAID files, particularly when using VD13D.  The way to distinguish
%them here is that a for a single RAID file, mapVBVD will output a struct, 
%whereas for a multi-RAID file, mapVBVD will output a cell array of structs.
%This code assumes that the data of interest is in the last element of the 
%cell array (possibly a bad assumption under some circumstances):

if nargin ==1
    imaLoopAssignment = {};
elseif iscell(imaLoopAssignment)&&isempty(imaLoopAssignment)
    
else
    if ~iscell(imaLoopAssignment)||~all(size(imaLoopAssignment)==[2 5])
        error('imaLoopAssignment must be a cell array of strings in the format {''t'',''coils'',''averages'',''subSpecs'',''extras'';'' '','' '','' '','' '','' ''}');
    end
end

twix_obj=mapVBVD(filename);
if isstruct(twix_obj)
    disp('single RAID file detected.');
    RaidLength=1;
elseif iscell(twix_obj)
    disp('multi RAID file detected.');
    RaidLength=length(twix_obj);
    %assume that the data of interest is in the last element of the cell.
    twix_obj=twix_obj{RaidLength};
end
dOut.data=twix_obj.image();
version=twix_obj.image.softwareVersion;
sqzSize=twix_obj.image.sqzSize; 
sqzDims=twix_obj.image.sqzDims;

data=dOut.data;

%Squeeze the data to remove singleton dims
fids=squeeze(data);

% Scanner version
isVB = strcmp(version,'vb');
isVD = strcmp(version,'vd');
isVE = strcmp(version,'ve');

% SEQUENCE IDENTIFICATION - for some special case handling and if the 
% second input is empty then check if the sequence is regonised from
% the preexisting list and use the information stored there. Otherwise use
% the default

%find out what sequence, the data were acquired with.  If this is a
%multi-raid file, then the header may contain multiple instances of
%'tSequenceFileName' for different scans (including a pre-scan).
%Therefore, if multi-raid file, we will need to do a bit of extra digging
%to find the correct sequence name.
sequence=twix_obj.hdr.Config.SequenceFileName;

%Try to find out what sequnece this is:
if ~isempty(strfind(sequence,'rm_special')) ||...  %Is this Ralf Mekle's SPECIAL sequence?
        ~isempty(strfind(sequence,'vq_special')) ||... %or the CIBM SPECIAL sequence?
        ~isempty(strfind(sequence,'jn_svs_special'));  %or Jamie Near's SPECIAL sequence?
    % isSpecial
    if isempty(imaLoopAssignment)
        if isVD || isVE
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Ave','Set',''};
        else
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Set','Ida',''};
        end
    end

    %If this is the SPECIAL sequence, it probably contains both inversion-on
    %and inversion-off subspectra on a single dimension, unless it is the VB
    %version of Jamie Near's SPECIAL sequence, in which case the subspecs are
    %already stored on separate dimensions.
    %Both Ralf Mekle's SPECIAL and the VD-VE version of Jamie Near's SPECIAL sequence
    %do not store the subspectra along a separate dimension of the data array,
    %so we will separate them artifically:
    isjnseq = ~isempty(strfind(sequence,'jn_'));
    if ~(isVB && isjnseq) %Catches any SPECIAL sequence except Jamie Near's VB version.
        squeezedData=squeeze(dOut.data);
        if twix_obj.image.NCol>1 && twix_obj.image.NCha>1
            data(:,:,:,1)=squeezedData(:,:,1:2:(end-1));
            data(:,:,:,2)=squeezedData(:,:,2:2:end);
            sqzSize=[sqzSize(1) sqzSize(2) sqzSize(3)/2 2];
        elseif twix_obj.NCol>1 && twixObj.image.NCha==1
            data(:,:,1)=squeezedData(:,1:2:(end-1));
            data(:,:,2)=squeezedData(:,2:2:end);
            sqzSize=[sqzSize(1) sqzSize(2)/2 2];
        end
        if isjnseq
            sqzDims{end+1}='Set';
        else
            sqzDims{end+1}='Ida';
        end
        
        %Rerun as data has been modified:
        %   Squeeze the data to remove singleton dims
        fids=squeeze(data);
    end
elseif ~isempty(strfind(sequence,'edit_529')); %Is this WIP 529 (MEGA-PRESS)?
    %isWIP529
    if isempty(imaLoopAssignment)
        if isVD || isVE
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Ave','Eco',''};
        else
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Set','Eco',''};
        end
    end
elseif ~isempty(strfind(sequence,'edit_859')); %Is this WIP 859 (MEGA-PRESS)?
    %isWIP859
    if isempty(imaLoopAssignment)
        if isVD || isVE
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Ave','Ide',''};
        else
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Set','Ide',''};
        end
    end
elseif ~isempty(strfind(sequence,'jn_')); %Is this any one of Jamie Near's sequences?
    %isjnseq
    if isempty(imaLoopAssignment)
        if isVD || isVE
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Ave','Set',''};
        else
            imaLoopAssignment = {'t','coils','averages','subSpecs','extras';'Col','Cha','Set','Ida',''};
        end
    end
elseif ~isempty(strfind(sequence,'eja_svs_mpress')); %Is this Eddie Auerbach's MEGA-PRESS?
    %isMinnMP
    if isempty(imaLoopAssignment)
        % Set is the averages dimension for all scanner baselines
        imaLoopAssignment ={'t','coils','averages','subSpecs','extras';'Col','Cha','Set','Eco',''};
    end
elseif ~isempty(strfind(sequence,'svs_se')) ||... %Is this the Siemens PRESS seqeunce?
        ~isempty(strfind(sequence,'svs_st'));    % or the Siemens STEAM sequence?
    %isSiemens
    if isempty(imaLoopAssignment)
        if isVB
            imaLoopAssignment ={'t','coils','averages','subSpecs','extras';'Col','Cha','Set','',''};
        else %VD/E
            imaLoopAssignment ={'t','coils','averages','subSpecs','extras';'Col','Cha','Ave','',''};
        end
    end
    
    %noticed that in the Siemens PRESS and STEAM sequences, there is sometimes
    %an extra dimension containing unwanted reference scans or something.  Remove them here.
    if (isVD || isVE) && strcmp(sqzDims{end},'Phs')
        sqzDims=sqzDims(1:end-1);
        sqzSize=sqzSize(1:end-1);
        if ndims(fids)==4
            fids=fids(:,:,:,2);
            fids=squeeze(fids);
        elseif ndims(fids)==3
            fids=fids(:,:,2);
            fids=squeeze(fids);
        elseif ndims(fids)==2
            fids=fids(:,2);
            fids=squeeze(fids);
        end
    end
    
else
    % Default behaviour
    if isempty(imaLoopAssignment)
        if isVB
            imaLoopAssignment ={'t','coils','averages','subSpecs','extras';'Col','Cha','Set','',''};
        else %VD/E
            imaLoopAssignment ={'t','coils','averages','subSpecs','extras';'Col','Cha','Ave','',''};
        end
    end
end

% Interpret imaLoopAssignment into a struct for more readible code below
for iDx = 1:size(imaLoopAssignment,2)
    imaStruct.(imaLoopAssignment{1,iDx}) = imaLoopAssignment{2,iDx};
end

%Make a pulse sequence identifier for the header (out.seq);
seq=sequence;

%Find the magnetic field strength:
Bo=twix_obj.hdr.Dicom.flMagneticFieldStrength;

%Find the number of averages:
Naverages=twix_obj.hdr.Meas.Averages;

%Find out if multiple coil elements were used:
Ncoils=twix_obj.hdr.Meas.iMaxNoOfRxChannels;  

%Find the TE:
TE = twix_obj.hdr.MeasYaps.alTE{1};  %Franck Lamberton

%Find the TR:
TR = twix_obj.hdr.MeasYaps.alTR{1};  %Franck Lamberton

%Now begin indexing the dimensions of the data array. ie. create the dims
%structure, which specifies which dimensions of the data array are being
%used to hold the time-domain data, the multiple coil channels, the
%average, the sub-spectra, and any additional dimensions.
sqzDims_update=sqzDims;
dimsToIndex=[1:length(sqzDims)];

%First index the dimension of the time-domain data
dims.t=find(strcmp(sqzDims,imaStruct.t));
if ~isempty(dims.t)
    %remove the time dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.t);
else
    dims.t=0;
    error('ERROR:  Spectrom contains no time domain information!!');
end

%Now index the dimension of the coil channels
dims.coils=find(strcmp(sqzDims,imaStruct.coils));
if ~isempty(dims.coils)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.coils);
else
    dims.coils=0;
end

%Now index the dimension of the averages
dims.averages=find(strcmp(sqzDims,imaStruct.averages));
if ~isempty(dims.averages)
    %remove the averages dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
else
    %If no Averages dimension was found, then check for a "Repetitions"
    %dimension.  If that is found, store it under "averages".  If both
    %"Averages" and "Repetitions" dimensions are found, "Repetitions" will
    %be indexed under "Extras", since "Repetitions is not currently an
    %option in FID-A.
    dims.averages=find(strcmp(sqzDims,'Rep'));
    if ~isempty(dims.averages)
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
    else
        %If neither an "Averages" or a "Repetitions" dimension is found,
        %then set the FID-A "Averages" dimension to zero.
        dims.averages=0;
    end
end

%Now we have indexed the dimensions containing the timepoints, the coil
%channels, and the averages.  As we indexed each dimension, we removed the
%corresponding index from the dimsToIndex vector.  At this point, if there
%are any values left in the dimsToIndex vector, then there must be some
%additional dimensions that need indexing. 
if ~isempty(dimsToIndex) 
    dims.subSpecs=find(strcmp(sqzDims,imaStruct.subSpecs));
    if isempty(dims.subSpecs)       
        % If nothing is prescribed then take the next dimension in the list
        dims.subSpecs=dimsToIndex(1);
    end
     %remove the sub-spectra dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.subSpecs);    
else
    dims.subSpecs=0;
end

%And if any further dimensions exist after indexing the sub-spectra, call
%these the 'extras' dimension. WTC. This doesn't handle the case where
%there are more than one more dimensions to index through.
if ~isempty(dimsToIndex)
    dims.extras=find(strcmp(sqzDims,imaStruct.extras));
    if isempty(dims.extras)       
        % If nothing is prescribed then take the next dimension in the list
        dims.extras=dimsToIndex(1);
    end
     %remove the sub-spectra dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.extras);    
else
    dims.extras=0;
end

%Now that we've indexed the dimensions of the data array, we now need to
%permute it so that the order of the dimensions is standardized:  we want
%the order to be as follows:  
%   1) time domain data.  
%   2) coils.
%   3) averages.
%   4) subSpecs.
%   5) extras.
permuteVector = [dims.t dims.coils dims.averages dims.subSpecs dims.extras];

% At the same time update the dims index to reflect the new order. This
% code works on the basis that the order is set in the line above. The call
% to fieldnames might not be that robust as I don't know what guarentees
% the order of the elements.
fn = fieldnames(dims);
count = 0;
for iDx = 1:numel(fn)
    if dims.(fn{iDx}) ~= 0
        count = count + 1;
        dims.(fn{iDx}) = count;
    end
end
permuteVector(permuteVector==0) = [];
fids=permute(fids,permuteVector);

%Now get the size of the data array:
sz=size(fids);

%Now take fft of time domain to get fid:
specs=fftshift(ifft(fids,[],dims.t),dims.t);
    

%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;  %Franck Lamberton
spectralwidth=1/dwelltime;
    
%Get TxFrq
txfrq=twix_obj.hdr.Meas.Frequency;

%Get Date
%date = getfield(regexp(twix_obj.hdr.MeasYaps.tReferenceImage0, ...
%'^".*\.(?<DATE>\d{8})\d*"$', 'names'), 'DATE');  %Franck Lamberton

date=''; %The above code for extracting the date from the header 
         %was causing problems.  Since date is not critical
         %for almost any applications, removing it now to be fixed at a
         %later date.

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    if dims.averages~=0
        averages=sz(dims.averages)*sz(dims.subSpecs);
        rawAverages=averages;
    else
        averages=sz(dims.subSpecs);
        rawAverages=1;
    end
else
    if dims.averages~=0
        averages=sz(dims.averages);
        rawAverages=averages;
    else
        averages=1;
        rawAverages=1;
    end
end

%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired 
%subspectra in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
ppm=ppm+4.65;

t=[0:dwelltime:(sz(1)-1)*dwelltime];


%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=seq;
out.te=TE/1000;
out.tr=TR/1000;
out.pointsToLeftshift=twix_obj.image.freeParam(1);


%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
if out.averages == 1
    out.flags.averaged=1;
else
    out.flags.averaged= 0;
end
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end



%DONE

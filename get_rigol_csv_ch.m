function dat = get_rigol_csv_ch(filename,varargin); 
% function dat = get_rigol_csv_ch(filename,varargin)
% 
% estrarre dati da un file filename.csv generato dall'oscilloscopio 
% Rigol MS02102A
% format [t, {V}] (unica matrice con tutti i dati, con t in colonna 1]
%    se presente 1 canale [t,v], se sono presenti 2 canali [t, v1, v2]
%
% opzione 'no_talk' spegne report a video
%
% wjwiv 20201020

talk = true; 

if length(varargin)>0
    for jj=1:length(varargin)
        if strcmp(varargin{jj},'no_talk')
            talk = false;
        end
    end
end

dat = csvread([filename '.csv'],2);
N_points = size(dat,1); 
if talk 
    disp(['load ' filename ', dati per ' int2str(N_points) ' punti']);
end

% capire che dati ci sono .. 
fid = fopen ([filename '.csv']);
if fid ~= -1
    % get first line of header
    s = fgetl(fid);
    s1 = textscan(s,'%s%s%s%s%s%s','delimiter',',');
    fclose ( fid );
end

CH1 = false; 
CH2 = false; 
N_ch = 0;

for jj=1:6
    if strcmp(s1{jj},'CH1') 
        CH1 = true; N_ch = N_ch +1;
        if talk
            disp([' dati per CH1']);
        end
    end
    if strcmp(s1{jj},'CH2') 
        CH2 = true; N_ch = N_ch +1;
        if talk
            disp([' dati per CH2']);
        end
    end
    if strcmp(s1{jj},'Start') 
        t0_col = jj; 
    end
end    

tdat = csvread([filename '.csv'],1,t0_col-1,[1, t0_col-1, 1, t0_col]);
toff = tdat(1); 
DT = tdat(2); 

dat(:,1) = dat(:,1)*DT + toff;
Ncol = size(dat,2); 

% levare ultima colonna (se solo zero)
if sum(dat(:,Ncol) == zeros(size(dat,1),1)) == size(dat,1)
    dat = dat(:,1:Ncol-1);
    if talk 
        disp(' rimuovere ultima colonna di zero');
    end
end


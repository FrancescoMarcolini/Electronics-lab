function [fig_nums] = bode_plot(f,H,varargin)
% function fig_nums = bode_plot(f,H,varargin)
% 
% funzione per generare un diagramma di Bode (modulo e fase)
%
% f: vettore di frequenze (colonna)
% H: vettore (colonna) funzione di trasferimento complesso (oppure due vettore, uno
%   ampiezza e uno con fase in radianti, [mod(H) angle(H)])
% fig_nums: numero figura + assi, necessari per aggiungere curve allo
%       stesso grafico
% 
% options (varargin):
%   'fig'     prossimo argomento  [fig_num axes_num_modulo axes_num_fase]
%           generato da precedente chiamata a bode_plot
%   'add_dB'  aggiungere scala dB al modulo
%   'Ohm'     unita' di misura  Ohm
%   'col'     prox arg e' il colore ('b', 'r', 'g', 'k', ecc)
%   'ylim'    prox arg e' limiti di scala y in modulo (es [1e-3 1])
%   'points'  usare punti singoli, non una linea (ad esempio per dati
%           sperimentali), prox arg marker ('.','o','*',ecc)
%   'err'     error bars ... solo per input [mod(H) angle(H)], prox argument
%           sono gli error bars [d(mod(H)) d(angle(H))]
% 
% wjwiv 20/3/2015
% 
% 9/2/2016
% modificato per gestire diverse versione di matlab (problema con la
% grafica con gli errori), grazie Mauro Hueller!

col = 'b';
f_label = 'f (Hz)';
pos_amp = [.15 .45 .7 .47];
pos_phase = [.15 .14 .7 .28];
add_dB = false;
lwid = 2;
mark = 'none';
msize = 20; 
line = '-';
use_errs = false;
amp_lab = '|H|';

if length(varargin)>0
    for jj=1:length(varargin)
        if strcmp(varargin{jj},'fig')
            fig_nums = varargin{jj+1};
            fig_num = fig_nums(1);
            amp_ax = fig_nums(2);
            phase_ax = fig_nums(3);
        end
        if strcmp(varargin{jj},'f_norm')
            if strcmp(varargin{jj+1},'3dB')
                f_label = ['f / f_{3dB}'];
            else 
                f_label = ['f \times ' varargin{jj+1}];
            end
        end
        if strcmp(varargin{jj},'add_dB')
            add_dB = true;
        end
        if strcmp(varargin{jj},'col')
            col = varargin{jj+1};
        end
        if strcmp(varargin{jj},'ylim')
            ylims=varargin{jj+1};
        end
        if strcmp(varargin{jj},'points')
            mark = varargin{jj+1};
            line = 'none';
        end
        if strcmp(varargin{jj},'err')
            use_errs = true;
            errs = varargin{jj+1};
            dH_mod = errs(:,1);
            dH_phase = errs(:,2);
        end
        if ( strcmp(varargin{jj},'ohm') | strcmp(varargin{jj},'Ohm') )
            amp_lab = 'Z (\Omega)';
        end
    end
end

if size(H,2) == 1
    mod_H = abs(H);
    phase_H = angle(H);
elseif size(H,2) == 2
    mod_H = H(:,1);
    phase_H = H(:,2);
end

if exist('fig_num')
    figure(fig_num);
else
    fig_num=figure;
    amp_ax=axes('position',pos_amp);
    phase_ax=axes('position',pos_phase);    
end

axes(amp_ax);
if use_errs
    erp = errorbar(f,mod_H,dH_mod);
    set(erp,'color',col,'linewidth',lwid,'linestyle',line,'marker',mark);
    set(gca,'xscale','log','yscale','log');
    setbarwidth(erp,f,0);
    if strcmp(mark,'.')
        set(erp, 'markersize',msize);
    end
else
    lin = loglog(f,mod_H,col,'linewidth',lwid,'linestyle',line,'marker',mark);
    if strcmp(mark,'.')
        set(lin, 'markersize',msize);
    end
end
hold on;
grid on;
set(amp_ax,'xticklabel','');
ylabel(amp_lab);

if exist('ylims')
    ylim(ylims);
else
    ylim([0.99*10^floor(log10(min(mod_H)*1.01)) 1.01*10^ceil(log10(max(mod_H)*.99))]);
end

xl = [10^floor(log10(min(abs(f)*1.01))) 10^ceil(log10(max(abs(f)*.99)))]; 
xlim(xl); 

if add_dB
    % check if there is an existing dB plot that needs to die ... 
    child = get(gcf,'Children'); 
    bad_child = child (child ~= amp_ax & child~= phase_ax);
    if length(bad_child) ==1 
        delete(bad_child);
    end
%    child = get(gcf,'Children')
    yl=get(amp_ax,'ylim');
    ytick = get(amp_ax,'ytick');
    ax_dB = axes('position',pos_amp);
    ylim_dB = [20*log10(yl(1)) 20*log10(yl(2))];
    set(ax_dB,'xtick',[],'ylim',[20*log10(yl(1)) 20*log10(yl(2))],'ylim',ylim_dB, ... 
        'ytick',20*log10(ytick),'color','none','yaxislocation','right');
    ylabel('|H| (dB)');
%    child = get(gcf,'Children')
    axes(amp_ax);

end


axes(phase_ax);
ang_deg = phase_H*180/pi;
if use_errs
    erp = errorbar(f,ang_deg,dH_phase*180/pi);
    set(erp,'color',col,'linewidth',lwid,'linestyle',line,'marker',mark);
    set(gca,'xscale','log');
    setbarwidth(erp,f,0);
    %disp('hello');    
    if strcmp(mark,'.')
        set(erp, 'markersize',msize);
    end
else
    sem  = semilogx(f,ang_deg,col,'linewidth',lwid,'linestyle',line,'marker',mark);
    if strcmp(mark,'.')
        set(sem, 'markersize',msize);
    end
end
hold on;
grid on;
xlabel(f_label);
ylabel('\phi (deg)');
xlim(xl);
if ~(max(ang_deg) == min(ang_deg))
    ylim([10*floor(min(ang_deg/10))-5 , 10*ceil(max(ang_deg/10))+5]);
end
hold on;

fig_nums=[fig_num amp_ax phase_ax];









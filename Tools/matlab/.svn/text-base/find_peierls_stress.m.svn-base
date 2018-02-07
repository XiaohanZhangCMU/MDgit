% find dislocation velocity from cfg files

if (~exist('tilt'))
%tilt = '48.53';
%tilt = '45.29';
tilt = '-90';
%tilt = '-29.5';
%tilt = '-11.42';
end

if (~exist('Nexpan'))
%Nexpan = 1;
Nexpan = 8;
end

if (~exist('runsdir'))
    %runsdir = '/mnt/MyBook/Meso_Backup/kwkang/Codes/MD++.07-28-2007/runs/ta-peierls/convergetest/YesPR0-tiltBox/';
    %runsdir = '/mnt/MyBook/Meso_Backup/kwkang/Codes/MD++.07-28-2007/runs/ta-peierls/harold_data/';
    runsdir = '/mnt/MyBook/Meso_Backup/kwkang/Codes/MD++.07-28-2007/runs/ta-peierls/wcr_data/';
    %runsdir = '../../../Codes/MD++.svn/runs/su-ahpcrc/';
end

if (~exist('folder_format'))
    folder_format = 'ta-peierls-tilt%s-expan%d-zxcgr/';
    %folder_format = 'ta-peierls-tilt%s-expan%d-zxcgr-s/';
end    

if (~exist('crit_dist'))
    crit_dist = 8;
end

if (~exist('crit_dist_back'))
    crit_dist_back = 5;
end

if (~exist('crit_ddr'))
    crit_ddr = 0.1;
    if strcmp(tilt,'0')
        crit_ddr = 1.2;
    end
end

if (~exist('figs'))
    figs=1;
end

if (~exist('octave'))
    octave=0;
end

if (~exist('write_html'))
    write_html=0;
end

if (~exist('html_folder'))
    html_folder = 'html';
end    
  
if (write_html)
    figs=1;
end


font_size = 18;

if (write_html)
    output_folder = sprintf('%s/tilt%s',html_folder,tilt);
    [status,message,messageid] = mkdir(output_folder);
end

peierlsstress = zeros(Nexpan,1);
for expan=1:Nexpan,
    cfgdir = sprintf(folder_format,tilt,expan);
    %[runsdir,cfgdir,'ta-stress*.cfg']
    cfgfiles = dir([runsdir,cfgdir,'ta-stress*.cfg*']);
    nfiles = length(cfgfiles);
    if (nfiles == 0)
        disp(sprintf('No cfg files in %s',cfgdir))
        break
    end

    % Load EPOT.dat
    epotdata = load([runsdir,cfgdir,'EPOT.dat']);
    stress = epotdata(2:end,1); dr = epotdata(2:end,4);
    [stress, idx] = sort(stress); dr = dr(idx);
    ddr = dr(2:end)-dr(1:end-1);
    % choose stress range to analyze cfg files
    id = find(ddr > crit_ddr) + 1;
    if(isempty(id))
        id1 = length(dr); id0 = 1;
    else
        id0 = max(id(1)-4,1); id1 = min(id(end),id0+10); 
    end
    
    % Find original position of dislocation 
    filename = [runsdir,cfgdir,'ta-PRrlx.cfg']; [np,data,h]=loadcfg(filename);
    % extract core atoms
    csd = data(:,4); ind = find(csd>0.5 & csd<8 );
    
    % make sure all atoms are nearest images
    s = data(ind,1:3); s(:,1)=s(:,1)-round(s(:,1)-s(1,1)); s(:,3)=s(:,3)-round(s(:,3)-0.5); r=s*h';
    % compute average dislocation position
    %x0 = mean(s(:,1)*h(1,1)); 
    x0 = mean(r(:,1));
    
    % plot dislocation position
    if (figs && id1-id0+2 < 15)
        figure(1); clf; subplot(1,id1-id0+2,1); plot(r(:,3),r(:,1),'ro'); %axis equal    
        xlim([0 h(3,3)]); ylim([0 h(1,1)]); 
        %set(gca,'FontSize',font_size); xlabel('z'), ylabel('x'), title('0');
        %xlabel('z','fontsize',font_size), ylabel('x','fontsize',font_size) , title('0','fontsize',font_size);
        set(gca,'xtickmode','manual'), set(gca,'xticklabelmode','manual')
        %set(gca,'xtick',[0 h(3,3)]), set(gca,'xticklabel',{'0','L_z'})
        xlabel('z'), title('0'); 
    end
    xavg = zeros(id1-id0+1,1); 
    for id=id0:id1,
        stress_id = stress(id);
        % load atom configuration at stress = id
        filename = [runsdir,cfgdir,sprintf('ta-stress%d.cfg',stress_id)]; 
        [np,data,h]=loadcfg(filename);

        % extract core atoms
        csd = data(:,4); ind = find(csd>0.5 & csd<8 ); 
        
        if (isempty(ind))
             % Screw dislocation has cross slip and escapes the cell. 
             xavg(id-id0+1) = x0+100;
        else            
            % make sure all atoms are nearest images
            s = data(ind,1:3); s(:,1)=s(:,1)-round(s(:,1)-s(1,1));  s(:,3)=s(:,3)-round(s(:,3)-0.5); r=s*h';
            % compute average dislocation position
            %xavg(id-id0+1) = mean(s(:,1)*h(1,1)); 
            xavg(id-id0+1) = mean(r(:,1));

            % plot dislocation position
            if (figs && id1-id0+2 < 15)
                figure(1); subplot(1,id1-id0+2,id-id0+2); p1=plot(r(:,3),r(:,1),'ro'); %axis equal
                xlim([0 h(3,3)]); ylim([0 h(1,1)]); set(gca,'YTickLabel',[]);
                %set(gca,'FontSize',font_size); xlabel('z'); title(num2str(stress_id));
                xlabel('z','fontsize',font_size), title(num2str(stress_id),'fontsize',font_size);
                set(gca,'xtickmode','manual'), set(gca,'xticklabelmode','manual')
                %set(gca,'xtick',[0 h(3,3)]), set(gca,'xticklabel',{'0','L_z'})
            end                
        end
    end

    if(figs)
        figure(2),
        subplot(2,1,1), p2=plot(stress(id0:id1), xavg, 'o-'); 
        set(p2,'MarkerSize',10); set(p2,'MarkerFaceColor','b'); set(p2,'LineWidth',2); ylim([0 h(1,1)]); 
        %set(gca,'FontSize',font_size); xlabel('stress (MPa)'); ylabel('x (Angstrom)'); legend('cfg',2);
        xlabel('stress (MPa)','FontSize',font_size); ylabel('x (Angstrom)','FontSize',font_size); legend('cfg',2);
        subplot(2,1,2), p2=plot(stress(id0:id1), dr(id0:id1), 'o-');
        set(p2,'MarkerSize',10); set(p2,'MarkerFaceColor','b'); set(p2,'LineWidth',2); 
        %set(gca,'FontSize',font_size); xlabel('stress (MPa)'); ylabel('dr');
        xlabel('stress (MPa)','FontSize',font_size); ylabel('dr','FontSize',font_size);
        legend('EPOT.dat',2);
    end

    % decide whether dislocation has moved by a significant distance
    %id_ps = find((xavg - x0 > crit_dist)+(xavg - x0 < -crit_dist_back),1,'first'); 
    %ep_ps = find(dr > 1, 1, 'first');
    id_ps = find((xavg - x0 > crit_dist)+(xavg - x0 < -crit_dist_back)); if(~isempty(id_ps)),id_ps = id_ps(1); end
    ep_ps = find(dr > 1); if(~isempty(ep_ps)), ep_ps = ep_ps(1); end
    
    if (~isempty(id_ps))
        % Peierls stress has been found
        ps = stress(id_ps+id0-1);
        disp(sprintf('tilt = %s expan = %d %d-%d Peierls stress = %d (MPa)',...
            tilt,expan,stress(id0),stress(id1),ps))
        if(figs)
            if (id1-id0+2 < 15)
            figure(1);
            subplot(1,id1-id0+2,id_ps+1); t1=title(num2str(ps)); 
            set(t1,'Color','r'); set(t1,'FontName','bold');
            end
            figure(2);
            subplot(2,1,1); hold on; p1=plot(ps, xavg(id_ps), 'ro');
            set(p1,'MarkerSize',16); set(p1,'LineWidth',3); hold off
        end
        peierlsstress(expan) = ps;        
    elseif (~isempty(ep_ps))
        % if dr in EPOT.dat is too large, we also declare dislocation has moved
        ps = stress(ep_ps);
        disp(sprintf('tilt = %s expan = %d %d-%d Peierls stress = %d (MPa)',...
            tilt,expan,stress(id0),stress(id1),ps))
        if(figs)
            figure(2);
            subplot(2,1,2); hold on; p2=plot(ps, dr(ep_ps), 'ro');
            set(p2,'MarkerSize',8); set(p2,'LineWidth',3); hold off
        end
        peierlsstress(expan) = ps;
    else
        disp(sprintf('Peierls stress for tilt = %s expan = %d is not found!',tilt,expan));
        continue;
    end
    
    if (write_html)
        if (id1-id0+2<15)
            figure(1); 
            plot_name_config = sprintf('%s/config_tilt%s_expan%d',output_folder,tilt,expan);
            if (~octave)
                print('-djpeg',[plot_name_config,'.jpg']);
            else
                print([plot_name_config,'.eps'],'-depsc');
                system(['convert -antialias -density 96x96 -enhance -quality 100 ',plot_name_config,'.eps ',plot_name_config,'.jpg']);
            end
        end
        figure(2);
        plot_name_pos = sprintf('%s/pos_tilt%s_expan%d',output_folder,tilt,expan);
        if (~octave)
            print('-djpeg',[plot_name_pos,'.jpg']);
        else
            print([plot_name_pos,'.eps'],'-depsc');
            system(['convert -antialias -density 96x96 -enhance -quality 100 ',plot_name_pos,'.eps ',plot_name_pos,'.jpg']);
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%% PS ERROR ESTIMATE %%%%%%%%%%%%%%%%%%%%%%%%
%
% Round-off err : Unbounded. Do not want to include this error,
% unless the image error is bigger than the round-off err.
%
% The trend of round-off err is (2MPa increase)/100,000atoms.
% Below 10^5 atoms, round-off err is less than 2 MPa.
% In this range, the image err determines the net err.
% If N > 10^5, both err contribute together. The round-off err
% can be estimated as 2MPa/10^5*N. If the round-off err is
% subtracted from the net err, we can get the image err.

% 1. Choose ps_final = ps(i) where N(i) ~ 10^5.
%    If ps(i) == 0, choose i such that N(i) is the 2nd nearest to 10^5.
%
% 2  If i == Nexpan, err = dps_L (purely image err)
%    where dps_L = |ps(i-1) - ps(i)|
%    If ps(i-1) == 0; dps_L = |ps(i-2) - ps(i)|/2
%
% 3. IF i < Nexpan, calculate differences with neighboring PS, 
%      dps_R = |ps(i+1) - ps(i)|.
%    Then, the error is 
%      err = min(err_img, err_num)
%    where err_num = 2/10^5*N(i+1)*(N(i+1)>10^5)
%          err_img = max(dps_R - err_num,0).
%
%    Note that the image error err_img is estimated by
%    subtracting the round-off err from the net error. And we
%    expect that the image error can not be less than zero.
%    In this criteria, err will not choose the round-off err,
%    unless round-off err is less than the image err.
%    If ps(i+1) == 0; Do the simliar error estimate with ps(i-1).
%  
% 5. Finally, err = max(err, DEFAULT_LOWER_BOUND_ERR)


% estimate PS from an array of PS values as a function of cell size
%ind2=find(np_array<=1e5); ind2 = ind2(end); ps_final = peierlsstress(ind2);
%if(ps_final==0 && ind2<Nexpan), ind2=ind2+1; ps_final = peierlsstress(ind2); end

[dN,idx_N]=sort(abs(np_array - 1e5));
ind2 = idx_N(1); ps_final = peierlsstress(ind2);
if(ps_final==0), ind2 = idx_N(2); ps_final=peierlsstress(ind2); end

% estimate PS error size
if(ps_final==0)
  % no need to find error
else

  if (ind2 == Nexpan)
    ps_L = peierlsstress(ind2-1);
    if (ps_L ~= 0)
      dps_L = abs(ps_L - ps_final);
    else
      dps_L = abs(peierlsstress(ind2-2) - ps_final)/2;
    end
    ps_error = dps_L;
  else
    ps_R = peierlsstress(ind2+1);
    if (ps_R ~= 0)
      dps_R = abs(ps_R - ps_final);
      err_num = 2/10^5*np_array(ind2+1)*(np_array(ind2+1)>10^5);
      err_img = max(dps_R - err_num,0);
      ps_error = ceil(min(err_img, err_num));
    elseif (ind2>1)
      dps_L = abs(peierlsstress(ind2-1) - ps_final);
      err_num = 2/10^5*np_array(ind2-1)*(np_array(ind2-1)>10^5);
      err_img = max(dps_L - err_num,0);
      ps_error = ceil(min(err_img, err_num));
    end
  end
end

ps_error = max(ps_error, DEFAULT_LOWER_BOUND_ERR);


if(figs)
    figure(3)
    p3=plot([1:Nexpan], peierlsstress,'bo-'); 
    set(p3,'MarkerSize',8); set(p3,'MarkerFaceColor','m'); set(p3,'LineWidth',2);
    %set(gca,'FontSize',font_size); xlabel('expan'); ylabel('Peierls stress (MPa)'); 
    xlabel('expan','FontSize',font_size); ylabel('Peierls stress (MPa)','FontSize',font_size); 
    ylim([0 max(peierlsstress)+(ps_error+1)]);
    
    figure(4)
    p4=plot(np_array, peierlsstress,'bo-'); hold on
    p4_2=errorbar(np_array(ind2), ps_final, ps_error); 
    set(p4,'MarkerSize',6); set(p4,'MarkerFaceColor','b'); set(p4,'LineWidth',2);     
    set(p4_2,'LineWidth',2); set(p4_2,'Color','r')
    %set(p4(2),'MarkerSize',8); set(p4(2),'MarkerFaceColor','m'); set(p4(2),'LineWidth',2);
    %set(gca,'FontSize',font_size); xlabel('expan'); ylabel('Peierls stress (MPa)'); 
    xlabel('NP0','FontSize',font_size); ylabel('Peierls stress (MPa)','FontSize',font_size); 
    ylim([0 max(peierlsstress)+(ps_error+1)]), hold off
end

%writing html files
if (write_html)
  html_name_config = sprintf('%s/config_tilt%s.html',output_folder,tilt);
  fid = fopen(html_name_config,'w');
  fprintf(fid,'<html>\n');
  fprintf(fid,'<title> dislocation core at tilt%s </title>\n',tilt);
  fprintf(fid,'<table border=1 width="900">\n');
  fprintf(fid,'<tbody>\n');
  for expan=1:Nexpan,
      fprintf(fid,'  <tr>\n');
      fprintf(fid,'    <td width="450">\n');
      fprintf(fid,'       <a name="expan%d">\n',expan);
      fprintf(fid,'       <font color=0x000000> tilt%s_expan%d </font>&nbsp;&nbsp;&nbsp;&nbsp;',tilt,expan);
      fprintf(fid,'       </a>\n');
      fprintf(fid,'       Peierls stress = %d MPa &nbsp;&nbsp;&nbsp;&nbsp;',peierlsstress(expan));
      fprintf(fid,'       <a href="index.html">view convergence</a><br>');
      fprintf(fid,'       <img src="pos_tilt%s_expan%d.jpg" width=400><hr>\n',tilt,expan);
      fprintf(fid,'    </td>\n');
      fprintf(fid,'    <td>\n');
      fprintf(fid,'       <img src="config_tilt%s_expan%d.jpg" width=400> \n',tilt,expan);
      fprintf(fid,'    </td>\n');
      fprintf(fid,'  </tr>\n');
  end
  fprintf(fid,'</tbody>\n');
  fprintf(fid,'</table>\n');
  fprintf(fid,'</html>\n');
  fclose(fid);
  
  figure(3);
  plot_name_psconv = sprintf('%s/ps_conv_tilt%s',output_folder,tilt);
  if (~octave)
      print('-djpeg',[plot_name_psconv,'.jpg']);
  else
      print([plot_name_psconv,'.eps'],'-depsc');
      system(['convert -antialias -density 96x96 -enhance -quality 100 ',plot_name_psconv,'.eps ',plot_name_psconv,'.jpg']);
  end
  figure(4);
  plot_name_psconv2 = sprintf('%s/ps_conv2_tilt%s',output_folder,tilt);
  if (~octave)
      print('-djpeg',[plot_name_psconv2,'.jpg']);
  else
      print([plot_name_psconv2,'.eps'],'-depsc');
      system(['convert -antialias -density 96x96 -enhance -quality 100 ',plot_name_psconv2,'.eps ',plot_name_psconv2,'.jpg']);
  end
  html_name_psconv = sprintf('%s/index.html',output_folder);
  fid = fopen(html_name_psconv,'w');
  fprintf(fid,'<html>\n');
  fprintf(fid,'<title> Dislocation Peierls Stress at tilt%s </title>\n',tilt);
  fprintf(fid,'<table border=1 width="600">\n');
  fprintf(fid,'<tbody>\n');
  fprintf(fid,'  <tr>\n');
  fprintf(fid,'    <td width="400">\n');
  fprintf(fid,'<H3>\n');
  fprintf(fid,'<font color=0x000000> Peierls stress at </font>');
  fprintf(fid,'<font color=0x0000FF>%s</font> degree = <font color=0x0000FF>%d MPa </font><br>',tilt,ps_final);
  fprintf(fid,'<img src="ps_conv_tilt%s.jpg" width=400>\n',tilt);
  fprintf(fid,'</H4>\n');
  fprintf(fid,'    </td>\n');
  fprintf(fid,'    <td>\n');
  fprintf(fid,'      <center>\n');
  fprintf(fid,'    go to <a href="../index.html">parent directory</a><br><br>\n');
  fprintf(fid,'    view dislocation position at <br><br>');
  for expan=1:Nexpan,
      fprintf(fid,'<a href="config_tilt%s.html#expan%d">expan%d</a><br>',tilt,expan,expan);
  end
  fprintf(fid,'      </center>\n');
  fprintf(fid,'    </td>\n');
  fprintf(fid,'  </tr>\n');
  fprintf(fid,'  <tr>\n');
  fprintf(fid,'    <td width="400">\n');
  fprintf(fid,'<img src="ps_conv2_tilt%s.jpg" width=400>\n',tilt);
  fprintf(fid,'    </td>\n');
  fprintf(fid,'    <td>\n');
  fprintf(fid,'      <center>\n');
  fprintf(fid,'      view <a href="ps_np.dat">data file</a>\n');
  fprintf(fid,'    </td>\n');
  fprintf(fid,'  </tr>\n');
  fprintf(fid,'</tbody>\n');
  fprintf(fid,'</table>\n');
  fprintf(fid,'</html>\n');
  fclose(fid);  

  % save Peierls stress data into file
  dat_name_psnp = sprintf('%s/ps_np.dat',output_folder);
  fid = fopen(dat_name_psnp,'w');
  fprintf(fid,'%% expan Peierls Stress (MPa)  NP (perfect crystal)\n');
  for expan=1:Nexpan,
      fprintf(fid,'     %d              %g               %d\n',expan, peierlsstress(expan), np_array(expan));
  end
  fclose(fid);    

end

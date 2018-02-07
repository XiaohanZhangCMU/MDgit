% Find Peierls stress for all the tilt angles

clear all; 

% Basic settings

write_html=1; 
Nexpan = 8;
%runsdir = '/mnt/MyBook/Meso_Backup/kwkang/Codes/MD++.07-28-2007/runs/ta-peierls/wcr_data/';
%runsdir = '/mnt/MyBook/Meso_Backup/kwkang/Codes/MD++.07-28-2007/runs/ta-peierls/su-ahpcrc_data/';
runsdir = 'H:/su-ahpcrc_data/runs/';
folder_format = 'ta-peierls-tilt%s-expan%d-zxcgr/';
%html_folder = 'html';
html_folder = 'D:/My Documents/Tex/A28_dis-core/data/html';
cellsizefile = 'D:/My Documents/Tex/A28_dis-core/data/tiltangle.dat';
octave=0;

DEFAULT_LOWER_BOUND_ERR = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to run this program, you should have a 'cellsize.out' file
% that contains the character angle and supercell geometry
cellsize = load(cellsizefile);
n_tilt_data = length(cellsize(:,1));
for i = 1:n_tilt_data,
    tilt_data{i} = sprintf('%.4g',cellsize(i,1));
    % this is a temporary fix, we shall remove these lines in the future
    if cellsize(i)==22,
        tilt_data{i} = '22.00';
    end
end

ps_data = zeros(n_tilt_data,3); % tilt ps error
for tilt_id=31:31,
%for tilt_id=1:n_tilt_data,

  tilt = tilt_data{tilt_id};
  ps_data(tilt_id,1) = str2num(tilt);
  
  chkdir = dir([runsdir,sprintf(folder_format(1:23),tilt),'*-zxcgr']);
  if isempty(chkdir)
    disp(['Directory ',sprintf(folder_format(1:23),tilt),'*-zxcgr Not Exist'])
    continue;
  end
  
  % calculate np0
  %lattice_i = latticedata(tilt_id,:);
  %cal_np0
  %tilt_structure(tilt_id).np = np_array;

  %ind=find(cellsize(:,1)==str2num(tilt));
  np0=det([cellsize(tilt_id,2:4);cellsize(tilt_id,6:8);cellsize(tilt_id,10:12)])*...
           cellsize(tilt_id,5)*cellsize(tilt_id,9)*cellsize(tilt_id,13) * 2 * 0.75;
  np_array=np0.*[1:Nexpan].^2;
  
  ps_final=0; ps_error=0;
  
  % calculate peierls stress
  find_peierls_stress;
  
  %tilt_structure(tilt_id).ps = peierlsstress';
  ps_data(tilt_id,2:3) = [ps_final ps_error];

end

if (write_html)
  html_name_alltilt = sprintf('%s/index.html',html_folder);
  %tilt_folders = dir([html_folder,'/tilt*']);
  tilt_folders = dir(html_folder);
  fid_mainhtml = fopen(html_name_alltilt,'w');
  fprintf(fid_mainhtml,'<html>\n');
  fprintf(fid_mainhtml,'<title> Dislocation Peierls Stress at different orientations</title>\n');
  fprintf(fid_mainhtml,'<H3> Dislocation Peierls Stress at different orientations</H3>\n');
  fprintf(fid_mainhtml,'<a href="ps_alltilt.jpg"><img src="ps_alltilt.jpg" width=600 border=0 align="top"></a>\n');
  fprintf(fid_mainhtml,'<table border=1 width="300">\n');
  fprintf(fid_mainhtml,'<tbody>\n');
  fprintf(fid_mainhtml,'  <tr>\n');
  fprintf(fid_mainhtml,'    <td>Tilt Angle (degree)</td><td>PS (MPa)</td><td>Error (MPa)</td>\n');
  fprintf(fid_mainhtml,'  </tr>\n');
  for tilt_id=1:n_tilt_data,
    tilt = tilt_data{tilt_id};
  
    fprintf(fid_mainhtml,'  <tr>\n');
    fprintf(fid_mainhtml,'    <td>\n');
    fprintf(fid_mainhtml,'      <center>\n');
    fprintf(fid_mainhtml,'        <a href="tilt%s/index.html">%s</a>',tilt,tilt);
    fprintf(fid_mainhtml,'      </center>\n');
    fprintf(fid_mainhtml,'    </td>\n');
    fprintf(fid_mainhtml,'    <td>%d</td><td>&plusmn;%d</td>\n',ps_data(tilt_id,2),ps_data(tilt_id,3));
    fprintf(fid_mainhtml,'  </tr>\n');
  end
  fprintf(fid_mainhtml,'</tbody>\n');
  fprintf(fid_mainhtml,'</table>\n');
  fprintf(fid_mainhtml,'</html>\n');
  fclose(fid_mainhtml);  
end



if (~octave)
    save ps_data.mat ps_data
else
    save -7 ps_data.mat ps_data
end

if (~octave)
    figure(5);
    %p1=stem(data1(:,1), data1(:,end), 'bo'); hold on
    p5=stem(ps_data(:,1),ps_data(:,2),'ro'); hold off
    %set(p1,'MarkerSize',3); set(p1,'Color',[0.8 0.8 0.8]);
    set(p5,'MarkerSize',3); set(p5,'MarkerFaceColor','r');
    xlim([-100 100]); ylim([0 180]); 
    t1 = text(48.53, 71, '[3 3 -1]');
    t2 = text(70.53,161, '[1 1 -1]');
    t3 = text(35.26,141, '[1 1 0]');
    t4 = text(12.28,121, '[5 5 3]');
    t5 = text(0.0,  171, '[1 1 1]');
    xlabel('\theta (degree)');
    ylabel('Peierls stress (MPa)');
    set(gcf,'PaperPosition', [1 1 10 6]);
    plot_name_psalltilt = sprintf('%s/ps_alltilt.jpg',html_folder);
    print('-djpeg',plot_name_psalltilt);
end

%nd = 52;
%ps_data = zeros(nd,2);
%for i=1:nd,
%  tilt = tilt_data{i};
%  ps_data(i,1) = str2num(tilt);
%  data = load(sprintf('%s/tilt%s/ps_np.dat',html_folder,tilt));
%  peierlsstress = data(:,2); NP=data(:,3); nexpan=length(NP);
%  peierlsstress
%  ind2=find(NP<=1e5,1,'last'); ps_final=peierlsstress(ind2);
%  if(ps_final==0 && ind2<nexpan) ind2=ind2+1; ps_final=peierlsstress(ind2); end
%  ps_data(i,2) = ps_final;
%end

disp(sprintf('cputime = %g',cputime()));



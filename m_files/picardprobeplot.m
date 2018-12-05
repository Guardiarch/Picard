% picardprobeplot plots data from the probes of the picard program

if exist('cometframe') ~= 0
  % Default to showing probe fields in the comet frame of reference
  cometframe = logical(1);
end

run('inputpicarda1.m')




fontname='utopia';

for ii= 1:Nprobes
  
  load(['outp/Eprobe' num2str(ii,'%4.4i') '.mat'])
  
  figure(10000+ii)
  clf
  set(gcf,'paperpositionmode','auto')
  plot(timestepsEprobe*dt,Eprobe)
  set(gca,'fontname',fontname,'fontsize',14)
  legend('E_x','E_y','E_z')
  grid on
  title(['$x=' num2str(probestruct(ii).rc(1)) ...
         '$\,m, $y=' num2str(probestruct(ii).rc(2)) ...
         '$\,m, $z=' num2str(probestruct(ii).rc(3)) '$\,m'], ...
        'interpreter','latex','fontname',fontname,'fontsize',18)
  
  dd=dir(['outp/Pprobe_iter*' num2str(ii,'%4.4i') '.mat']);
  if length(dd)>0
    load(['outp/' dd(end).name])
    figure(20000+ii)
    clf
    set(gcf,'paperpositionmode','auto')
    for jj=1:length(Pprobe)
      if ~isempty(Pprobe(jj).r)
        subplot(length(Pprobe),4,(jj-1)*4+1)
        plot3(Pprobe(jj).r(:,1),Pprobe(jj).r(:,2),Pprobe(jj).r(:,3),'k.')
        xlabel('$x$','interpreter','latex','fontname',fontname,'fontsize',18)
        ylabel('$y$','interpreter','latex','fontname',fontname,'fontsize',18)
        zlabel('$z$','interpreter','latex','fontname',fontname,'fontsize',18)
        title(['$x=' num2str(probestruct(ii).rc(1)) ...
               '$\,m, $y=' num2str(probestruct(ii).rc(2)) ...
               '$\,m, $z=' num2str(probestruct(ii).rc(3)) '$\,m'], ...
              'interpreter','latex','fontname',fontname,'fontsize',14)
        axis equal
        subplot(length(Pprobe),4,(jj-1)*4+2)
        minv = min(Pprobe(jj).v(:,1));
        minbin = floor(minv*10^(-floor(log10(abs(minv))))) * ...
                 10^(floor(log10(abs(minv))));
        maxv = max(Pprobe(jj).v(:,1));
        maxbin = ceil(maxv*10^(-floor(log10(abs(maxv))))) * ...
                 10^(floor(log10(abs(maxv))));
        bins=minbin+(maxbin-minbin)*(0.5+[0:9])/10;
        hist(Pprobe(jj).v(:,1),bins)
        aa=axis;axis([minbin maxbin aa(3:4)])
        xlabel('$v_{x}$','interpreter','latex', ...
               'fontname',fontname,'fontsize',18)
        subplot(length(Pprobe),4,(jj-1)*4+3)
        minv = min(Pprobe(jj).v(:,2));
        minbin = floor(minv*10^(-floor(log10(abs(minv))))) * ...
                 10^(floor(log10(abs(minv))));
        maxv = max(Pprobe(jj).v(:,2));
        maxbin = ceil(maxv*10^(-floor(log10(abs(maxv))))) * ...
                 10^(floor(log10(abs(maxv))));
        bins=minbin+(maxbin-minbin)*(0.5+[0:9])/10;
        hist(Pprobe(jj).v(:,2),bins)
        aa=axis;axis([minbin maxbin aa(3:4)])
        xlabel('$v_{y}$','interpreter','latex', ...
               'fontname',fontname,'fontsize',18)
        title(['Species ' num2str(jj)], ...
              'interpreter','latex','fontname',fontname,'fontsize',14)
        subplot(length(Pprobe),4,(jj-1)*4+4)
        minv = min(Pprobe(jj).v(:,3));
        minbin = floor(minv*10^(-floor(log10(abs(minv))))) * ...
                 10^(floor(log10(abs(minv))));
        maxv = max(Pprobe(jj).v(:,3));
        maxbin = ceil(maxv*10^(-floor(log10(abs(maxv))))) * ...
                 10^(floor(log10(abs(maxv))));
        bins=minbin+(maxbin-minbin)*(0.5+[0:9])/10;
        hist(Pprobe(jj).v(:,3),bins)
        aa=axis;axis([minbin maxbin aa(3:4)])
        xlabel('$v_{z}$','interpreter','latex', ...
               'fontname',fontname,'fontsize',18)
      end
    end
  end
end


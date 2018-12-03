% picardprobeplot plots data from the probes of the picard program

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
        hist(Pprobe(jj).r(:,1))
        xlabel('$v_{x}$','interpreter','latex', ...
               'fontname',fontname,'fontsize',18)
        subplot(length(Pprobe),4,(jj-1)*4+3)
        hist(Pprobe(jj).r(:,2))
        xlabel('$v_{y}$','interpreter','latex', ...
               'fontname',fontname,'fontsize',18)
        title(['Species ' num2str(jj)], ...
              'interpreter','latex','fontname',fontname,'fontsize',14)
        subplot(length(Pprobe),4,(jj-1)*4+4)
        hist(Pprobe(jj).r(:,3))
        xlabel('$v_{z}$','interpreter','latex', ...
               'fontname',fontname,'fontsize',18)
      end
    end
  end
end


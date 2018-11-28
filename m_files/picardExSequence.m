% picardExSequence plots a sequence of Ex components
%
% HG 2018-11-13

run('outp/parameters.m')

fontname = 'utopia';

dd = dir('outp/Efield*.mat');

figure(300)
clf
set(gcf,'paperpositionmode','auto','position',[41 68 1550 900])
clear ha;
cax=zeros(length(dd),2);

for ii=1:length(dd)
  load(['outp/' dd(ii).name])
  ha(ii) = subplot(3,5,ii);
  [a iy0]=min(abs(ycorn-0));
  pp(:,:) = mean(Ex(:,[iy0-1:iy0],:),2);pp=pp.';
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(xcorn,zcorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([xmin xmax zmin zmax])
  shading flat
  xlabel('x','fontname',fontname,'fontsize',18)
  ylabel('z','fontname',fontname,'fontsize',18)
  title(['iteration=' num2str(str2num(dd(ii).name(7:12)))], ...
        'fontname',fontname,'fontsize',14)
  cax(ii,:)=caxis(ha(ii));
end
hapos=get(ha,'position');
caxmin=min(cax(:,1));
caxmax=max(cax(:,2));
for ii=1:length(dd);
  caxis(ha(ii),[caxmin caxmax])
end
for jj=1:3
  for ii=1:5
    splno=(jj-1)*5+ii; % subplot number
    set(ha(splno),'position', ...
                  [(0.06+(ii-1)*0.185) (0.08+(3-jj)*0.33) hapos{splno}(3:4)])
  end
end
hcol=colorbar(ha(end));
set(hcol,'position',[0.94 0.1 0.02 0.85])




if 3==4
  h=get(0,'children')
  ccc = pwd;
  cd ~
  for ii =1:length(h)
    filename = ['fig' num2str(ii,'%0.2d') '.png'];
    print(h(ii),'-r600','-dpng',filename)
  end
  cd(ccc)
end

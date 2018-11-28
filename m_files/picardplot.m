% picardplot plots outputs from the picard program.
%
% HG 2018-11-03

run('inputpicarda1.m')

dd = dir('outp/density*.mat');
load(['outp/' dd(end).name])

fontname = 'utopia';

for ii = 1:Nspecies
  figure(20+ii)
  clf
  set(gcf,'paperpositionmode','auto','position',[5 247-10*ii 1331 420])
  subplot(1,3,1)
  [a ix0]=min(abs(xcorn-0));
  pp(:,:) = log10(mean(particle(ii).density([ix0-1:ix0],:,:),1));pp=pp.';
  pp(isnan(pp) | isinf(pp)) = -15;
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(ycorn,zcorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([ymin ymax zmin zmax])
  shading flat
  xlabel('y','fontname',fontname,'fontsize',18)
  ylabel('z','fontname',fontname,'fontsize',18)
  colorbar
  cax=caxis;
  caxis([5 cax(2)])
  title(['n_{' num2str(ii) '}'],'fontname',fontname,'fontsize',14)
  
  subplot(1,3,2)
  [a iy0]=min(abs(ycorn-0));
  pp(:,:) = log10(mean(particle(ii).density(:,[iy0-1:iy0],:),2));pp=pp.';
  pp(isnan(pp) | isinf(pp)) = -15;
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
  colorbar
  cax=caxis;
  caxis([5 cax(2)])
  title(['n_{' num2str(ii) '}'],'fontname',fontname,'fontsize',14)

  subplot(1,3,3)
  [a iz0]=min(abs(zcorn-0));
  pp(:,:) = log10(mean(particle(ii).density(:,:,[iz0-1:iz0]),3));pp=pp.';
  pp(isnan(pp) | isinf(pp)) = -15;
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(xcorn,ycorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([xmin xmax ymin ymax])
  shading flat
  xlabel('x','fontname',fontname,'fontsize',18)
  ylabel('y','fontname',fontname,'fontsize',18)
  colorbar
  cax=caxis;
  caxis([5 cax(2)])
  title(['n_{' num2str(ii) '}'],'fontname',fontname,'fontsize',14)
end


dd = dir('outp/Efield*.mat');
load(['outp/' dd(end).name])

figure(30)
clf
set(gcf,'paperpositionmode','auto','position',[5 65 1017 620])
subplot(3,3,1)
[a ix0]=min(abs(xcorn-0));
pp(:,:) = mean(Ex([ix0-1:ix0],:,:),1);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(ycorn,zcorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([ymin ymax zmin zmax])
shading flat
xlabel('y','fontname',fontname,'fontsize',18)
ylabel('z','fontname',fontname,'fontsize',18)
colorbar
title('E_{x}','fontname',fontname,'fontsize',14)

subplot(3,3,2)
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
colorbar
title('E_{x}','fontname',fontname,'fontsize',14)

subplot(3,3,3)
[a iz0]=min(abs(zcorn-0));
pp(:,:) = mean(Ex(:,:,[iz0-1:iz0]),3);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(xcorn,ycorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([xmin xmax ymin ymax])
shading flat
xlabel('x','fontname',fontname,'fontsize',18)
ylabel('y','fontname',fontname,'fontsize',18)
colorbar
title('E_{x}','fontname',fontname,'fontsize',14)

subplot(3,3,4)
[a ix0]=min(abs(xcorn-0));
pp(:,:) = mean(Ey([ix0-1:ix0],:,:),1);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(ycorn,zcorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([ymin ymax zmin zmax])
shading flat
xlabel('y','fontname',fontname,'fontsize',18)
ylabel('z','fontname',fontname,'fontsize',18)
colorbar
title('E_{y}','fontname',fontname,'fontsize',14)

subplot(3,3,5)
[a iy0]=min(abs(ycorn-0));
pp(:,:) = mean(Ey(:,[iy0-1:iy0],:),2);pp=pp.';
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
colorbar
title('E_{y}','fontname',fontname,'fontsize',14)

subplot(3,3,6)
[a iz0]=min(abs(zcorn-0));
pp(:,:) = mean(Ey(:,:,[iz0-1:iz0]),3);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(xcorn,ycorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([xmin xmax ymin ymax])
shading flat
xlabel('x','fontname',fontname,'fontsize',18)
ylabel('y','fontname',fontname,'fontsize',18)
colorbar
title('E_{y}','fontname',fontname,'fontsize',14)

subplot(3,3,7)
[a ix0]=min(abs(xcorn-0));
pp(:,:) = mean(Ez([ix0-1:ix0],:,:),1);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(ycorn,zcorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([ymin ymax zmin zmax])
shading flat
xlabel('y','fontname',fontname,'fontsize',18)
ylabel('z','fontname',fontname,'fontsize',18)
colorbar
title('E_{z}','fontname',fontname,'fontsize',14)

subplot(3,3,8)
[a iy0]=min(abs(ycorn-0));
pp(:,:) = mean(Ez(:,[iy0-1:iy0],:),2);pp=pp.';
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
colorbar
title('E_{z}','fontname',fontname,'fontsize',14)

subplot(3,3,9)
[a iz0]=min(abs(zcorn-0));
pp(:,:) = mean(Ez(:,:,[iz0-1:iz0]),3);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(xcorn,ycorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([xmin xmax ymin ymax])
shading flat
xlabel('x','fontname',fontname,'fontsize',18)
ylabel('y','fontname',fontname,'fontsize',18)
colorbar
title('E_{z}','fontname',fontname,'fontsize',14)

dd = dir('outp/potential*.mat');
load(['outp/' dd(end).name])

figure(40)
clf
set(gcf,'paperpositionmode','auto','position',[5 247 1331 420])
subplot(1,3,1)
[a ix0]=min(abs(xcorn-0));
pp(:,:) = mean(U([ix0-1:ix0],:,:),1);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(ycorn,zcorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([ymin ymax zmin zmax])
shading flat
xlabel('y','fontname',fontname,'fontsize',18)
ylabel('z','fontname',fontname,'fontsize',18)
colorbar
title('potential','fontname',fontname,'fontsize',14)

subplot(1,3,2)
[a iy0]=min(abs(ycorn-0));
pp(:,:) = mean(U(:,[iy0-1:iy0],:),2);pp=pp.';
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
colorbar
title('potential','fontname',fontname,'fontsize',14)

subplot(1,3,3)
[a iz0]=min(abs(zcorn-0));
pp(:,:) = mean(U(:,:,[iz0-1:iz0]),3);pp=pp.';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
surf(xcorn,ycorn,pp)
clear pp
set(gca,'fontname',fontname,'fontsize',14)
view(2)
axis square
axis([xmin xmax ymin ymax])
shading flat
xlabel('x','fontname',fontname,'fontsize',18)
ylabel('y','fontname',fontname,'fontsize',18)
colorbar
title('potential','fontname',fontname,'fontsize',14)


for ii = 1:Nspecies
  dd = dir(['outp/flux' num2str(ii,'%2.2i') '*.mat']);
  load(['outp/' dd(end).name])

  figure(50+ii)
  clf
  set(gcf,'paperpositionmode','auto','position',[5 65+ii*10 1017 620])
  subplot(3,3,1)
  [a ix0]=min(abs(xcorn-0));
  pp(:,:) = mean(Fx([ix0-1:ix0],:,:),1);pp=pp.';
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(ycorn,zcorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([ymin ymax zmin zmax])
  shading flat
  xlabel('y','fontname',fontname,'fontsize',18)
  ylabel('z','fontname',fontname,'fontsize',18)
  colorbar
  title(['\Gamma_{x}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,2)
  [a iy0]=min(abs(ycorn-0));
  pp(:,:) = mean(Fx(:,[iy0-1:iy0],:),2);pp=pp.';
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
  colorbar
  title(['\Gamma_{x}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,3)
  [a iz0]=min(abs(zcorn-0));
  pp(:,:) = mean(Fx(:,:,[iz0-1:iz0]),3);pp=pp.';
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(xcorn,ycorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([xmin xmax ymin ymax])
  shading flat
  xlabel('x','fontname',fontname,'fontsize',18)
  ylabel('y','fontname',fontname,'fontsize',18)
  colorbar
  title(['\Gamma_{x}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,4)
  [a ix0]=min(abs(xcorn-0));
  pp(:,:) = mean(Fy([ix0-1:ix0],:,:),1);pp=pp.';
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(ycorn,zcorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([ymin ymax zmin zmax])
  shading flat
  xlabel('y','fontname',fontname,'fontsize',18)
  ylabel('z','fontname',fontname,'fontsize',18)
  colorbar
  title(['\Gamma_{y}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,5)
  [a iy0]=min(abs(ycorn-0));
  pp(:,:) = mean(Fy(:,[iy0-1:iy0],:),2);pp=pp.';
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
  colorbar
  title(['\Gamma_{y}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,6)
  [a iz0]=min(abs(zcorn-0));
  pp(:,:) = mean(Fy(:,:,[iz0-1:iz0]),3);pp=pp.';
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(xcorn,ycorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([xmin xmax ymin ymax])
  shading flat
  xlabel('x','fontname',fontname,'fontsize',18)
  ylabel('y','fontname',fontname,'fontsize',18)
  colorbar
  title(['\Gamma_{y}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,7)
  [a ix0]=min(abs(xcorn-0));
  pp(:,:) = mean(Fz([ix0-1:ix0],:,:),1);pp=pp.';
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(ycorn,zcorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([ymin ymax zmin zmax])
  shading flat
  xlabel('y','fontname',fontname,'fontsize',18)
  ylabel('z','fontname',fontname,'fontsize',18)
  colorbar
  title(['\Gamma_{z}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,8)
  [a iy0]=min(abs(ycorn-0));
  pp(:,:) = mean(Fz(:,[iy0-1:iy0],:),2);pp=pp.';
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
  colorbar
  title(['\Gamma_{z}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

  subplot(3,3,9)
  [a iz0]=min(abs(zcorn-0));
  pp(:,:) = mean(Fz(:,:,[iz0-1:iz0]),3);pp=pp.';
  ss=size(pp);
  pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  surf(xcorn,ycorn,pp)
  clear pp
  set(gca,'fontname',fontname,'fontsize',14)
  view(2)
  axis square
  axis([xmin xmax ymin ymax])
  shading flat
  xlabel('x','fontname',fontname,'fontsize',18)
  ylabel('y','fontname',fontname,'fontsize',18)
  colorbar
  title(['\Gamma_{z}, species ' num2str(ii)],'fontname',fontname,'fontsize',14)

end


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

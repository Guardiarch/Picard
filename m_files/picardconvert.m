% picardconvert converts outputs from the picard program to matlab's file
% format.
%
% HG 2018-11-22

ccc=pwd;
[a0,a1]=system('uname -n'); [a0,a2]=system('ps |grep -i matlab');
dummy=[a1 a2 datestr(now)];
dummy=dummy(double(dummy)~=10); dummy=dummy(double(dummy)~=32);

run('inputpicarda1.m')

Nx=Nx_local*iprocs;
Ny=Ny_local*jprocs;
Nz=Nz_local*kprocs;
dxyz=[(xmax-xmin)/Nx (ymax-ymin)/Ny (zmax-zmin)/Nz];
xcorn = xmin + [0:Nx]*dxyz(1); x = 0.5*(xcorn(1:end-1)+xcorn(2:end));
ycorn = ymin + [0:Ny]*dxyz(2); y = 0.5*(ycorn(1:end-1)+ycorn(2:end));
zcorn = zmin + [0:Nz]*dxyz(3); z = 0.5*(zcorn(1:end-1)+zcorn(2:end));
Nprocs = iprocs*jprocs*kprocs;

% Get probe data and put it in a structure
probestruct = struct();
fid=fopen('inputpicarda1.m');
for ii=1:Nprobes
  rprobe=0;
    theline=fgetl(fid);
  if length(theline)<5,
    theline=[theline '     '];
  end
  while ~(strcmp(theline(1:5),'%PROB') | strcmp(theline(1:5),'%prob'))
    theline=fgetl(fid);
    if length(theline)<5,
      theline=[theline '     '];
    end
  end
  while ~(strcmp(theline(1:4),'%END') | strcmp(theline(1:4),'%end'))
    eval(theline)
    theline=fgetl(fid);
    if length(theline)<4,
      theline=[theline '    '];
    end
  end
  probestruct(ii).xc=xc; % probe positions specified in input file
  probestruct(ii).yc=yc;
  probestruct(ii).zc=zc;
  [a b]=min(abs(x-xc));probestruct(ii).rc(1)=x(b); % real probe position
  [a b]=min(abs(y-yc));probestruct(ii).rc(2)=x(b);
  [a b]=min(abs(z-zc));probestruct(ii).rc(3)=x(b);
  probestruct(ii).rprobe=rprobe;
end

% Density
% Prevent two processes from performing simultaneous conversions
if exist([pwd '/outp/datfiles/lock.density'])
  disp('Another process is already working on this directory.')
  disp('If this is not the case, remove the file')
  disp([pwd '/outp/datfiles/lock.density'])
else
  all_is_fine=logical(0);
  try
    dlmwrite('outp/datfiles/lock.density',dummy,'');
    fid=fopen('outp/datfiles/lock.density','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    % The action is put inside this try block. The purpose of this is that if
    % an error happens after the writing of the lock file, this lock file
    % shall be removed so that future attempts are not blocked.
    try
      dd = dir('outp/datfiles/density/n*p00000.picard.dat');
      Ntimesteps = length(dd);
      if Ntimesteps>0
        disp('Working on density.')
      end
        for gg = 1:Ntimesteps
        for hh = 1:Nspecies
          particle(hh).density=zeros(Nx,Ny,Nz);
        end
        iterationstring = dd(gg).name(3:9);
        iteration = str2num(iterationstring);
        dd1 = dir(['outp/datfiles/density/n_' iterationstring 'p*.picard.dat']);
        if length(dd1) ~= Nprocs
          error('Some file or other would seem to be missing.')
        end
        for kk = 0:kprocs-1
          for jj = 0:jprocs-1
            for ii = 0:iprocs-1
              process = ii+iprocs*jj+iprocs*jprocs*kk;
              data = load(['outp/datfiles/density/' dd1(process+1).name]);
              for hh = 1:Nspecies
                nlocal=reshape(data(:,hh),Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
                particle(hh).density(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                                     jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                                     kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                    nlocal(2:end-1,2:end-1,2:end-1);
                clear nlocal
              end
            end
          end
        end
        filename = ['outp/density' iterationstring '.mat'];
        save(filename, '-v7.3', 'particle', 'Nprocs', ...
             'x', 'y', 'z', 'xcorn', 'ycorn', 'zcorn') 
        for ii = 1:Nprocs
          delete(['outp/datfiles/density/' dd1(ii).name]);
        end
      end
      delete('outp/datfiles/lock.density')
      if Ntimesteps>0
        disp('Density: done!')
      end
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('outp/datfiles/lock.density')
    end % end try
  end %  end if all_is_fine
end % end if exist('datfiles/lock.density')


% Efield
% Prevent two processes from performing simultaneous conversions
if exist([pwd '/outp/datfiles/lock.Efield'])
  disp('Another process is already working on this directory.')
  disp('If this is not the case, remove the file')
  disp([pwd '/outp/datfiles/lock.Efield'])
else
  all_is_fine=logical(0);
  try
    dlmwrite('outp/datfiles/lock.Efield',dummy,'');
    fid=fopen('outp/datfiles/lock.Efield','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    % The action is put inside this try block. The purpose of this is that if
    % an error happens after the writing of the lock file, this lock file
    % shall be removed so that future attempts are not blocked.
    try
      dd = dir('outp/datfiles/Efield/E*p00000.picard.dat');
      Ntimesteps = length(dd);
      if Ntimesteps>0
        disp('Working on Efield.')
      end
      for gg = 1:Ntimesteps
        Ex=zeros(Nx,Ny,Nz);
        Ey=zeros(Nx,Ny,Nz);
        Ez=zeros(Nx,Ny,Nz);
        iterationstring = dd(gg).name(3:9);
        iteration = str2num(iterationstring);
        dd1 = dir(['outp/datfiles/Efield/E_' iterationstring 'p*.picard.dat']);
        if length(dd1) ~= Nprocs
          error('Some file or other would seem to be missing.')
        end
        for kk = 0:kprocs-1
          for jj = 0:jprocs-1
            for ii = 0:iprocs-1
              process = ii+iprocs*jj+iprocs*jprocs*kk;
              data = load(['outp/datfiles/Efield/' dd1(process+1).name]);
              Exlocal = reshape(data(:,1),Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
              Eylocal = reshape(data(:,2),Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
              Ezlocal = reshape(data(:,3),Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
              Ex(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                 jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                 kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                  Exlocal(2:end-1,2:end-1,2:end-1);
              Ey(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                 jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                 kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                  Eylocal(2:end-1,2:end-1,2:end-1);
              Ez(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                 jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                 kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                  Ezlocal(2:end-1,2:end-1,2:end-1);
              clear Exlocal Eylocal Ezlocal
            end
    end
        end
        filename = ['outp/Efield' iterationstring '.mat'];
        save(filename, '-v7.3', 'Ex','Ey','Ez', 'Nprocs', ...
             'x', 'y', 'z', 'xcorn', 'ycorn', 'zcorn') 
        for ii = 1:Nprocs
          delete(['outp/datfiles/Efield/' dd1(ii).name]);
        end
      end
      delete('outp/datfiles/lock.Efield')
      if Ntimesteps>0
        disp('Efield: done!')
      end
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('outp/datfiles/lock.Efield')
    end % end try
  end %  end if all_is_fine
end % end if exist('datfiles/lock.Efield')

% Potential
% Prevent two processes from performing simultaneous conversions
if exist([pwd '/outp/datfiles/lock.potential'])
  disp('Another process is already working on this directory.')
  disp('If this is not the case, remove the file')
  disp([pwd '/outp/datfiles/lock.potential'])
else
  all_is_fine=logical(0);
  try
    dlmwrite('outp/datfiles/lock.potential',dummy,'');
    fid=fopen('outp/datfiles/lock.potential','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    % The action is put inside this try block. The purpose of this is that if
    % an error happens after the writing of the lock file, this lock file
    % shall be removed so that future attempts are not blocked.
    try
      dd = dir('outp/datfiles/potential/UE*p00000.picard.dat');
      Ntimesteps = length(dd);
      if Ntimesteps>0
        disp('Working on potential.')
      end
      for gg = 1:Ntimesteps
        U=zeros(Nx,Ny,Nz);
        iterationstring = dd(gg).name(4:10);
        iteration = str2num(iterationstring);
        dd1 = dir(['outp/datfiles/potential/UE_' iterationstring ...
                   'p*.picard.dat']);
        if length(dd1) ~= Nprocs
          error('Some file or other would seem to be missing.')
        end
        for kk = 0:kprocs-1
          for jj = 0:jprocs-1
            for ii = 0:iprocs-1
              process = ii+iprocs*jj+iprocs*jprocs*kk;
              data = load(['outp/datfiles/potential/' dd1(process+1).name]);
              Ulocal = reshape(data,Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
              U(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                  Ulocal(2:end-1,2:end-1,2:end-1);
              clear Ulocal
            end
          end
        end
        filename = ['outp/potential' iterationstring '.mat'];
        save(filename, '-v7.3', 'U', 'Nprocs', ...
             'x', 'y', 'z', 'xcorn', 'ycorn', 'zcorn') 
        for ii = 1:Nprocs
          delete(['outp/datfiles/potential/' dd1(ii).name]);
        end
      end
      delete('outp/datfiles/lock.potential')
      if Ntimesteps>0
        disp('Potential: done!')
      end
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('outp/datfiles/lock.potential')
    end % end try
  end %  end if all_is_fine
end % end if exist('datfiles/lock.potential')

% Flux
% Prevent two processes from performing simultaneous conversions
if exist([pwd '/outp/datfiles/lock.flux'])
  disp('Another process is already working on this directory.')
  disp('If this is not the case, remove the file')
  disp([pwd '/outp/datfiles/lock.flux'])
else
  all_is_fine=logical(0);
  try
    dlmwrite('outp/datfiles/lock.flux',dummy,'');
    fid=fopen('outp/datfiles/lock.flux','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    % The action is put inside this try block. The purpose of this is that if
    % an error happens after the writing of the lock file, this lock file
    % shall be removed so that future attempts are not blocked.
    try
      disp('Working on flux.')
      for hh = 1:Nspecies
        dd = dir(['outp/datfiles/flux/F' num2str(hh,'%2.2i') ...
                  '*p00000.picard.dat']);
        Ntimesteps = length(dd);
        for gg = 1:Ntimesteps
          Fx=zeros(Nx,Ny,Nz);
          Fy=zeros(Nx,Ny,Nz);
          Fz=zeros(Nx,Ny,Nz);
          iterationstring = dd(gg).name(5:11);
          iteration = str2num(iterationstring);
          dd1 = dir(['outp/datfiles/flux/F' num2str(hh,'%2.2i') '_' ...
                     iterationstring 'p*.picard.dat']);
          if length(dd1) ~= Nprocs
            error('Some file or other would seem to be missing.')
          end
          for kk = 0:kprocs-1
            for jj = 0:jprocs-1
              for ii = 0:iprocs-1
                process = ii+iprocs*jj+iprocs*jprocs*kk;
                data = load(['outp/datfiles/flux/' dd1(process+1).name]);
                Fxlocal = ...
                    reshape(data(:,1),Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
                Fylocal = ...
                    reshape(data(:,2),Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
                Fzlocal = ...
                    reshape(data(:,3),Nx/iprocs+2,Ny/jprocs+2,Nz/kprocs+2);
                Fx(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                   jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                   kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                    Fxlocal(2:end-1,2:end-1,2:end-1);
                Fy(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                   jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                   kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                    Fylocal(2:end-1,2:end-1,2:end-1);
                Fz(ii*Nx/iprocs+1:(ii+1)*Nx/iprocs, ...
                   jj*Ny/jprocs+1:(jj+1)*Ny/jprocs, ...
                   kk*Nz/kprocs+1:(kk+1)*Nz/kprocs) = ...
                    Fzlocal(2:end-1,2:end-1,2:end-1);
                clear Fxlocal Fylocal Fzlocal
              end
            end
          end
          filename = ['outp/flux' num2str(hh,'%2.2i') '_' ...
                      iterationstring '.mat'];
          save(filename, '-v7.3', 'Fx','Fy','Fz', 'Nprocs', ...
               'x', 'y', 'z', 'xcorn', 'ycorn', 'zcorn') 
          for ii = 1:Nprocs
            delete(['outp/datfiles/flux/' dd1(ii).name]);
          end
        end
      end
      delete('outp/datfiles/lock.flux')
      disp('Flux: done!')
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('outp/datfiles/lock.flux')
    end % end try
  end %  end if all_is_fine
end % end if exist('datfiles/lock.flux')

% Eprobe
% Prevent two processes from performing simultaneous conversions
if exist([pwd '/outp/datfiles/lock.Eprobe'])
  disp('Another process is already working on this directory.')
  disp('If this is not the case, remove the file')
  disp([pwd '/outp/datfiles/lock.Eprobe'])
else
  all_is_fine=logical(0);
  try
    dlmwrite('outp/datfiles/lock.Eprobe',dummy,'');
    fid=fopen('outp/datfiles/lock.Eprobe','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    % The action is put inside this try block. The purpose of this is that if
    % an error happens after the writing of the lock file, this lock file
    % shall be removed so that future attempts are not blocked.
    try
      disp('Working on Eprobe.')
      for probe = 1:Nprobes
        dd=dir(['outp/datfiles/Eprobe/Eprobe*p' ...
                num2str(probe,'%4.4i') '.picard.dat']);
        if length(dd)>0
          EprobeInExistence=logical(0);
          matfilename=[pwd '/outp/Eprobe' num2str(probe,'%4.4i') '.mat'];
          if exist(matfilename)
            EprobeInExistence=logical(1);
          end
          if EprobeInExistence
            load(matfilename)
            hightimes=[];
            for ii=1:length(dd)
              hightimes(ii)=str2num(dd(ii).name(16:22));
            end
            maxtimeinfile=timestepsEprobe(end);
            for ii=find(hightimes>maxtimeinfile)
              data=load(['outp/datfiles/Eprobe/' dd(ii).name]);
              timesteps=data(:,1); inte = timesteps>maxtimeinfile;
              timestepsEprobe = [timestepsEprobe;timesteps(inte,1)];
              Eprobe = [Eprobe;data(inte,2:4)];
              maxtimeinfile=timestepsEprobe(end);
            end
          else
            Eprobe = []; timestepsEprobe=[];
            for ii=1:length(dd)
              data=load(['outp/datfiles/Eprobe/' dd(ii).name]);
              timestepsEprobe = [timestepsEprobe;data(:,1)];
              Eprobe = [Eprobe;data(:,2:4)];
            end
          end % end if EprobeInExistence
          save(matfilename,'-v7.3','probestruct','timestepsEprobe', ...
               'Eprobe');
          disp(['saved ' matfilename])
          clear Eprobe timestepsEprobe
          for ii=1:length(dd)
            delete(['outp/datfiles/Eprobe/' dd(ii).name])
          end
        end % end if length(dd)>0
      end % end for probe = 1:Nprobes
      delete('outp/datfiles/lock.Eprobe')
      disp('Eprobe: done!')
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('outp/datfiles/lock.Eprobe')
    end % end try
  end % end if all_is_fine
end % end if exist([pwd '/outp/datfiles/lock.Eprobe'])


% Pprobe
% Prevent two processes from performing simultaneous conversions
if exist([pwd '/outp/datfiles/lock.Pprobe'])
  disp('Another process is already working on this directory.')
  disp('If this is not the case, remove the file')
  disp([pwd '/outp/datfiles/lock.Pprobe'])
else
  all_is_fine=logical(0);
  try
    dlmwrite('outp/datfiles/lock.Pprobe',dummy,'');
    fid=fopen('outp/datfiles/lock.Pprobe','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    % The action is put inside this try block. The purpose of this is that if
    % an error happens after the writing of the lock file, this lock file
    % shall be removed so that future attempts are not blocked.
    try
      disp('Working on Pprobe.')
      for probe = 1:Nprobes
        dd=dir(['outp/datfiles/Pprobe/Pprobe*p' ...
                num2str(probe,'%4.4i') '.picard.dat']);

        for ii=1:length(dd)
          matfilename=[pwd '/outp/' dd(ii).name(1:end-11) '.mat'];
          data=load(['outp/datfiles/Pprobe/' dd(ii).name]);
          Pprobe=struct([]);
          for sp=1:Nspecies; Pprobe(sp).r=[]; Pprobe(sp).v=[]; end
          ss=size(data);
          for jj=1:ss(1)
            Pprobe(data(jj,7)).r = [Pprobe(data(jj,7)).r; data(jj,1:3)];
            Pprobe(data(jj,7)).v = [Pprobe(data(jj,7)).v; data(jj,4:6)];
          end
          save(matfilename,'-v7.3','probestruct','Pprobe');
          disp(['saved ' matfilename])
          clear Pprobe data
        end

        for ii=1:length(dd)
          delete(['outp/datfiles/Pprobe/' dd(ii).name])
        end
      end % end for probe = 1:Nprobes
      delete('outp/datfiles/lock.Pprobe')
      disp('Pprobe: done!')
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('outp/datfiles/lock.Pprobe')
    end % end try
  end % end if all_is_fine
end % end if exist([pwd '/outp/datfiles/lock.Pprobe'])

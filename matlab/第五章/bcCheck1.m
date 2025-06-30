%--------------------------------------------------------------------------
function [kk2,mm2,ff2,sdof2]=bcCheck1(kk,mm,ff,bcdof,bcval)
%--------------------------------------------------------------------------
%  Purpose:
%     Check the boundary conditions and eliminate the columns and rows
%     associated with the boundary conditions in the system equation
%  Synopsis:
%     [kk2,mm2,ff2,sdof2]=bcChec1k(kk,mm,ff,,bcdof,bcval)
%  Variable Description:
%     kk, mm = system stiffness matrix and system mass matrix
%     ff = system force vector
%     bcdof = a vector containing dofs associated with boundary conditions
%     bcval = a vector containing boundary condition values associated with the dofs in 'bcdof'
%--------------------------------------------------------------------------
%  (1) check the boundary conditions
%--------------------------------------------------------------------------
[sdof,n1]=size(kk);

jk=0;

for ii=1:sdof                             % loop for check of the boundary conditions

  check=bcdof(ii);
  if check==1
    jk=jk+1;                                   % location of the zero main element
    bci(jk)=ii;                         % storing the location of the zero main element
  end

end
%--------------------------------------------------------------------------
%  (2) eliminate the columns and rows associated with bcs
%--------------------------------------------------------------------------
if jk~=0

  for jj=jk:-1:1                              % loop for moving the columns and rows
                                            % associated with boundary conditions
    hh=sdof-bci(jj);
    hr=bci(jj);

    for i=1:hh                              % exchanging the columns of the kk, etc.
      kt=kk(hr+i-1,:);
      mt=mm(hr+i-1,:);
      ft=ff(hr+i-1);
      kk(hr+i-1,:)=kk(hr+i,:);
      mm(hr+i-1,:)=mm(hr+i,:);
      ff(hr+i-1)=ff(hr+i);
      kk(hr+i,:)=kt;
      mm(hr+i,:)=mt;
      ff(hr+i)=ft;
    end

    for j=1:hh                              % exchanging the rows of the kk, mm, ff
      kt=kk(:,hr+j-1);
      mt=mm(:,hr+j-1);
      kk(:,hr+j-1)=kk(:,hr+j);
      mm(:,hr+j-1)=mm(:,hr+j);
      kk(:,hr+j)=kt;
      mm(:,hr+j)=mt;
    end

  end

end

sdof2=sdof-jk;
kk2=[kk(1:sdof2,1:sdof2)];                 % eliminating the columns and rows of the kk
mm2=[mm(1:sdof2,1:sdof2)];             % eliminating the columns and rows of the mm
ff2=[ff(1:sdof2)];                                    % eliminating the rows of the ff
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------

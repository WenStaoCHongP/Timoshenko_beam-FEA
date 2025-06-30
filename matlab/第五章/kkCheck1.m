%--------------------------------------------------------------------------
function [kk1,mm1,ff1,bcdof1,bcval1,sdof1]=kkCheck1(kk,mm,ff,bcdof,bcval)
%--------------------------------------------------------------------------
%  Purpose:
%     Check whether the main elements of system stiffness matrix are zeros and eliminate
%     the columns and rows associated with the zero main element in the system matrix
%     equation.
%  Synopsis:
%     [kk1,mm1,ff1,bcdof1,bcval1]=kkcheck1(kk,mm,ff,bcdof,bcval)
%  Variable Description:
%     kk = system stiffness matrix
%     mm = system mass vector
%     ff = system force vector
%     sdof = total degrees of freedom of the system
%     bcdof = a vector containing dofs associated with boundary conditions
%     bcval = a vector containing boundary condition values associated with the dofs in 'bcdof'
%--------------------------------------------------------------------------
%  (1) check the zero main element in the kk
%--------------------------------------------------------------------------
[sdof,n1]=size(kk);

jk=0;

for ii=1:sdof                         % loop for check of the zero main elements in kk
  check=kk(ii,ii);
  if check==0
    jk=jk+1;                                   % location of the zero main element
    kki(jk)=ii;                        % storring the location of the zero main element
  end
end

if jk~=0
  disp('main elements of kk:')                         % display the zero main element
  kkii=kki
  disp('equal zero')
end
%--------------------------------------------------------------------------
%  (2) eliminate the columns and rows in kk, mm, ff, bcdof, bcval
%--------------------------------------------------------------------------
if jk~=0
  for jj=jk:-1:1                 % loop for moving the columns and rows associated with
                               % the zero main element in the system matrix equation
    hh=sdof-kki(jj);
    hr=kki(jj);

    for i=1:hh                                % exchanging the rows in the equation
      kt=kk(hr+i-1,:);
      mt=mm(hr+i-1,:);
      ft=ff(hr+i-1);
      bdt=bcdof(hr+i-1);
      bvt=bcval(hr+i-1,:);

      kk(hr+i-1,:)=kk(hr+i,:);
      mm(hr+i-1,:)=mm(hr+i,:);
      ff(hr+i-1)=ff(hr+i);
      bcdof(hr+i-1)=bcdof(hr+i);
      bcval(hr+i-1,:)=bcval(hr+i,:);

      kk(hr+i,:)=kt;
      mm(hr+i,:)=mt;
      ff(hr+i)=ft;
      bcdof(hr+i,:)=bdt;
      bcval(hr+i,:)=bvt;
    end

    for j=1:hh                             % exchanging the columns in the equation
      kt=kk(:,hr+j-1);
      mt=mm(:,hr+j-1);
      kk(:,hr+j-1)=kk(:,hr+j);
      mm(:,hr+j-1)=mm(:,hr+j);
      kk(:,hr+j)=kt;
      mm(:,hr+j)=mt;
    end

  end

end

sdof1=sdof-jk;
kk1=[kk(1:sdof1,1:sdof1)];                % eliminating the columns and rows of the kk
mm1=[mm(1:sdof1,1:sdof1)];             % eliminating the columns and rows of the mm
ff1=[ff(1:sdof1)];                                    % eliminating the rows of the ff
bcdof1=[bcdof(1:sdof1)];                          % eliminating the rows of the bcdof
bcval1=[bcval(1:sdof1,:)];                          % eliminating the rows of the bcval
%--------------------------------------------------------------------------
%     The end
%-------------------------------------------------------------------------- 

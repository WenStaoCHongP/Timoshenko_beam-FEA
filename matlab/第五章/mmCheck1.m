%--------------------------------------------------------------------------
function [kk2,mm2,ff2,bcdof2,bcval2,sdof2]=mmCheck1(kk,mm,ff,bcdof,bcval)
%--------------------------------------------------------------------------
%  Purpose:
%     Check whether the main elements of system mass matrix are zeros and eliminate
%     the columns and rows associated with the zero main element in the system matrices
%     kk, mm, ff, bcdof, bcval
%  Synopsis:
%     [kk2,mm2,ff2,bcdof2,bcval2,sdof2]=mmCheck1(kk,mm,ff,bcdof,bcval)
%  Variable Description:
%     kk, mm = system stiffness matrix and system mass matrix
%     ff = system force vector
%     bcdof = a vector containing dofs associated with boundary conditions
%     bcval = a vector containing boundary condition values associated with the dofs in 'bcdof'
%     sdof = total degrees of freedom of the system
%--------------------------------------------------------------------------
%  (1) check the zero main element in the mm
%--------------------------------------------------------------------------
[sdof,n2]=size(kk);

jk=0;

for ii=1:sdof                         % loop for check of the zero main elements in mm
  check=mm(ii,ii);
  if check==0
    jk=jk+1;                                   % location of the zero main element
    mmi(jk)=ii;                        % storing the location of the zero main element
  end
end
%--------------------------------------------------------------------------
%  (2) eliminate the columns and rows in kk, mm, ff, bcdof, bcval
%--------------------------------------------------------------------------
if jk~=0
  for jj=jk:-1:1                 % loop for moving the columns and rows associated with
                               % the zero main element in the system matrix equation
    hh=sdof-mmi(jj);
    hr=mmi(jj);

    for i=1:hh                                % exchanging the rows in the equation
      kt=kk(hr+i-1,:);
      mt=mm(hr+i-1,:);
      ft=ff(hr+i-1);
      bdt=bcdof(hr+i-1);
      bvt=bcval(hr+i-1);
      kk(hr+i-1,:)=kk(hr+i,:);
      mm(hr+i-1,:)=mm(hr+i,:);
      ff(hr+i-1)=ff(hr+i);
      bcdof(hr+i-1)=bcdof(hr+i);
      bcval(hr+i-1)=bcval(hr+i);
      kk(hr+i,:)=kt;
      mm(hr+i,:)=mt;
      ff(hr+i)=ft;
      bcdof(hr+i)=bdt;
      bcval(hr+i)=bvt;
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
 
sdof2=sdof-jk;
kk2=[kk(1:sdof2,1:sdof2)];                 % eliminating the columns and rows of the kk
mm2=[mm(1:sdof2,1:sdof2)];             % eliminating the columns and rows of the mm
ff2=[ff(1:sdof2)];                                    % eliminating the rows of the ff
bcdof2=[bcdof(1:sdof2)];                          % eliminating the rows of the bcdof
bcval2=[bcval(1:sdof2)];                           % eliminating the rows of the bcval
%--------------------------------------------------------------------------
%     The end
%--------------------------------------------------------------------------

% usage: 
%   1. call as a regular function in Matlab
%   2. matlab -r "extract_gridline(<initial file number>,
%                                  <final file number>,
%                                  <file number between two consecutive files>,
%                                  <block ID at which grid line is extracted>,
%                                  <tuple describing how grid line is extracted>,
%                                  '<coordinate name along with grid line is extracted>',
%                                  '<which Plot3D format>;)"
%      e.g. matlab -r "preproc(2,100,4,3,[0 5 1],'x','solution)"
%             Extract grid lines starting from XXXXXXX.000002.q to XXXXXXX.000100.q at every 4 file; 
%             Plot3D solution format is assumed in those q files; 
%             Grid lines are extracted along XI direction at ETA = 5 and ZETA = 1 at the block 3; 
%             XI direction is aligned with the x axis
% for batch jobs, usage 2 is recommended; make sure to encircle the function name by " " and to use ' ' for string arguments

% Written by Jeonglae Kim, August 2017

function extract_gridline(count_begin,count_end,count_jump, ...
                          which_block,which_gridline,which_xyz,which_plot3d)

% be more tidy
%clc; clear all; close all;
format compact;

% constants
XDIR = 1; YDIR = 2; ZDIR = 3;
          RDIR = 2;
XI = 1; ETA = 2; ZETA = 3;
TRUE = 1; FALSE = 0;

% read grid
[numBlocks,numPoints,xyz] = read_grid_plot3d('../camille.xyz');
assert(which_block <= numBlocks);
dir_extract_not = find(which_gridline); % e.g. if which_gridline = [10 0 3], dir_extract_not = [1 3]
dir_extract = find(~which_gridline); % e.g. if which_gridline = [10 0 3],  dir_extract = [2]
if (length(dir_extract_not) ~= 2)
  fprintf('Only two elements of which_gridline can be non-zero, so that line data can be extracted along the third one; e.g. [31 0 1], [0 2 3].\n');
  exit
end % length(dir_extract_not)
for idir = 1:2 % check whether two non-zero grid-point indices are in range in their directions, respectively
  assert(which_gridline(dir_extract_not(idir)) <= numPoints{which_block}(dir_extract_not(idir)));
end % idir
if (strcmp(which_xyz,'x') == TRUE || strcmp(which_xyz,'X') == TRUE)
  tmp = XDIR;
elseif (strcmp(which_xyz,'y') == TRUE || strcmp(which_xyz,'Y') == TRUE || strcmp(which_xyz,'r') == TRUE || strcmp(which_xyz,'R') == TRUE)
  tmp = YDIR;
elseif (strcmp(which_xyz,'z') == TRUE || strcmp(which_xyz,'Z') == TRUE)
  tmp = ZDIR;
end % strcmp(which_xyz,'x')
which_xyz = tmp;

fid_out = fopen('contour_xt.dat','wt');
for count = count_begin:count_jump:count_end

  % read solution data
  fname_solution = sprintf('../camille.%6.6d.q',count);
  fprintf('Working on %s\n',fname_solution);
  silence = TRUE;
  if (which_plot3d == 'solution')
    [numBlocks_solution,numPoints_solution,timestep,time,q] = read_solution_plot3d(fname_solution,silence);
    numVars = 5;

  elseif (which_plot3d == 'function')
    [numBlocks_solution,numPoints_solution,numVars,q] = read_function_plot3d(fname_solution,silence);

  end % which_plot3d

  % file header
  if (count == count_begin)
    fprintf(fid_out,'variables = loc,count');
    %
    % name file has variable names
    varname = cell(numVars,1);
    fid_in = fopen('../camille.nam','rt');
    for ivar = 1:numVars
      varname{ivar} = fgetl(fid_in);
    end % ivar
    fclose(fid_in);
    %
    for ivar = 1:numVars
      fprintf(fid_out,',%s',varname{ivar});
    end % ivar
    fprintf(fid_out,'\n');
    fprintf(fid_out,'zone i = %d,j = %d,f = point\n',numPoints{which_block}(dir_extract),floor((count_end-count_begin)/count_jump)+1);
  end % count

  % consistency check between grid and solution
  assert(numBlocks == numBlocks_solution);
  for ib = 1:numBlocks
    for idir = XI:ZETA
      assert(numPoints{ib}(idir) == numPoints_solution{ib}(idir));
    end % idir
  end % ib

  for index = 1:numPoints{which_block}(dir_extract)
    if (dir_extract == XI)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(index,which_gridline(ETA),which_gridline(ZETA)),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(index,which_gridline(ETA),which_gridline(ZETA)));
        end % ivar
      fprintf(fid_out,'\n');
    elseif (dir_extract == ETA)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(which_gridline(XI),index,which_gridline(ZETA)),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(which_gridline(XI),index,which_gridline(ZETA)));
        end % ivar
      fprintf(fid_out,'\n');
    elseif (dir_extract == ZETA)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(which_gridline(XI),which_gridline(ETA),index),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(which_gridline(XI),which_gridline(ETA),index));
        end % ivar
      fprintf(fid_out,'\n');
    end % dir_extract
  end % index

  clear numBlocks_solution numPoints_solution q;

end % count
fclose(fid_out);

% problem specifics
[numBlocks_function,numPoints_function,numVars,q] = read_function_plot3d('../mean_camille.q',silence);
%
fid_out = fopen('contour_xt_mean.dat','wt');
fprintf(fid_out,'variables = loc,count');
%
% name file has variable names
clear varname;
varname = cell(numVars,1);
fid_in = fopen('../mean_camille.nam','rt');
for ivar = 1:numVars
  varname{ivar} = fgetl(fid_in);
end % ivar
fclose(fid_in);
%
for ivar = 1:numVars
  fprintf(fid_out,',%s',varname{ivar});
end % ivar
fprintf(fid_out,'\n');
fprintf(fid_out,'zone i = %d,j = %d,f = point\n',numPoints{which_block}(dir_extract),floor((count_end-count_begin)/count_jump)+1);
%
for count = count_begin:count_jump:count_end
  for index = 1:numPoints{which_block}(dir_extract)
    if (dir_extract == XI)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(index,which_gridline(ETA),which_gridline(ZETA)),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(index,which_gridline(ETA),which_gridline(ZETA)));
        end % ivar
      fprintf(fid_out,'\n');
    elseif (dir_extract == ETA)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(which_gridline(XI),index,which_gridline(ZETA)),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(which_gridline(XI),index,which_gridline(ZETA)));
        end % ivar
      fprintf(fid_out,'\n');
    elseif (dir_extract == ZETA)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(which_gridline(XI),which_gridline(ETA),index),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(which_gridline(XI),which_gridline(ETA),index));
        end % ivar
      fprintf(fid_out,'\n');
    end % dir_extract
  end % index
end % count
clear q;

[numBlocks_function,numPoints_function,numVars,q] = read_function_plot3d('../aux_camille.000000.q',silence);
%
fid_out = fopen('contour_xt_aux.dat','wt');
fprintf(fid_out,'variables = loc,count');
%
% name file has variable names
clear varname;
numVars = 7;
varname = cell(numVars,1);
fid_in = fopen('../aux_camille.nam','rt');
for ivar = 1:numVars
  varname{ivar} = fgetl(fid_in);
end % ivar
fclose(fid_in);
%
for ivar = 1:numVars
  fprintf(fid_out,',%s',varname{ivar});
end % ivar
fprintf(fid_out,'\n');
fprintf(fid_out,'zone i = %d,j = %d,f = point\n',numPoints{which_block}(dir_extract),floor((count_end-count_begin)/count_jump)+1);
%
for count = count_begin:count_jump:count_end
  for index = 1:numPoints{which_block}(dir_extract)
    if (dir_extract == XI)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(index,which_gridline(ETA),which_gridline(ZETA)),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(index,which_gridline(ETA),which_gridline(ZETA)));
        end % ivar
      fprintf(fid_out,'\n');
    elseif (dir_extract == ETA)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(which_gridline(XI),index,which_gridline(ZETA)),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(which_gridline(XI),index,which_gridline(ZETA)));
        end % ivar
      fprintf(fid_out,'\n');
    elseif (dir_extract == ZETA)
        fprintf(fid_out,'%15.6e %15d',xyz{which_xyz,which_block}(which_gridline(XI),which_gridline(ETA),index),count);
        for ivar = 1:numVars
          fprintf(fid_out,' %15.6e',q{ivar,which_block}(which_gridline(XI),which_gridline(ETA),index));
        end % ivar
      fprintf(fid_out,'\n');
    end % dir_extract
  end % index
end % count
clear q;

exit

end

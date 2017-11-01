% This script reads a Plot3D grid file. 
% To be consistent with the flow solver, the grid file assumes a xyz file with 
%   MULTIPLE grids
%   WHOLE format
%   IBLANK
% See doc/plot3d/plot3d_manual_ch8.pdf to see what they indicate.

% Written by Jeonglae Kim, August 2017

function [numBlocks,numPoints,xyz] = read_grid_plot3d(fname_grid,silence)

% be more tidy
%clc; clear all; close all;
format compact;

% constants
XDIR = 1; YDIR = 2; ZDIR = 3;
XI = 1; ETA = 2; ZETA = 3;
TRUE = 1; FALSE = 0;
numVars = 3;

if (nargin == 1)
  silence = FALSE;
end % nargin

% read grid
fid_in = fopen(fname_grid,'rb');
if (silence == FALSE)
  fprintf('Reading a grid file named %s\n',fname_grid);
end % silence
%
numBlocks = fread(fid_in,1,'int'); % number of blocks
if (silence == FALSE)
  fprintf('You have %d block(s) in your grid file.\n\n',numBlocks);
end % silence
%
numPoints = cell(numBlocks,1);
if (silence == FALSE)
  fprintf('Block | # of points in XI | # of points in ETA | # of points in ZETA\n');
end % silence
for ib = 1:numBlocks
  numPoints{ib} = fread(fid_in,3,'int'); % number of points in \xi, \eta, & \zeta directions
  if (silence == FALSE)
    fprintf('%5d %19d %20d %21d\n',ib,numPoints{ib}(XI),numPoints{ib}(ETA),numPoints{ib}(ZETA));
  end % silence
end % ib
%
xyz = cell(numVars,numBlocks);
for ib = 1:numBlocks
  numPoints_thisBlock = prod(numPoints{ib}(XI:ZETA));
  for idir = XDIR:ZDIR
    xyz{idir,ib} = fread(fid_in,numPoints_thisBlock,'double');
    xyz{idir,ib} = reshape(xyz{idir,ib},[numPoints{ib}(XI) numPoints{ib}(ETA) numPoints{ib}(ZETA)]);
  end % idir
end % ib
if (silence == FALSE)
  fprintf('\nxyz point information has been all read and reshaped.\n\n');
end % silence

fclose(fid_in);

end

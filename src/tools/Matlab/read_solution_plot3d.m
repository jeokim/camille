% This script reads a Plot3D solution file. 
% To be consistent with the flow solver, the solution file assumes a q file with 
%   MULTIPLE grids
%   WHOLE format
%   without Jacobian
% See doc/plot3d/plot3d_manual_ch8.pdf to see what they indicate.

% Written by Jeonglae Kim, August 2017

function [numBlocks,numPoints,timestep,time,q] = read_solution_plot3d(fname_solution,silence)

% be more tidy
%clc; clear all; close all;
format compact;

% constants
XDIR = 1; YDIR = 2; ZDIR = 3;
XI = 1; ETA = 2; ZETA = 3;
TRUE = 1; FALSE = 0;
numVars = 5;

if (nargin == 1)
  silence = FALSE;
end % nargin

% read solution
fid_in = fopen(fname_solution,'rb');
if (silence == FALSE)
  fprintf('Reading a solution file named %s\n',fname_solution);
end % silence
%
numBlocks = fread(fid_in,1,'int'); % number of blocks
if (silence == FALSE)
  fprintf('You have %d block(s) in your solution file.\n\n',numBlocks);
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
q = cell(numVars,numBlocks);
for ib = 1:numBlocks
  tau = fread(fid_in,4,'double');
  timestep = tau(3);
  time = tau(4);

  numPoints_thisBlock = prod(numPoints{ib}(XI:ZETA));
  for ivar = 1:numVars
    q{ivar,ib} = fread(fid_in,numPoints_thisBlock,'double');
    q{ivar,ib} = reshape(q{ivar,ib},[numPoints{ib}(XI) numPoints{ib}(ETA) numPoints{ib}(ZETA)]);
  end % ivar
end % ib
if (silence == FALSE)
  fprintf('\nSolution fields have been all read and reshaped.\n\n');
end % silence

fclose(fid_in);

end
clc; clear;

LASTN = maxNumCompThreads(1);

% some assumptions all over the code:
% -- tol for solves is always 1.0e-8
% -- the W-cycle and smoothing iters are set in two_level_mg.m
% -- when singular vectors are mentioned, it's of D and not D^{-1}

% -- when computing 4D traces, we work with Ptilde only i.e. there's no
%    extra 'front factor' B
% -- when computing 3D traces, Ptilde=identity but B!=identity

%%
% CHANGEABLE PARAMS

% the shift on the lattice .. as we're doing here a lattice with a
% displacement. The displacement is on the second direction, that's why
% the 16 extra factor
nr_displ_sites = 0;
% number of iterations within Block Power Iteration
bpi_iters = 0;
% CASE=1 is Hutchinson, CASE=2 is MGMLMC
CASE = 1;
% for CASE=2, choose the level to compute the variance of
% not needed anymore
level_nr = 1;
% number of deflation vectors
k = 16;
% size of the sample to use to estimate the variance
sample_size = 100;
%if CASE==2 && level_nr>1 && k>0
% disabling deflated MGMLMC completely
% if CASE==2 && k>0
%   error("Deflated MGMLMC has been disabled completely for now\n");
% end
% set global param indicating whether we're dealing with
% 3D traces (within FOR5269) or not
global do_3D_traces;
do_3D_traces = 1;

% this is meant for 3D traces only!
use_W_identity = 1;

if do_3D_traces==1 && nr_displ_sites~=0
  error("We have disabled displacements in the lattice when computing ..." + ...
      "3D traces, for now, but this can be easily re-introduced\n")
end

% IMPORTANT : this is a set-able parameter that needs to be included soon,
%             which at the moment we assume as 1 i.e. P3D corresponds to
%             t=1
%t = 1;

%%
% DO NOT CHANGE

addpath("src/");
filename = "matrices/D16/Dreord";
D = get_matrix(filename);

% multigrid hierarchy
nr_levels = 3;
mgh = mg_setup(D,nr_levels,nr_displ_sites,use_W_identity);

% options for defl_type : "EVs", "RSVs", "LSVs"
% do not change this parameter, we're assuming always using RSVs
defl_type = "RSVs";

%% -----------------------------------------------

if CASE==1
  % use Hutchinson
  % options for algorithm : "Hutch", "mgmlmc"
  alg_type = "Hutch";

  fprintf("Case #1 : Hutchinson\n");

  fprintf("\n");
  fprintf("k = %d\n",k);
  % compute the vectors used in inexact deflation
  mgh = compute_deflation_vectors(defl_type,k,mgh,alg_type,bpi_iters);
  % compute the trace
  [tracex,variance,~] = compute_trace(k,mgh,alg_type,1.0e-2,sample_size);
  fprintf("Trace = %f+i%f\n",real(tracex),imag(tracex));
  fprintf("Variance = %f\n",variance);

  fprintf("\n");

  % when running with -nodisplay -nosplash, re-enable this exit
  %exit;
end

%% -----------------------------------------------

if CASE==2
  % use MGMLMC
  % options for algorithm : "Hutch", "mgmlmc"
  alg_type = "mgmlmc";

  fprintf("Case #2 : MGMLMC\n");

  fprintf("\n");
  fprintf("k = %d\n",k);
  % compute the vectors used in inexact deflation
  mgh = compute_deflation_vectors(defl_type,k,mgh,alg_type,bpi_iters);
  % compute the trace
  total_trace = 0.0;
  for level_nr=1:length(mgh.D)-1
    [tracex,variance,~] = compute_trace(k,mgh,alg_type,1.0e-2,sample_size,level_nr);
    fprintf("Trace = %f+i%f\n",real(tracex),imag(tracex));
    fprintf("Variance = %f\n",variance);
    total_trace = total_trace + tracex;
  end
  fprintf("Total trace = %f+i%f\n",real(total_trace),imag(total_trace));
  fprintf("\n");

  % when running with -nodisplay -nosplash, re-enable this exit
  %exit;
end

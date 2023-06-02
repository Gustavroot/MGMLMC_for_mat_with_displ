clc; clear;

LASTN = maxNumCompThreads(1);

% some assumptions all over the code:
% -- tol for solves is always 1.0e-8
% -- the W-cycle and smoothing iters are set in two_level_mg.m
% -- when singular vectors are mentioned, it's of D and not D^{-1}

addpath("src/");

filename = "matrices/D16/Dreord";
D = get_matrix(filename);

% the shift on the lattice .. as we're doing here a lattice with a
% displacement. The displacement is on the second direction, that's why
% the 16 extra factor
nr_displ_sites = 5;

% multigrid hierarchy
nr_levels = 3;
mgh = mg_setup(D,nr_levels,nr_displ_sites);

% options for defl_type : "EVs", "RSVs", "LSVs"
% do not change this parameter, we're assuming always using RSVs
defl_type = "RSVs";
% number of iterations within Block Power Iteration
bpi_iters = 5;

% four cases, as seen in the following if statements
CASE = 1;

%% -----------------------------------------------

if CASE==1
  % options for trace_type : "invD", "invD*perm"
  trace_type = "invD";

  % use Hutchinson
  % options for algorithm : "Hutch", "mgmlmc"
  alg_type = "Hutch";

  fprintf("Set of numerical experiments #1 : Hutchinson, lattice without displacement\n");

  k = 10;
  fprintf("\n");
  fprintf("k = %d\n",k);
  % compute the vectors used in inexact deflation
  mgh = compute_deflation_vectors(defl_type,k,mgh,alg_type,bpi_iters);
  % compute the trace
  [tracex,variance,~] = compute_trace(trace_type,k,mgh,alg_type,1.0e-2,1000);
  fprintf("Trace = %f+i%f\n",real(tracex),imag(tracex));
  fprintf("Variance = %f\n",variance);

  fprintf("\n");

  exit;
end

%% -----------------------------------------------

if CASE==2
  % options for trace_type : "invD", "invD*perm"
  trace_type = "invD";

  % use MGMLMC
  % options for algorithm : "Hutch", "mgmlmc"
  alg_type = "mgmlmc";

  fprintf("Set of numerical experiments #2 : MGMLMC, lattice without displacement\n");

  %error("FIXME : MGMLMC is only computings its first level difference, and GPM is not taken into account\n");

  k = 5; % 28;
  fprintf("\n");
  fprintf("k = %d\n",k);
  % compute the vectors used in inexact deflation
  mgh = compute_deflation_vectors(defl_type,k,mgh,alg_type,bpi_iters);
  % compute the trace
  [tracex,variance,~] = compute_trace(trace_type,k,mgh,alg_type,1.0e-2,1000);
  fprintf("Trace = %f+i%f\n",real(tracex),imag(tracex));
  fprintf("Variance = %f\n",variance);

  fprintf("\n");

  exit;
end

%% -----------------------------------------------

% % options for trace_type : "invD", "invD*perm"
% trace_type = "invD*perm";
% 
% % use Hutchinson
% 
% % options for algorithm : "Hutch", "mgmlmc"
% alg_type = "Hutch";
% 
% fprintf("Set of numerical experiments #3 : Hutchinson, lattice with displacement\n");
% 
% for k=0:50:200
%   % compute the vectors used in inexact deflation
%   X = compute_deflation_vectors(D,defl_type,k,mgh,alg_type);
% 
%   % compute the trace
%   compute_trace(D,trace_type,k,X,mgh,alg_type);
% end
% 
% fprintf("\n");
% 
% % use MGMLMC
% 
% % options for algorithm : "Hutch", "mgmlmc"
% alg_type = "mgmlmc";
% 
% fprintf("Set of numerical experiments #4 : MGMLMC, lattice with displacement\n");
% 
% for k=0:50:200
%   % compute the vectors used in inexact deflation
%   X = compute_deflation_vectors(D,defl_type,k,mgh,alg_type);
% 
%   % compute the trace
%   compute_trace(D,trace_type,k,X,mgh,alg_type);
% end
% 
% fprintf("\n");

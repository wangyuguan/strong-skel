function stats = test_empty_proxy()
%TEST_EMPTY_PROXY Check exact far-field sampling and asymmetric factor use.
% Add FLAM to the MATLAB path, then run this test from the test directory.
  root = fileparts(fileparts(mfilename('fullpath')));
  addpath(root);
  assert(exist('id','file') ~= 0,'Add FLAM to the MATLAB path first.');
  assert(strcmp(which('srskelf_sv_nn'),fullfile(root,'srskelf_sv_nn.m')));
  assert(strcmp(which('srskelf_mv_nn'),fullfile(root,'srskelf_mv_nn.m')));

  [xx,yy] = ndgrid(linspace(0,1,20));
  x = [xx(:)'; yy(:)'];
  N = size(x,2);
  a = x(1,:)'; b = x(2,:)';
  U = [ones(N,1), a+1i*b, exp(a-0.4*b)];
  V = [exp(-b), b-0.3i*a, cos(2*a+b)];
  % Different left/right factors give a complex, nonsymmetric matrix.
  % Local coupling also exercises elimination with uncompressed neighbors.
  A = 2*eye(N) + U*V.'/N + diag(0.03*ones(N-1,1),1);
  B = [sin(3*a)+1i*cos(2*b), exp(1i*(a+2*b))];
  occ = 16;
  t = shypoct(x,occ,Inf,[]);
  assert(t.nlvl >= 4,'Test requires distant boxes beyond the interaction list.');
  first = t.lvp(t.nlvl)+1;
  slf = sort(t.nodes(first).xi);
  nbr = [t.nodes(t.nodes(first).nbor).xi];
  expected_far = setdiff(1:N,[slf nbr]);
  interaction_list = [t.nodes(t.nodes(first).ilist).xi];
  assert(~isempty(setdiff(expected_far,interaction_list)));
  checked_far = false;

  rng_before = rng;
  F = srskelf_asym(@block,x,occ,1e-10,[],struct('symm','n','verb',1));
  assert(checked_far,'No far-field block was checked.');
  assert(isequal(rng_before,rng),'Factorization changed the random state.');
  assert(F.nlvl >= 4);
  assert(F.lvp(2) > 0,'Test must eliminate DOFs at the deepest level.');
  X = srskelf_sv_nn(F,B);
  solve_error = norm(A*X-B,'fro')/norm(B,'fro');
  multiply_error = norm(srskelf_mv_nn(F,B)-A*B,'fro')/norm(A*B,'fro');
  assert(solve_error < 1e-8,'Solve residual too large.');
  assert(multiply_error < 1e-8,'Multiply error too large.');
  stats = struct('N',N,'levels',F.nlvl,'deepest_factors',F.lvp(2), ...
      'solve_residual',solve_error,'multiply_error',multiply_error);
  fprintf('PASS: N=%d, levels=%d, deepest factors=%d\n',N,F.nlvl,F.lvp(2));
  fprintf('solve residual %.3e; multiply error %.3e\n',solve_error,multiply_error);

  function K = block(I,J)
    if ~checked_far
      assert(isequal(I(:),expected_far(:)), ...
          'Empty proxy must sample every active non-neighbor DOF.');
      assert(isequal(J(:),slf(:)));
      checked_far = true;
    end
    K = A(I,J);
  end
end

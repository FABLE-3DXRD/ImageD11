% Dump MTEX direction colour keys (sector, centre, references) in lattice terms.
% Every vector v is written as uvw (v = u a + v b + w c) and hkl (h = v.a, k = v.b, l = v.c).
addpath('~/mtex/mtex-latest'); startup_mtex;
tri = {[5 6 7],[80 95 105]}; mono = {[5 6 7]};
G = {'1',tri; '-1',tri; ...
 '211',{[5 6 7],[100 90 90]}; 'm11',{[5 6 7],[100 90 90]}; '2/m11',{[5 6 7],[100 90 90]}; ...
 '121',{[5 6 7],[90 100 90]}; '1m1',{[5 6 7],[90 100 90]}; '12/m1',{[5 6 7],[90 100 90]}; ...
 '112',{[5 6 7],[90 90 100]}; '11m',{[5 6 7],[90 90 100]}; '112/m',{[5 6 7],[90 90 100]}; ...
 '222',{[5 6 7],[90 90 90]}; '2mm',{[5 6 7],[90 90 90]}; 'm2m',{[5 6 7],[90 90 90]}; 'mm2',{[5 6 7],[90 90 90]}; 'mmm',{[5 6 7],[90 90 90]}};
for g = {'4','-4','4/m','422','4mm','-42m','-4m2','4/mmm'}, G(end+1,:) = {g{1},{[5 5 9],[90 90 90]}}; end
for g = {'3','-3','321','312','3m1','31m','-3m1','-31m','6','-6','6/m','622','6mm','-6m2','-62m','6/mmm'}, G(end+1,:) = {g{1},{[5 5 9],[90 90 120]}}; end
for g = {'23','m-3','432','-43m','m-3m'}, G(end+1,:) = {g{1},{[5 5 5],[90 90 90]}}; end
out = {};
for i = 1:size(G,1)
  nm = G{i,1}; cell = G{i,2};
  cs0 = crystalSymmetry(nm, cell{1}, cell{2}*degree);
  A = [cs0.aAxis.xyz; cs0.bAxis.xyz; cs0.cAxis.xyz];   % rows a,b,c (Cartesian, MTEX frame)
  R = inv(A).';                                       % rows a*,b*,c*
  axs = {'a',A(1,:);'b',A(2,:);'c',A(3,:)}; rec = {'a*',R(1,:);'b*',R(2,:);'c*',R(3,:)};
  % plotting conventions: out-of-screen = +-reciprocal axis, east = +-real axis normal to it, plus the default
  pcs = {'default', cs0.how2plot};
  for k = 1:3, for j = setdiff(1:3,k), for so = [1 -1], for se = [1 -1]
    pcs(end+1,:) = {sprintf('%s%s|%s%s', sgn(so), rec{k,1}, sgn(se), axs{j,1}), ...
      plottingConvention(vector3d(so*R(k,:)), vector3d(se*A(j,:)))};
  end, end, end, end
  % symmetry operations as integer matrices acting on uvw columns: uvw' = M uvw
  ro = cs0.rot; ops = zeros(numel(ro),3,3);
  for k = 1:numel(ro)
    w = ro(k) * vector3d(A.'); img = reshape(w.xyz,[],3);   % images of a,b,c
    ops(k,:,:) = round((img / A).', 9);
  end
  out{end+1} = struct('group',nm,'ops',ops,'A',A,'laue',char(cs0.Laue.pointGroup), ...
    'lauerows',find(ismember(round(reshape(ops,[],9),6), round(reshape(lauemats(cs0,A),[],9),6),'rows')));
  rng(0); X = randn(500,3); X = X ./ vecnorm(X,2,2);   % fixed test directions (Cartesian, then uvw)
  out{end}.dirs_uvw = X / A;
  for p = 1:size(pcs,1)
    cs = crystalSymmetry(nm, cell{1}, cell{2}*degree);
    cs.frame.how2plot = pcs{p,2};
    e = struct('group',nm,'id',cs.id,'laue',char(cs.Laue.pointGroup),'abc',cell{1},'angles',cell{2}, ...
      'how2plot',pcs{p,1},'outOfScreen',L(pcs{p,2}.outOfScreen,A),'east',L(pcs{p,2}.east,A));
    % TSL: Laue group sector, centre, hue reference
    try
      k = TSLDirectionKey(cs.Laue); sR = k.sR;
      v = sR.vertices; if ~isempty(v), ref = v(1); else, ref = vector3d(cs.Laue.aAxisRec); end
      e.tsl = struct('N',L(sR.N,A),'vertices',L(sR.vertices,A),'center',L(sR.center,A), ...
        'ref',L(ref,A),'rx',L(rxOf(sR,sR.center,ref),A),'rgb',rgbOf(k,Miller(vector3d(X.'),cs.Laue)));
    catch err, e.tsl = err.message; end
    % HSV: point group sector + extra white/black planes, white centre, hue reference east
    try
      k = HSVDirectionKey(cs); sR = k.sR;
      n0 = numel(sR.N) - numel(k.refl);
      wC = vector3d(project2FundamentalRegion(k.whiteCenter,cs));
      e.hsv = struct('N',L(sR.N(1:n0),A),'refl',L(k.refl,A),'vertices',L(cs.fundamentalSector.vertices,A), ...
        'whiteCenter',L(wC,A),'ref',L(cs.how2plot.east,A),'rx',L(rxOf(sR,wC,vector3d(cs.how2plot.east)),A), ...
        'rgb',rgbOf(k,Miller(vector3d(X.'),cs)));
    catch err, e.hsv = err.message; end
    out{end+1} = e;
  end
  fprintf('%s done\n', nm);
end
fid = fopen('mtex72_keys.json','w'); fprintf(fid,'%s',jsonencode(out)); fclose(fid);

function M = lauemats(cs,A)
  L = cs.Laue.rot; M = zeros(numel(L),3,3);
  for k = 1:numel(L), w = L(k) * vector3d(A.'); M(k,:,:) = round((reshape(w.xyz,[],3) / A).', 9); end
end
function c = rgbOf(k,h)   % colours of h, or a message when MTEX returns the wrong shape
  c = k.direction2color(h);
  if numel(c) ~= 3*length(h), c = sprintf('direction2color returned %s', mat2str(size(c))); end
end
function s = sgn(x), if x > 0, s = '+'; else, s = '-'; end, end
function S = L(v,A)  % lattice components of vector3d v, one row per vector
  if isempty(v), S = struct('uvw',zeros(0,3),'hkl',zeros(0,3)); return; end
  x = reshape(v.xyz,[],3);
  S = struct('uvw', x / A, 'hkl', x * A.');
end
function rx = rxOf(sR,center,ref)   % as in polarCoordinates.m
  center = center.normalize;
  rx = sR.how2plot.outOfScreen - center;
  if tnull(rx,center), rx = ref - center; end
  if tnull(rx,center), rx = orth(center); end
  rx = normalize(rx - dot(rx,center) .* center);
end
function tf = tnull(rx,c), tf = norm(rx - dot(rx,c).*c) < 1e-10; end

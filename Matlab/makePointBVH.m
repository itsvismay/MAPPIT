 % Format:
% root: 0 [bboxMin] [bboxMax] rightSubIdx
% leaf: vidx [bboxMin] [bboxMax] 0
function T = makePointBVH(V, vidx)
  assert(numel(vidx) > 0);

  bboxMin = min(V(vidx, :),[],1);
  bboxMax = max(V(vidx, :),[],1);
  if numel(vidx) == 1
    % TODO: generalize for 3d later
    T = [vidx bboxMin bboxMax 0];
    return;
  end

  [~,largestAxis] = max(bboxMax - bboxMin);

  T = [0 bboxMin bboxMax 0];

  split = (bboxMin(largestAxis) + bboxMax(largestAxis)) / 2;
  inleft = find(arrayfun(@(i) any(V(i, largestAxis) < split), vidx));
  inright = find(arrayfun(@(i) any(V(i, largestAxis) >= split), vidx));
  straddlers = intersect(inleft, inright);
  inleftonly = setdiff(inleft, straddlers);
  inrightonly = setdiff(inright, straddlers);

  if numel(inleftonly) == 0 && numel(inrightonly) == 0
    % All triangles straddle - split arbitrarily between halves
    ns = floor(numel(straddlers) / 2);
    inleft = straddlers(1:ns);
    inright = straddlers(ns+1:end);
  elseif numel(inleftonly) == 0
    % No triangles in left subtree - give it the straddlers
    inleft = straddlers;
    inright = inrightonly;
  elseif numel(inrightonly) == 0
    % No triangles in right subtree - give it the straddlers
    inleft = inleftonly;
    inright = straddlers;
  else
    % Both sides nonempty - give the straddlers to the left subtree
    inright = inrightonly;
  end

  T = [ T ; makePointBVH(V, vidx(inleft)) ];
  T(1, end) = size(T, 1) + 1;
  T = [ T ; makePointBVH(V, vidx(inright)) ];
end

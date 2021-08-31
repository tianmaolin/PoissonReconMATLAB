function v = isoValueGrid(imp, samples, maxDepth)
% Maolin Tian, 2020
n = 2^maxDepth;
scale = 1 / n;
v = 0;
for s = 1:samples.Count
  pos = round(samples.Location(s, :)*n);
  if (pos(1) < 1 || pos(1) > n || pos(2) < 1 || pos(2) > n)
    warning("pos outrange in isoValueGrid!");
    continue;
  end
  v = v + imp(pos(1), pos(2));
end
v = v / samples.Count;

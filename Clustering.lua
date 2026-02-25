label = "Clustering" 

about = [[
  - K-means: clustering for selected mark symbols
  - K-means++: improved initialization for better clustering
  - DBScan:  Clustering method to group points based on how closely they are packed together.
  - HDBScan: Automatically detects and selects stable clusters using the Excess of Mass algorithm.
]]

revertOriginal = _G.revertOriginal


-- ========================================================================================
--                                      Support
-- ========================================================================================
-- Returns a list of vectors of all the selected points
function get_unique_selected_points(model)
  local page = model:page()
  local points = {}

  -- Loop through all objects on the current page
  for i, obj, sel, _ in page:objects() do
    -- Check if the object is selected and is a reference
    if sel and obj:type() == "reference" then
      local point = obj:matrix() * obj:position() 
      local replaced = false
      
      -- Check if there is already a point at the same coordinates (need the newest one to be able to change color of the updated reference)
      for j, check in ipairs(points) do 
        if check.point.x == point.x and check.point.y == point.y then 
          -- If duplicate, replace with newest selection
          points[j] = { idx = i, point = point, obj = obj} 
          replaced = true
          break
        end
      end
      
      -- If it's a new coordinate, add it to our list
      if not replaced then
        table.insert(points, { idx = i, point = point, obj = obj}) 
      end
    end
  end
  return points
end
--! CONVEX HULL (GRAHAM SCAN) -- from the library
-- https://www.codingdrills.com/tutorial/introduction-to-divide-and-conquer-algorithms/convex-hull-graham-scan

-- Function to calculate the squared distance between two points
function squared_distance(p1, p2)
    return (p1.x - p2.x)^2 + (p1.y - p2.y)^2
end

-- Function to compare two points with respect to a given 'lowest' point
-- Closure over the lowest point to create a compare function
function create_compare_function(lowest, model)
    return function(p1, p2) -- anonymous function

        -- Determine the orientation of the triplet (lowest, p1, p2)
        local o = orientation(lowest, p1, p2, model)

        -- If p1 and p2 are collinear with lowest, choose the farther one to lowest
        if o == 0 then
            return squared_distance(lowest, p1) < squared_distance(lowest, p2)
        end

        -- For non-collinear points, choose the one that forms a counterclockwise turn with lowest
        return o == 2
    end
end

-- Function to find the orientation of ordered triplet (p, q, r).
-- The function returns the following values:
-- 0 : Collinear points
-- 1 : Clockwise points
-- 2 : Counterclockwise  
function orientation(p, q, r, model)
    -- print the vectors and val
    -- print_vertices({p, q, r}, "Orientation", model)
    local val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    -- print(val, "Orientation", model)
    if val == 0 then return 0  -- Collinear
    elseif val > 0 then return 2  -- Counterclockwise
    else return 1  -- Clockwise
    end
end

function convex_hull(points, model)
    local n = #points
    if n < 3 then return {} end  -- Less than 3 points cannot form a convex hull

    -- Find the point with the lowest y-coordinate (or leftmost in case of a tie)
    local lowest = 1
    for i = 2, n do
        if points[i].y < points[lowest].y or (points[i].y == points[lowest].y and points[i].x < points[lowest].x) then
            lowest = i
        end
    end

    -- Swap the lowest point to the start of the array
    points[1], points[lowest] = points[lowest], points[1]

    -- Sort the rest of the points based on their polar angle with the lowest point
    local compare = create_compare_function(points[1], model) -- closure over the lowest point
    table.sort(points, compare)

    -- Sorted points are necessary but not sufficient to form a convex hull.
    --! The stack is used to maintain the vertices of the convex hull in construction.

    -- Initializing stack with the first three sorted points
    -- These form the starting basis of the convex hull.
    local stack = {points[1], points[2], points[3]}
    local non_stack = {}

    -- Process the remaining points to build the convex hull
    for i = 4, n do
        -- Check if adding the new point maintains the convex shape.
        -- Remove points from the stack if they create a 'right turn'.
        -- This ensures only convex shapes are formed.
        while #stack > 1 and orientation(stack[#stack - 1], stack[#stack], points[i]) ~= 2 do
            table.remove(stack)
        end
        table.insert(stack, points[i])  -- Add the new point to the stack
    end

    -- The stack now contains the vertices of the convex hull in counterclockwise order.
    return stack
end

function create_shape_from_vertices(v, model)
	local shape = {type="curve", closed=true;}
	for i=1, #v-1 do
		table.insert(shape, {type="segment", v[i], v[i+1]})
	end
  	table.insert(shape, {type="segment", v[#v], v[1]})
	return shape
end

function point_on_segment(p, s)
    local cross_product = (p.x - s[1].x) * (s[2].y - s[1].y) - (p.y - s[1].y) * (s[2].x - s[1].x)
    if cross_product ~= 0 then return false end

    local dot_product = (p.x - s[1].x) * (s[2].x - s[1].x) + (p.y - s[1].y) * (s[2].y - s[1].y)
    if dot_product < 0 then return false end
    if dot_product > squared_distance(s[1], s[2]) then return false end

    return true
end

function create_segments_from_vertices(vertices)
	local segments_start_finish = {}
	for i=1, #vertices-1 do
		table.insert( segments_start_finish, {vertices[i],vertices[i+1]} )
	end

	table.insert( segments_start_finish, {vertices[#vertices], vertices[1]} )
	return segments_start_finish
end


local function copy_list(list)
  local out = {}
  for i = 1, #list do
    out[i] = list[i]
  end
  return out
end


local function euclidean_distance(a, b)
  local sum = 0.0
  for i = 1, #a do
    local d = a[i] - b[i]
    sum = sum + d * d
  end
  return math.sqrt(sum)
end

local function validate_labels(n, labels)
  if labels == nil then
    local out = {}
    for i = 1, n do
      out[i] = i
    end
    return out
  end
  return copy_list(labels)
end

local function precompute_distances(points, distance)
  local n = #points
  local dist = {}
  for i = 1, n do
    dist[i] = {}
    dist[i][i] = 0.0
  end

  for i = 1, n do
    for j = i + 1, n do
      local d = distance(points[i], points[j])
      dist[i][j] = d
      dist[j][i] = d
    end
  end
  return dist
end

local function remove_pair(ids, i, j)
  table.remove(ids, j)
  table.remove(ids, i)
end
-- ========================================================================================
--                                      K-Means
-- ========================================================================================

-- K-means algo on plain 2d points arrays
function kmeans(points, k, max_iters, initial_centroids)
  max_iters = max_iters or 100
  local n = #points
  if n == 0 or k < 1 then return nil, nil end
  if k > n then k = n end

  local function dist2(a, b)
    local dx = a[1] - b[1]
    local dy = a[2] - b[2]
    return dx * dx + dy * dy
  end

  -- initialize centroids as random unique points
  local centroids = {}
  if initial_centroids then
    for i = 1, k do
      centroids[i] = { initial_centroids[i][1], initial_centroids[i][2] }
    end
  else
    for i = 1, k do
      centroids[i] = { points[i][1], points[i][2] }
    end
  end


  -- assign each point to the nearest centroid
  local assign = {}
  for iter = 1, max_iters do
    local changed = false
    -- assign
    for i = 1, n do
      local best, bestj = 1e309, 1
      for j = 1, k do
        local d = dist2(points[i], centroids[j])
        if d < best then best, bestj = d, j end
      end
      if assign [i] ~= bestj then
        assign[i] = bestj
        changed = true
      end
    end
    if not changed and iter > 1 then break end

    -- recalc centroids as the mean of each cluster
    local sumx, sumy, cnt = {}, {}, {}
    for j = 1, k do sumx[j], sumy[j], cnt[j] = 0, 0, 0 end
    for i = 1, n do
      local c = assign[i]
      sumx[c] = sumx[c] + points[i][1]
      sumy[c] = sumy[c] + points[i][2]
      cnt[c] = cnt[c] + 1
    end
    for j = 1, k do
      if cnt[j] > 0 then
        centroids[j][1] = sumx[j] / cnt[j]
        centroids[j][2] = sumy[j] / cnt[j]
      end
    end
  end
  return assign, centroids
end


-- ========================================================================================
--                                      K-Means++
-- ========================================================================================

-- K-means++ initialization: choose centroids that are far apart
function kmeanspp_init(points, k)
  local n = #points
  if n == 0 or k < 1 then return nil end
  if k > n then k = n end
  
  local function dist2(a, b)
    local dx = a[1] - b[1]
    local dy = a[2] - b[2]
    return dx * dx + dy * dy
  end
  
  local centroids = {}
  local used = {}
  
  -- choose first centroid randomly
  local first_idx = math.random(1, n)
  centroids[1] = { points[first_idx][1], points[first_idx][2] }
  used[first_idx] = true
  
  -- choose remaining centroids with probability proportional to D(x)^2 whatever that means
  for c = 2, k do
    local distances = {}
    local sum_dist = 0
    
    -- for each point find distance to nearest centroid
    for i = 1, n do
      if not used[i] then
        local min_dist = 1e309
        for j = 1, #centroids do
          local d = dist2(points[i], centroids[j])
          if d < min_dist then min_dist = d end
        end
        distances[i] = min_dist
        sum_dist = sum_dist + min_dist
      else
        distances[i] = 0
      end
    end
    
    -- choose next centroid with a weighted probability
    if sum_dist > 0 then
      local threshold = math.random() * sum_dist
      local cumsum = 0
      for i = 1, n do
        if not used[i] then
          cumsum = cumsum + distances[i]
          if cumsum >= threshold then
            centroids[c] = { points[i][1], points[i][2] }
            used[i] = true
            break
          end
        end
      end
    else
      -- if all remaining points are at the same location
      for i = 1, n do
        if not used[i] then
          centroids[c] = { points[i][1], points[i][2] }
          used[i] = true
          break
        end
      end
    end
  end
  -- there is prolly a way to optimize this, but ill do it later
  return centroids
end



-- simple prompt for k using ipe's UI
local function prompt_k(model, default_k)
  local s = model:getString("Enter number of clusters (k)")
  if not s or s:match("^%s*$") then return nil end
  local k = tonumber(s)
  if not k or k < 1 then
    model:warning("Please enter a positive integer for k.")
    return nil
  end
  return math.floor(k)
end

-- extract positions of selected mark symbols
local function collect_selected_points(page)
  local points, indices = {}, {}
  
  for i, obj, sel, layer in page:objects() do
    if sel and obj:type() == "reference" then
     local v = obj:matrix() * obj:position()
      points[#points + 1] = { v.x, v.y }
      indices[#indices + 1] = i 
    end
  end
  return points, indices
end

-- Action: K-means clustering and recolor selected marks
local function cmd_kmeans(model)
  local page = model:page()
  local points, sel_indices = collect_selected_points(page)

  if #points == 0 then
    model:warning("Select mark symbols to cluster!")
    return
  end

  local k = prompt_k(model, math.min(3, #points))
  if not k then return end
  if k > #points then
    model:warning("K must be <= number of selected marks!")
    return
  end

  local assign, centroids = kmeans(points, k, 100)
  if not assign then
    model:warning("K-means failed.")
    return
  end

  -- yay colors (I have brain damage) 
  local colors = { "red", "blue", "green", "orange", "purple", "cyan", "brown", "magenta", "navy", "gray" }
  -- according to the images I saw, its supposed to be a cricle, but thats really hard
  -- also this part was by far the biggest headache, it just did not want to work
  local t = { label = "K-means: color selected marks by cluster",
              pno = model.pno,
              vno = model.vno,
              original = model:page():clone(),
              final = model:page():clone(),
              undo = revertOriginal,
              redo = _G.revertFinal }
  
  local idx = 1
  for i, obj, sel, layer in t.final:objects() do
    if sel and obj:type() == "reference" then
      local c = assign[idx]
      local color = colors[((c - 1) % #colors) + 1]
      obj:set("stroke", color)
      obj:set("fill", color)
      idx = idx + 1
    end
  end
  
  model:register(t) 
  
  
  for j = 1, k do
    local color = colors[((j - 1) % #colors) + 1]
    local stroke = ipe.Reference(model.attributes, "mark/cross(sx)", ipe.Vector(centroids[j][1] ,centroids[j][2]))
    stroke:set("stroke", color)
    model:creation("Centroid", stroke)
  end
  
  local clusterPoints = {}
  
  for i=1, k do
    clusterPoints[i]={}
    for j, p in ipairs(points) do
      if assign[j]==i then
        table.insert(clusterPoints[i], { x = p[1], y = p[2] })
      end
    end
  end
  
  for j=1, k do
    local color = colors[((j - 1) % #colors) + 1]
    if #clusterPoints[j]>= 3 then
      local convexHull= convex_hull(clusterPoints[j],model)
      local vertices = {}
      for _, p in ipairs(convexHull) do
        table.insert(vertices, ipe.Vector(p.x, p.y))
      end
      local shape = create_shape_from_vertices(vertices, model)
      local path = ipe.Path(model.attributes, { shape })
      path:set("stroke", color)
      path:set("dashstyle", "dashed")
      model:creation("Cluster hull", path)
    end
  end
end

-- Action: K-means++ or something im tried
local function cmd_kmeanspp(model)
  local page = model:page()
  local points, sel_indices = collect_selected_points(page)

  if #points == 0 then
    model:warning("Select mark symbols to cluster!")
    return
  end

  local k = prompt_k(model, math.min(3, #points))
  if not k then return end
  if k > #points then
    model:warning("K must be <= number of selected marks!")
    return
  end

  -- use K-means++ initialization
  local initial_centroids = kmeanspp_init(points, k)
  if not initial_centroids then
    model:warning("K-means++ initialization failed.")
    return
  end

  -- run K-means with the improved initialization
  local assign, centroids = kmeans(points, k, 100, initial_centroids)
  if not assign then
    model:warning("K-means++ failed.")
    return
  end

  -- yay colors again (its 2am and I have CS exam at 9)
  local colors = { "red", "blue", "green", "orange", "purple", "cyan", "brown", "magenta", "navy", "gray" }
  -- too lazy to code a srko, you get colorfull dots intead, have a feild day
  local t = { label = "K-means++: color selected marks by cluster",
              pno = model.pno,
              vno = model.vno,
              original = model:page():clone(),
              final = model:page():clone(),
              undo = revertOriginal,
              redo = _G.revertFinal }
  
  local idx = 1
  for i, obj, sel, layer in t.final:objects() do
    if sel and obj:type() == "reference" then
      local c = assign[idx]
      local color = colors[((c - 1) % #colors) + 1]
      obj:set("stroke", color)
      obj:set("fill", color)
      idx = idx + 1
    end
  end
  
  model:register(t)
  
  for j = 1, k do
    local color = colors[((j - 1) % #colors) + 1]
    local stroke = ipe.Reference(model.attributes, "mark/cross(sx)", ipe.Vector(centroids[j][1] ,centroids[j][2]))
    stroke:set("stroke", color)
    model:creation("Centroid", stroke)
  end
  
  local clusterPoints = {}
  
  for i=1, k do
    clusterPoints[i]={}
    for j, p in ipairs(points) do
      if assign[j]==i then
        table.insert(clusterPoints[i], { x = p[1], y = p[2] })
      end
    end
  end
  
  for j=1, k do
    local color = colors[((j - 1) % #colors) + 1]
    if #clusterPoints[j]>= 3 then
      local convexHull= convex_hull(clusterPoints[j],model)
      local vertices = {}
      for _, p in ipairs(convexHull) do
        table.insert(vertices, ipe.Vector(p.x, p.y))
      end
      local shape = create_shape_from_vertices(vertices, model)
      local path = ipe.Path(model.attributes, { shape })
      path:set("stroke", color)
      path:set("dashstyle", "dashed")
      model:creation("Cluster hull", path)
    end
  end
end



-- ========================================================================================
--                                      K-Mediods
-- ========================================================================================



-- Helper function to calculate Euclidean distance between two points
local function dist(p1, p2)
  local dx = p1[1] - p2[1]
  local dy = p1[2] - p2[2]
  return math.sqrt(dx * dx + dy * dy)
end

-- based on the algorithm by Park and Jun
-- described here https://www.sciencedirect.com/science/article/pii/S095741740800081X
function kmedoids(points, k, max_iters)
  max_iters = max_iters or 100
  local n = #points
  if n == 0 or k < 1 then return nil, nil end
  if k > n then k = n end

  -- calculate the distance matrix D
  local D = {}
  for i = 1, n do
    D[i] = {}
    for j = 1, n do
      D[i][j] = dist(points[i], points[j])
    end
  end

  -- Calculate vj for object j
  local v = {}
  for j = 1, n do
    local vj_sum = 0
    for i = 1, n do
      local sum_dil = 0
      for l = 1, n do
        sum_dil = sum_dil + D[i][l]
      end
      -- avoid division by zero
      if sum_dil > 0 then
        vj_sum = vj_sum + (D[i][j] / sum_dil)
      end
    end
    table.insert(v, {index = j, value = vj_sum})
  end

  -- Sort vj's in ascending order and select first k smallest for initial medoids
  table.sort(v, function(a, b) return a.value < b.value end)
  
  local medoids = {}
  for idx = 1, k do
    table.insert(medoids, v[idx].index)
  end

  -- Helper closure to assign points to clusters and sum costs
  local function assign_and_calculate_cost(current_medoids)
    local assign = {}
    local clusters = {}
    for _, m in ipairs(current_medoids) do
      clusters[m] = {}
    end
    
    local total_cost = 0
    
    for i = 1, n do
      local min_dist = math.huge
      local nearest_medoid = -1
      local cluster_idx = 1
      
      for c, m in ipairs(current_medoids) do
        if D[i][m] < min_dist then
          min_dist = D[i][m]
          nearest_medoid = m
          cluster_idx = c -- Store 1 to K for color mapping later
        end
      end
      table.insert(clusters[nearest_medoid], i)
      assign[i] = cluster_idx 
      total_cost = total_cost + min_dist
    end
    
    return assign, clusters, total_cost
  end

  local assign, clusters, previous_cost = assign_and_calculate_cost(medoids)

  -- Steps 2 and 3 of Park and Jun's algorithm
  for iter = 1, max_iters do

    -- Update medoids
    local new_medoids = {}
    
    for _, m in ipairs(medoids) do
      local cluster_objects = clusters[m]
      local best_new_medoid = m
      local min_total_cluster_dist = math.huge
      
      if cluster_objects then
        for _, candidate in ipairs(cluster_objects) do
          local candidate_dist_sum = 0
          for _, other_obj in ipairs(cluster_objects) do
            candidate_dist_sum = candidate_dist_sum + D[candidate][other_obj]
          end
          
          if candidate_dist_sum < min_total_cluster_dist then
            min_total_cluster_dist = candidate_dist_sum
            best_new_medoid = candidate
          end
        end
      end
      table.insert(new_medoids, best_new_medoid)
    end
    
    medoids = new_medoids
    
    -- reassign objects to medoids
    local current_assign, current_clusters, current_cost = assign_and_calculate_cost(medoids)
    
    if current_cost == previous_cost then
      assign = current_assign
      break
    else
      assign = current_assign
      clusters = current_clusters
      previous_cost = current_cost
    end
  end

  -- Map medoid indices back to actual point coordinates for Ipe drawing
  local medoids_coords = {}
  for i = 1, k do
    medoids_coords[i] = { points[medoids[i]][1], points[medoids[i]][2] }
  end

  return assign, medoids_coords
end

local function cmd_kmedoids(model)
  local page = model:page()
  local points, sel_indices = collect_selected_points(page)

  if #points == 0 then
    model:warning("Select mark symbols to cluster!")
    return
  end

  local k = prompt_k(model, math.min(3, #points))
  if not k then return end
  if k > #points then
    model:warning("K must be <= number of selected marks!")
    return
  end

  local assign, medoids = kmedoids(points, k, 100)
  if not assign then
    model:warning("K-medoids failed")
    return
  end

  local colors = { "red", "blue", "green", "orange", "purple", "cyan", "brown", "magenta", "navy", "gray" }
  
  local t = { label = "K-medoids: color selected marks by cluster",
              pno = model.pno,
              vno = model.vno,
              original = model:page():clone(),
              final = model:page():clone(),
              undo = revertOriginal,
              redo = _G.revertFinal }
  
  local idx = 1
  for i, obj, sel, layer in t.final:objects() do
    if sel and obj:type() == "reference" then
      local c = assign[idx]
      local color = colors[((c - 1) % #colors) + 1]
      obj:set("stroke", color)
      obj:set("fill", color)
      idx = idx + 1
    end
  end
  
  model:register(t) 
  
  for j = 1, k do -- Draws a cross at the point of each cluster's mediod
    local color = colors[((j - 1) % #colors) + 1]
    local stroke = ipe.Reference(model.attributes, "mark/cross(sx)", ipe.Vector(medoids[j][1] ,medoids[j][2]))
    stroke:set("stroke", color)
    model:creation("Medoid", stroke)
  end
  
  -- Gather points for convex hulls of clusters
  local clusterPoints = {}
  for i=1, k do
    clusterPoints[i]={}
  end
  
  for j, p in ipairs(points) do
    local cluster_idx = assign[j]
    if cluster_idx then
      table.insert(clusterPoints[cluster_idx], { x = p[1], y = p[2] })
    end
  end
  
  for j=1, k do -- Draw Convex Hulls
    local color = colors[((j - 1) % #colors) + 1]
    if #clusterPoints[j] >= 3 then
      local convexHull = convex_hull(clusterPoints[j], model)
      local vertices = {}
      for _, p in ipairs(convexHull) do
        table.insert(vertices, ipe.Vector(p.x, p.y))
      end
      
      if #vertices > 0 then -- Failsafe for degenerate hulls
        local shape = create_shape_from_vertices(vertices, model)
        local path = ipe.Path(model.attributes, { shape })
        path:set("stroke", color)
        path:set("dashstyle", "dashed")
        model:creation("Cluster hull", path)
      end
    end
  end
end

-- ========================================================================================
--                                      DBScan
-- ========================================================================================

-- Finding distance between two points
local function distance(p1, p2)
  local dx = p1.point.x - p2.point.x
  local dy = p1.point.y - p2.point.y
  return math.sqrt(dx*dx + dy*dy)
end

-- Find neighbors within eps
local function RangeQuery(points, corePoint, eps)
  local neighbors = {}

  for i, point in ipairs(points) do
    if distance(corePoint, point) <= eps then
      table.insert(neighbors, point)
    end
  end

  return neighbors
end

-- Core DBSCAN Logic Function
local function DBSCAN(points, eps, minPts)
  -- Counter for clusters
  local clusterCount = 0  
  
  -- Array to label each point as its cluster ID or "Noise"
  local labels = {}                      
  
  -- Table of the clusters
  local clusters = {}                      

  -- Initialize every point in labels array 
  for i = 1, #points do
    labels[i] = nil
  end

  for i, P in ipairs(points) do
    -- Skip if already processed
    if labels[i] ~= nil then
      goto continue
    end

    -- Find neighbors of P
    local neighbors = RangeQuery(points, P, eps)

    -- Checking if neighbors can be a core point via its density
    if #neighbors < minPts then
      labels[i] = "Noise"
      goto continue
    end

    -- Start a new cluster
    clusterCount = clusterCount + 1
    clusters[clusterCount] = {}
    labels[i] = clusterCount
    table.insert(clusters[clusterCount], P)

    -- All neighbors in this the cluster of this core point are inserted into the seedSet
    local seedSet = {}
    for idx, Q in ipairs(neighbors) do
      if Q ~= P then
        table.insert(seedSet, Q)
      end
    end

    -- Expand cluster
    local j = 1

    -- Iterating over the seedSet
    while j <= #seedSet do
      local Q = seedSet[j]

      -- Find index of Q in points
      local q_index = nil
      for k, pt in ipairs(points) do
        if pt == Q then
          q_index = k
          break
        end
      end

      if labels[q_index] == "Noise" then
        labels[q_index] = clusterCount  -- convert Noise to border point
      end

      if labels[q_index] == nil then
        labels[q_index] = clusterCount
        table.insert(clusters[clusterCount], Q)

        -- Find neighbors of Q
        local neighborsQ = RangeQuery(points, Q, eps)

        -- If Q is a core point, then add its neighbors to the seedSet
        if #neighborsQ >= minPts then
          -- Add new neighbors to seed set (union)
          for _, npt in ipairs(neighborsQ) do
            local found = false
            for _, spt in ipairs(seedSet) do
              if spt == npt then
                found = true
                break
              end
            end
            -- Making sure that the neighbors of Q being added are not duplicates (already neighbors of seedSet)
            if not found then
              table.insert(seedSet, npt)
            end
          end
        end
      end

      j = j + 1
    end

    ::continue::
  end

  return labels, clusters
end



function dbscan(model)
  local points = get_unique_selected_points(model)
  if #points == 0 then
    model:warning("No reference points selected!")
    return
  end


  local stringEps = model:getString("Enter epsilon (each square has length 16)")
  local stringMinPts = model:getString("Enter min Pts")
  local eps = math.tointeger(stringEps)
  local minPts = math.tointeger(stringMinPts)

   -- Run DBSCAN to get cluster assignments
  local labels, clusters = DBSCAN(points, eps, minPts)
  local colors = {"red","blue","green","orange","purple","brown","navy","pink"}

  model:warning("Clusters found: " .. #clusters)

  -- Applying colors to points in their respective clusters
  for i, pt in ipairs(points) do
    if labels[i] ~= "Noise" and labels[i] > 0 then
      -- Pick a color based on Cluster ID
      pt.obj:set("stroke", colors[(labels[i]-1) % #colors + 1])
      model:page():replace(pt.idx, pt.obj)
    else
      -- Mark outliers as gray
      pt.obj:set("stroke", "gray")
      model:page():replace(pt.idx, pt.obj)
    end
  end
  
  local clusterPoints = {}
  
  for i=1, #clusters do
    clusterPoints[i]={}
    for _, p in ipairs(clusters[i]) do
      if p and p.point and p.point.x and p.point.y then
        table.insert(clusterPoints[i], { x = p.point.x, y =  p.point.y})
      end
    end
  end
  
  for j=1, #clusters do
    local color = colors[((j - 1) % #colors) + 1]
    if #clusterPoints[j]>= 3 then
      local convexHull= convex_hull(clusterPoints[j],model)
      local vertices = {}
      for _, p in ipairs(convexHull) do
        table.insert(vertices, ipe.Vector(p.x, p.y))
      end
      local shape = create_shape_from_vertices(vertices, model)
      local path = ipe.Path(model.attributes, { shape })
      path:set("stroke", color)
      path:set("dashstyle", "dashed")
      model:creation("Cluster hull", path)
    end
  end

end

-- ========================================================================================
--                                      HDBScan
-- ========================================================================================
--[[ 
    ALGORITHM PSEUDOCODE (HDBSCAN):
    1. Compute the core distance w.r.t. mpts for all data objects in X. 
    2. Compute an MST of Gmpts, the Mutual Reachability Graph. 
    3. Extend the MST to obtain MSText, by adding for each vertex a “self edge” with the core distance.
    4. Extract the HDBSCAN hierarchy as a dendrogram from MSText:
        4.1 For the root of the tree assign all objects the same label.
        4.2 Iteratively remove all edges from MSText in decreasing order of weights.
          4.2.1 Before removal, set dendrogram scale (lambda) as the weight of the edge.
         4.2.2 Assign labels to components; if it has at least one edge, it's a cluster, else "noise".
    5. (Post-Processing) Use "Excess of Mass" to select the most stable clusters from the hierarchy.
--]]


-- Core HDBSCAN Logic Function
-- Parameters of points (all selected points), and min_size (min # of points to form a cluster)
local function HDBSCAN(points, min_size)
  -- This table will store the 'core distance' for every point
  local core = {}
  
  -- [Step 1]: Compute the core distance for all data objects
  for i=1, #points do

    -- Temp list to store distances to all other points
    local d = {} 
    for j=1, #points do 
      -- Not calculating distance to itself
      if i ~= j then table.insert(d, distance(points[i], points[j])) 
      end 
    end

    -- Sort distances from smallest to largest
    table.sort(d)

    -- The core distance is the distance to the neighbor at the 'min_size' position.
    -- If there aren't enough points, we fall back to the last point or 0.
    core[i] = d[min_size] or d[#d] or 0
  end

  -- [Step 2 & 3 Prep]: Mutual Reachability Distance
  local function mr(i, j) 
    -- Return max of Point A core distance, Point B core distance, and actual distance between the two points
    return math.max(core[i], core[j], distance(points[i], points[j])) 
  end

  -- [Step 2]: Minimum Spanning Tree 
  local edges, in_tree, min_w, close = {}, {}, {}, {}
  
  -- Initialize distances to 'infinity' for all points
  for i=2, #points do min_w[i]=math.huge 
  end

  -- Building the tree from the first point
  in_tree[1] = true
  local curr = 1
  
  -- Loop until every point is connected to the tree
  for _=1, #points-1 do
    local best, val = -1, math.huge
    for v=1, #points do
      -- If the current point is not yet in the growing tree
      if not in_tree[v] then
        -- Calculate Mutual Reachability distance from our current point to 'v'
        local w = mr(curr, v)

        -- If this path is shorter than the shortest known path to 'v', update it
        if w < (min_w[v] or math.huge) then 
          min_w[v]=w; 
          close[v]=curr 
        end

        -- Keep track of the absolute shortest edge available to any point not in the tree
        if min_w[v] < val then 
          val=min_w[v]
          best=v 
        end
      end
    end

    -- If no path is found, the points are disconnected
    if best == -1 then 
      break 
    end
    
    -- Add the best point found to the tree
    in_tree[best] = true; 
    curr = best

    -- Store the edge details: point 'u' to point 'v' with weight 'w'
    table.insert(edges, {u=close[best], v=best, w=val})
  end
  
  -- Sort all MST edges by weight (ascending) 
  table.sort(edges, function(a,b) return a.w < b.w end) 

  -- [Step 4]: Extract the HDBSCAN hierarchy (Dendrogram)

  local nodes = {} 
  -- Initialize every point as its own individual leaf node
  for i=1, #points do 
    nodes[i] = {
      -- Initila size of 1 point
      size=1,           
      stab=0,         
      -- Lamda starts at the largest number since it is 1 / distance
      lambda=math.huge, 
      pts={i}         
    } 
  end

  -- Tracks which cluster a point belongs to
  local parent = {}

  -- Follows the parent chain to find the root of a point's cluster
  local function find(i) 
    if parent[i] then 
      return find(parent[i]) 
    end 
    return i 
  end

  -- Iterate through sorted edges to merge clusters
  for _, e in ipairs(edges) do
    local r1, r2 = find(e.u), find(e.v)
    -- If the two points of the edge are currently in different clusters:
    if r1 ~= r2 then
      -- lambda is the inverse of distance. Higher lambda = more dense/closer.
      local lambda = 1/(e.w + 1e-9)
      local n1 = nodes[r1]
      local n2 = nodes[r2]
      
      -- Calculate stability (how much density a cluster)
      local function update_stab(n, d_lambda)
         -- Only count stability if the cluster meets the minimum size requirement
         return (n.size >= min_size) and ((n.lambda - d_lambda) * n.size) or 0
      end
      
      -- Update stability for the two branches being merged
      n1.stab = n1.stab + update_stab(n1, lambda)
      n2.stab = n2.stab + update_stab(n2, lambda)

      -- Create a new Parent node representing the merged cluster
      local new_idx = #nodes + 1
      nodes[new_idx] = { 
        -- Combined size of both clusters
        size = n1.size + n2.size, 
        stab = 0, 
        lambda = lambda, 
        children = {r1, r2} 
      }

      -- Update the table so both children now point to the new parent
      parent[r1], parent[r2] = new_idx, new_idx
    end
  end

  -- [Step 5]: Stability Optimization (finding the real clusters by comparing a parents stability with the childrens)
  local function select_clusters(idx)
    local n = nodes[idx]
    -- If it's a leaf node (individual point), its stability is as is.
    if not n.children then 
      n.is_selected_stab = n.stab; 
      return 
    end
    
    -- Recursively check the children first
    select_clusters(n.children[1])
    select_clusters(n.children[2])
    
    -- Sum of the stability of the children
    local child_sum = nodes[n.children[1]].is_selected_stab + nodes[n.children[2]].is_selected_stab
    
    -- If the children are collectively more stable than the parent, keep the children.
    if child_sum > n.stab then 
      n.is_selected_stab = child_sum
      -- Don't pick this paretn
      n.picked = false 
    else 
      -- If the parent is more stable than its children, this parent is a cluster.
      n.is_selected_stab = n.stab
      n.picked = true  
    end
  end
  
  select_clusters(#nodes)

  -- Assign Labels
  -- Prepare a list to store the final ID for every point
  local labels = {}
  for i=1,#points do labels[i] = "Noise" end
  local clusters = {}

  -- Helper function to drill down from a "picked" cluster node to find all individual points
  local function collect_leaves(idx, cluster_id)
    local n = nodes[idx]
    if not n.children then
      -- If we hit a leaf point, assign it the cluster ID
      local p_idx = n.pts[1]
      labels[p_idx] = cluster_id
      table.insert(clusters[cluster_id], points[p_idx])
    else
      -- Keep going down the tree
      collect_leaves(n.children[1], cluster_id)
      collect_leaves(n.children[2], cluster_id)
    end
  end

  -- Traverse the tree to find nodes we flagged as picked during optimization
  local function traverse(idx)
    local n = nodes[idx]
    if n.picked and n.size >= min_size then
      -- If this node is picked create a new cluster and collect all its points
      table.insert(clusters, {})
      collect_leaves(idx, #clusters)
    elseif n.children then
      -- If not picked, look at its children to see if they were picked
      traverse(n.children[1]); 
      traverse(n.children[2])
    end
  end
  
  -- If the root itself is picked, start there, otherwise, start from its children
  if nodes[#nodes].picked then 
    traverse(#nodes)
  elseif nodes[#nodes].children then 
    traverse(nodes[#nodes].children[1]); 
    traverse(nodes[#nodes].children[2]) 
  end

  return labels, clusters
end
  

function hdbscan(model)
  local points = get_unique_selected_points(model)
  
  if #points < 2 then
    model:warning("Select at least 2 points!")
    return
  end

  local minSize = 4

   -- Run HDBSCAN to get cluster assignments
  local labels, clusters = HDBSCAN(points, minSize)
  local colors = {"red","blue","green","orange","purple","brown","navy","pink"}

  model:warning("Clusters found: " .. #clusters)

  -- Applying colors to points in their respective clusters
  for i, pt in ipairs(points) do
    if labels[i] ~= "Noise" and labels[i] > 0 then
      -- Pick a color based on Cluster ID
      pt.obj:set("stroke", colors[(labels[i]-1) % #colors + 1])
      model:page():replace(pt.idx, pt.obj)
    else
      -- Mark outliers as gray
      pt.obj:set("stroke", "gray")
      model:page():replace(pt.idx, pt.obj)
    end
  end
  
  local clusterPoints = {}
  
  for i=1, #clusters do
    clusterPoints[i]={}
    for _, p in ipairs(clusters[i]) do
      if p and p.point and p.point.x and p.point.y then
        table.insert(clusterPoints[i], { x =  p.point.x, y = p.point.y })
      end
    end
  end
  
  for j=1, #clusters do
    local color = colors[((j - 1) % #colors) + 1]
    if #clusterPoints[j]>= 3 then
      local convexHull= convex_hull(clusterPoints[j],model)
      local vertices = {}
      for _, p in ipairs(convexHull) do
        table.insert(vertices, ipe.Vector(p.x, p.y))
      end
      local shape = create_shape_from_vertices(vertices, model)
      local path = ipe.Path(model.attributes, { shape })
      path:set("stroke", color)
      path:set("dashstyle", "dashed")
      model:creation("Cluster hull", path)
    end
  end
end

-- ========================================================================================
--                                      Complete Linkage
-- ========================================================================================
local function cmd_hierarchical(model, method)
  local page = model:page()
  local points, sel_indices = collect_selected_points(page)

  if #points == 0 then
    model:warning("Select mark symbols to cluster!")
    return
  end

  local k = prompt_k(model, math.min(3, #points))
  if not k then return end
  if k > #points then
    model:warning("K must be <= number of selected marks!")
    return
  end

  local tree = method.cluster(points)
  local assign = method.cut(tree, k)

  local colors = { "red", "blue", "green", "orange", "purple", "cyan", "brown", "magenta", "navy", "gray" }
  local t = { label = "Hierarchical",
              pno = model.pno,
              vno = model.vno,
              original = model:page():clone(),
              final = model:page():clone(),
              undo = revertOriginal,
              redo = _G.revertFinal }
  
  local idx = 1
  for i, obj, sel in t.final:objects() do
    if sel and obj:type() == "reference" then
      obj:set("stroke", colors[((assign[idx] - 1) % #colors) + 1])
      obj:set("fill", colors[((assign[idx] - 1) % #colors) + 1])
      idx = idx + 1
    end
  end
  model:register(t)
  
  local clusterPoints = {}
  for i = 1, k do clusterPoints[i] = {} end
  for j, p in ipairs(points) do
    table.insert(clusterPoints[assign[j]], { x = p[1], y = p[2] })
  end
  for j = 1, k do
    if #clusterPoints[j] >= 3 then
      local hull = convex_hull(clusterPoints[j], model)
      if #hull >= 3 then
        local verts = {}
        for _, p in ipairs(hull) do table.insert(verts, ipe.Vector(p.x, p.y)) end
        local path = ipe.Path(model.attributes, { create_shape_from_vertices(verts, model) })
        path:set("stroke", colors[((j - 1) % #colors) + 1])
        path:set("dashstyle", "dashed")
        model:creation("Cluster hull", path)
      end
    end
  end
end

local CompleteLinkage = {}

local function cluster_distance_complete(clusters, point_dist, id_a, id_b)
  local members_a = clusters[id_a].members
  local members_b = clusters[id_b].members
  local best = -math.huge
  for i = 1, #members_a do
    local pa = members_a[i]
    for j = 1, #members_b do
      local pb = members_b[j]
      local d = point_dist[pa][pb]
      if d > best then
        best = d
      end
    end
  end
  return best
end

function CompleteLinkage.cluster(points, opts)
  opts = opts or {}

  local n = #points
  local distance = opts.distance or euclidean_distance
  local labels = validate_labels(n, opts.labels)
  

  local point_dist = precompute_distances(points, distance)
  local clusters = {}
  local active_ids = {}
  local merges = {}

  for i = 1, n do
    clusters[i] = { id = i, members = { i } }
    active_ids[i] = i
  end

  local next_cluster_id = n + 1
  
  while #active_ids > 1 do
    local best_i = nil
    local best_j = nil
    local best_a = nil
    local best_b = nil
    local best_d = math.huge

    for i = 1, #active_ids - 1 do
      for j = i + 1, #active_ids do
        local id_a = active_ids[i]
        local id_b = active_ids[j]
        if id_a > id_b then
          id_a, id_b = id_b, id_a
        end

        local d = cluster_distance_complete(clusters, point_dist, id_a, id_b)
        if d < best_d - 1e-12 or
          (math.abs(d - best_d) <= 1e-12 and better_pair(id_a, id_b, best_a, best_b)) then
          best_i = i
          best_j = j
          best_a = id_a
          best_b = id_b
          best_d = d
        end
      end
    end

    local left_members = clusters[best_a].members
    local right_members = clusters[best_b].members
    local merged_members = {}
    for i = 1, #left_members do
      merged_members[#merged_members + 1] = left_members[i]
    end
    for i = 1, #right_members do
      merged_members[#merged_members + 1] = right_members[i]
    end
    table.sort(merged_members)

    local new_id = next_cluster_id
    next_cluster_id = next_cluster_id + 1
    clusters[new_id] = { id = new_id, members = merged_members }
    merges[#merges + 1] = {
      left = best_a,
      right = best_b,
      distance = best_d,
      size = #merged_members,
      id = new_id,
    }

    remove_pair(active_ids, best_i, best_j)
    active_ids[#active_ids + 1] = new_id
  end

  return {
    method = "complete",
    n = n,
    labels = labels,
    merges = merges,
    root_id = active_ids[1],
  }
end

function CompleteLinkage.cut(tree, k)

  local n = tree.n
  local merges_to_apply = n - k
  local clusters = {}

  for i = 1, n do
    clusters[i] = { i }
  end

  for step = 1, merges_to_apply do
    local m = tree.merges[step]
    local left = clusters[m.left]
    local right = clusters[m.right]

    local merged = {}
    for i = 1, #left do
      merged[#merged + 1] = left[i]
    end
    for i = 1, #right do
      merged[#merged + 1] = right[i]
    end
    table.sort(merged)

    clusters[m.left] = nil
    clusters[m.right] = nil
    clusters[m.id] = merged
  end

  local ordered = {}
  for _, members in pairs(clusters) do
    ordered[#ordered + 1] = { min_index = members[1], members = members }
  end

  table.sort(ordered, function(a, b)
    return a.min_index < b.min_index
  end)

  local assignments = {}
  local groups = {}
  for cluster_id = 1, #ordered do
    local members = ordered[cluster_id].members
    groups[cluster_id] = {}
    for i = 1, #members do
      local idx = members[i]
      assignments[idx] = cluster_id
      groups[cluster_id][#groups[cluster_id] + 1] = tree.labels[idx]
    end
  end
  return assignments, groups
end



-- ========================================================================================
--                                      Single Linkage
-- ========================================================================================

local SingleLinkage = {}

local function better_pair(a1, a2, b1, b2)
  return a1 < b1 or (a1 == b1 and a2 < b2)
end


local function cluster_distance_single(clusters, point_dist, id_a, id_b)
  local members_a = clusters[id_a].members
  local members_b = clusters[id_b].members
  local best = math.huge
  
  for i = 1, #members_a do
    local pa = members_a[i]
    for j = 1, #members_b do
      local pb = members_b[j]
      local d = point_dist[pa][pb]
      if d < best then
        best = d
      end
    end
  end
  return best
end

function SingleLinkage.cluster(points, opts)
  opts = opts or {}
  local n = #points
  local distance = opts.distance or euclidean_distance
  local labels = validate_labels(n, opts.labels)
  local point_dist = precompute_distances(points, distance)
  local clusters = {}
  local active_ids = {}
  local merges = {}

  for i = 1, n do
    clusters[i] = { id = i, members = { i } }
    active_ids[i] = i
  end

  local next_cluster_id = n + 1

  while #active_ids > 1 do
    local best_i = nil
    local best_j = nil
    local best_a = nil
    local best_b = nil
    local best_d = math.huge

    for i = 1, #active_ids - 1 do
      for j = i + 1, #active_ids do
        local id_a = active_ids[i]
        local id_b = active_ids[j]
        if id_a > id_b then
          id_a, id_b = id_b, id_a
        end

        local d = cluster_distance_single(clusters, point_dist, id_a, id_b)
        if d < best_d - 1e-12 or
          (math.abs(d - best_d) <= 1e-12 and better_pair(id_a, id_b, best_a, best_b)) then
          best_i = i
          best_j = j
          best_a = id_a
          best_b = id_b
          best_d = d
        end
      end
    end

    local left_members = clusters[best_a].members
    local right_members = clusters[best_b].members
    local merged_members = {}
    for i = 1, #left_members do
      merged_members[#merged_members + 1] = left_members[i]
    end
    for i = 1, #right_members do
      merged_members[#merged_members + 1] = right_members[i]
    end
    table.sort(merged_members)

    local new_id = next_cluster_id
    next_cluster_id = next_cluster_id + 1
    clusters[new_id] = { id = new_id, members = merged_members }
    merges[#merges + 1] = {
      left = best_a,
      right = best_b,
      distance = best_d,
      size = #merged_members,
      id = new_id,
    }

    remove_pair(active_ids, best_i, best_j)
    active_ids[#active_ids + 1] = new_id
  end

  return {
    method = "single",
    n = n,
    labels = labels,
    merges = merges,
    root_id = active_ids[1],
  }
end

function SingleLinkage.cut(tree, k)

  local n = tree.n
  local merges_to_apply = n - k
  local clusters = {}

  for i = 1, n do
    clusters[i] = { i }
  end

  for step = 1, merges_to_apply do
    local m = tree.merges[step]
    local left = clusters[m.left]
    local right = clusters[m.right]
    local merged = {}
    for i = 1, #left do
      merged[#merged + 1] = left[i]
    end
    for i = 1, #right do
      merged[#merged + 1] = right[i]
    end
    table.sort(merged)

    clusters[m.left] = nil
    clusters[m.right] = nil
    clusters[m.id] = merged
  end

  local ordered = {}
  for _, members in pairs(clusters) do
    ordered[#ordered + 1] = { min_index = members[1], members = members }
  end
  table.sort(ordered, function(a, b)
    return a.min_index < b.min_index
  end)

  local assignments = {}
  local groups = {}
  for cluster_id = 1, #ordered do
    local members = ordered[cluster_id].members
    groups[cluster_id] = {}
    for i = 1, #members do
      local idx = members[i]
      assignments[idx] = cluster_id
      groups[cluster_id][#groups[cluster_id] + 1] = tree.labels[idx]
    end
  end
  return assignments, groups
end

-- ========================================================================================
--                                      Mean Shift
-- ========================================================================================
--[[
    PSEUDOCODE:
    1. For each point, repeat:
        a. Find all points within bandwidth h.
        b. Move the point to the mean of these neighbors.
        c. Stop if the point moves less than a small threshold or after max iterations.
    2. After shifting, merge points that are close together into clusters.
    3. Assign each original point to a cluster based on its shifted position.
]]
--https://www.tutorialspoint.com/machine_learning/machine_learning_mean_shift_clustering.htm
--https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html

-- quick bandwidth estimate
local function estimate_bandwidth(points)
  local n = #points
  -- fewer than 2 points just say bandwidth 10 idk kinda random
  if n < 2 then return 10 end

  --find and store nearest neighbor 
  local nn = {}
  for i = 1, n do
    local best = math.huge
    for j = 1, n do
      if i ~= j then
        local dx = points[i][1] - points[j][1]
        local dy = points[i][2] - points[j][2]
        local d = dx*dx + dy*dy
        if d < best then best = d end
      end
    end
    nn[#nn+1] = math.sqrt(best)
  end

  table.sort(nn)
  return nn[math.floor(#nn/2)]  -- median distance = estimate
end

function meanshift(points, h, max_iters)
  max_iters = max_iters or 20
  local n = #points
  if n == 0 then return nil, nil end

  local h2 = h*h

  --squared distance between two points
  local function dist2(a, b)
    local dx = a[1] - b[1]
    local dy = a[2] - b[2]
    return dx*dx + dy*dy
  end

  -- shift points toward local density peaks
  local modes = {}
  for i = 1, n do
    local y = { points[i][1], points[i][2] }
    for _ = 1, max_iters do
      local mx, my, c = 0, 0, 0
      -- find neighbors within bandwidth and compute mean 
      for j = 1, n do
        if dist2(y, points[j]) < h2 then
          mx = mx + points[j][1]
          my = my + points[j][2]
          c = c + 1
        end
      end
      --if no neighbors, stop
      if c == 0 then break end
      -- shifted position (mean of neighbors)
      local ny = { mx/c, my/c }
      -- stop if small
      if dist2(y, ny) < 1e-4 then break end
      y = ny
    end
    modes[i] = y
  end

  -- merge modes that are very close (cluster centers)
  local centers, assign = {}, {}
  local tol2 = h2 * 0.2

--check if near a closeby center, assign it or otherwise make a new center
  for i = 1, n do
    local m = modes[i]
    local done = false
    for j = 1, #centers do
      if dist2(m, centers[j]) < tol2 then
        assign[i] = j
        done = true
        break
      end
    end
    if not done then
      centers[#centers+1] = { m[1], m[2] }
      assign[i] = #centers
    end
  end

  return assign, centers, #centers
end

local function cmd_meanshift(model)
  local page = model:page()

local raw_points, sel_indices = collect_selected_points(page)

if #raw_points == 0 then
  model:warning("Select mark symbols to cluster!")
  return
end

local points = raw_points   
  -- bandwidth suggestion
  local suggested = estimate_bandwidth(points) * 3
  model:warning("Suggested bandwidth: " .. string.format("%.2f", suggested))

  -- ask user for bandwidth
  local s = model:getString("Enter bandwidth h (blank = suggested " .. string.format("%.2f", suggested) .. ")")

  if s == nil then
    return
  end

  -- if blank input, use suggested
  local h
  if s:match("^%s*$") then
    h = suggested
  else
    h = tonumber(s)
    if not h or h <= 0 then
      model:warning("Please enter a positive number for bandwidth.")
      return
    end
  end

  -- warning for very small h
  if h < suggested * 0.3 then
    model:warning("Warning: Bandwidth is very small — may produce many tiny clusters.")
  end

  -- run mean shift
  local assign, centers, k_est = meanshift(points, h, 20)
  if not assign then
    model:warning("Mean Shift failed.")
    return
  end

  -- give estimated clusters
  model:warning("Estimated clusters: " .. k_est)

  local colors = { "red","blue","green","orange","purple","cyan","brown","magenta","navy","gray" }

  local t = {
    label = "Mean Shift",
    pno = model.pno,
    vno = model.vno,
    original = page:clone(),
    final = page:clone(),
    undo = revertOriginal,
    redo = _G.revertFinal,
  }

  -- recolor points
  local idx = 1
  for i, obj, sel in t.final:objects() do
    if sel and obj:type() == "reference" then
      local c = assign[idx]
      local col = colors[((c - 1) % #colors) + 1]
      obj:set("stroke", col)
      obj:set("fill", col)
      idx = idx + 1
    end
  end

  model:register(t)

  -- draw cluster centers
  for j = 1, #centers do
    local col = colors[((j - 1) % #colors) + 1]
    local ref = ipe.Reference(
      model.attributes,
      "mark/cross(sx)",
      ipe.Vector(centers[j][1], centers[j][2])
    )
    ref:set("stroke", col)
    model:creation("Mean Shift Center", ref)
  end
  
  -- Convex hulls
  local clusterPoints = {}
  for j, p in ipairs(points) do
    local c = assign[j]
    if c then
      clusterPoints[c] = clusterPoints[c] or {}
      table.insert(clusterPoints[c], { x = p[1], y = p[2] })
    end
  end

  for j = 1, #centers do
    local color = colors[((j - 1) % #colors) + 1]
    if #clusterPoints[j] >= 3 then
      local hull = convex_hull(clusterPoints[j], model)
      if #hull >= 3 then
        local verts = {}
        for _, p in ipairs(hull) do
          table.insert(verts, ipe.Vector(p.x, p.y))
        end
        local shape = create_shape_from_vertices(verts, model)
        local path = ipe.Path(model.attributes, { shape })
        path:set("stroke", color)
        path:set("dashstyle", "dashed")
        model:creation("Cluster hull", path)
      end
    end
  end
end

-- methods that get passed to IPE
methods = {
  { label = "K-means", run = cmd_kmeans },
  { label = "K-means++", run = cmd_kmeanspp },
  { label = "K-medoids", run = cmd_kmedoids },
  { label = "DBSCAN", run = dbscan },
  { label = "HDBSCAN", run = hdbscan },
  { label = "Complete Linkage", run = function(model) cmd_hierarchical(model, CompleteLinkage) end },
  { label = "Single Linkage", run = function(model) cmd_hierarchical(model, SingleLinkage) end },
  { label = "Mean Shift", run = cmd_meanshift },
}


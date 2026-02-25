label = "Neighborhood Graphs"
revertOriginal = _G.revertOriginal
about = [[
    This ipelet makes the neighborhood graphs of different types for a set of sites
]]

function incorrect(title, model)
  model:warning(title)
end

function boundingTriangle(points,model)
	local minX, minY, maxX, maxY = math.huge, math.huge, -math.huge, -math.huge
	
	for _, p in ipairs(points) do
		if p.x < minX then minX = p.x end
		if p.x > maxX then maxX = p.x end
		if p.y < minY then minY = p.y end
		if p.y > maxY then maxY = p.y end
	end
	
	local difference = math.max(maxX-minX , maxY-minY)
	
	local p1 = { x = minX - 10*difference, y = minY - 10*difference}
	
	local p2 = { x = maxX + 10*difference, y = minY - 10*difference}
	
	local p3 = { x = minX + 10*difference, y = maxY + 10*difference}
	
	local boundingTriangle = { a = p1, b = p2, c = p3 }
	
	return boundingTriangle
end

function Circumcircle(A, B, C, D)
    -- Difference in the x values from points to new point
    local a = A.x - D.x
    local d = B.x - D.x
    local g = C.x - D.x
    -- Difference in the y values from points to new point
    local b = A.y - D.y
    local e = B.y - D.y
    local h = C.y - D.y
    -- Squared Euclidean distances
    local c = a*a + b*b
    local f = d*d + e*e
    local i = g*g + h*h
    -- Determinant: det = a(ei - fh) - b(di - fg) + c(dh - eg)
    local det = a * (e*i - f*h) - b * (d*i - f*g) + c * (d*h - e*g)
    return det
end

function get_unique_selected_points(points, model)
  local uniquePoints = {}
  for i = 1, #points do
    if not_in_table(uniquePoints, points[i]) then
      table.insert(uniquePoints, points[i])
    end
  end
  return uniquePoints
end


function create_shape_from_vertices_open(v, model)
  local shape = { type = "curve", closed = false }
  for i = 1, #v - 1 do 
    table.insert(shape, { type = "segment", v[i], v[i+1] })
  end
  return shape
end


function not_in_table(vectors, vector_comp)
    local flag = true
    for _, vertex in ipairs(vectors) do
        if vertex == vector_comp then
            flag = false
        end
    end
    return flag
end


function checkEdges(triangle, edge)
    local p1 = edge[1]
    local p2 = edge[2]

    local flag1 = (triangle.a == p1) or (triangle.b == p1) or (triangle.c == p1)
    local flag2 = (triangle.a == p2) or (triangle.b == p2) or (triangle.c == p2)

    if flag1 and flag2 then
        return true
    end
    return false
end

function edgeLength(p1,p2)
  -- We don't actually need to take the sqaure root.
    return (p2.x-p1.x)*(p2.x-p1.x)+(p2.y-p1.y)*(p2.y-p1.y)
end

function nearest_neighbor(model)
  run(model,"nearest")
end

function furthest_neighbor(model)
  run(model,"furthest")
end

function mutual_nearest_neighbor(model)
  run(model,"mutual")
end

function relative_neighbor(model)
  run(model,"relative")
end

function sphere_of_influence(model)
  run(model,"influence")
end

function gabriel(model)
  run(model,"gabriel")
end


function Urquhart(model)
  Ur(model)
end


function k_mutual_neighbor(model)
  local k = model:getString("Input k: ")
  if not k then return end
  k = tonumber(k)
  if k < 1 then
   incorrect("Value should be a positive integer",model)
   return
  end
  k = math.floor(tonumber(k))
  run(model,"k",k)
end

function kth_mutual_neighbor(model)
  local k = model:getString("Input k: ")
  if not k then return end
  k = tonumber(k)
  if k < 1 then
   incorrect("Value should be a positive integer",model)
   return
  end
  k = math.floor(tonumber(k))
  run(model,"kthmutual",k)
end

function symmetric_k_nearest_neighbor(model)
  local k = model:getString("Input k: ")
  if not k then return end
  k = tonumber(k)
  if k < 1 then
   incorrect("Value should be a positive integer",model)
   return
  end
  k = math.floor(tonumber(k))
  run(model,"symmetric_k",k)
end

function symmetric_nearest_neighbor(model)
  run(model,"symmetric_k", 1)
end

function symmetric_kth_nearest_neighbor(model)
  local k = model:getString("Input k: ")
  if not k then return end
  k = tonumber(k)
  if k < 1 then
   incorrect("Value should be a positive integer",model)
   return
  end
  k = math.floor(tonumber(k))
  run(model,"symmetric_kth",k)
end


function kth_nearest_neighbor(model)
  local k = model:getString("Input k: ")
  if not k then return end
  k = tonumber(k)
  if k < 1 then
   incorrect("Value should be a positive integer",model)
   return
  end
  k = math.floor(tonumber(k))
  run(model,"kth",k)
end

function k_nearest_neighbor(model)
  local k = model:getString("Input k: ")
  if not k then return end
  k = tonumber(k)
  if k < 1 then
   incorrect("Value should be a positive integer",model)
   return
  end
  k = math.floor(tonumber(k))
  run(model,"knearest",k)
end

function yao_graph(model)
  local k = model:getString("Number of sectors: ")
  if not k then return end
  k = tonumber(k)
  if k < 1 then
   incorrect("Value should be a positive integer",model)
   return
  end
  k = math.floor(tonumber(k))
  run(model, "yao", k)
end


function epsilon_neighbor(model)
  local epsi = model:getString("Input epsilon: ")
  if not epsi then return end
  epsi = tonumber(epsi)
  if epsi <= 0 then
   incorrect("Value should be positive", model)
   return
  end
  run(model,"epsilon", nil, epsi)
end

function findNearest(point,vertices)
  local minDist = math.huge
  local minPoints = {}
  for _, point2 in ipairs(vertices) do
    if point ~= point2 and edgeLength(point,point2)<minDist then
      minPoints = {}
      minDist = edgeLength(point,point2)
      table.insert(minPoints, point2)
    elseif point ~= point2 and edgeLength(point,point2)==minDist then
      table.insert(minPoints, point2)
    end
  end
  return minPoints
end

function inTable(point,points)
  for _, point2 in ipairs(points) do
      if point2 == point then return true end
  end
  return false
end

function run(model, version, k, epsi )
  local p = model:page()
  if not p:hasSelection() then incorrect("Please select at least 2 points", model) return end
  
  local points = {}
  -- Count the number of points
  local count = 0
  for _, obj, sel, _ in p:objects() do
    if sel then
      count = count + 1
      -- Make sure they are points
      if obj:type() ~= "reference" then
        incorrect("One or more selections are not points", model)
        return
      else
        -- Stick them in a table
        table.insert(points, obj:matrix() * obj:position())
      end
    end
  end
  
  if count < 2 then incorrect("Please select at least 2 points", model) return end
  
  if k and count <= k then
    incorrect("Value for k is larger than or equal to the amount of points", model) 
    return
  end
  
  local vertices = get_unique_selected_points(points,model)
  
  local edges = {}

  if version == "nearest" then
    for _, point1 in ipairs(vertices) do
      local near = findNearest(point1,vertices)
      for _, point2 in ipairs(near) do
        table.insert(edges, { a = point1, b = point2 })
      end
    end
  end
  
  if version == "furthest" then
    for _, point1 in ipairs(vertices) do
      local maxDist = 0
      local maxPoints = {}
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 and edgeLength(point1,point2)>maxDist then
          maxPoints = {point2}
          maxDist = edgeLength(point1,point2)
        elseif point1 ~= point2 and edgeLength(point1,point2)==maxDist then
          maxDist = edgeLength(point1,point2)
          table.insert(maxPoints,point2)
        end
      end
      for _, maxPoint in ipairs(maxPoints) do
        table.insert(edges, { a = point1, b = maxPoint })
      end
    end
  end
  
  if version == "kth" then
    for _, point1 in ipairs(vertices) do
      local Dist = {}
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 then
          table.insert(Dist, { point = point2, dist = edgeLength(point1, point2) })
        end
      end
      table.sort(Dist, function(a, b) return a.dist < b.dist end)
      if #Dist >= k then
        local kthDistance = Dist[k].dist
        for _, entry in ipairs(Dist) do 
          if entry.dist==kthDistance then
            table.insert(edges, { a = point1, b = entry.point})
          end
        end
      end
    end
  end
  
  if version == "knearest" then
    for _, point1 in ipairs(vertices) do
      local Dist = {}
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 then
          table.insert(Dist, { point = point2, dist = edgeLength(point1, point2) })
        end
      end
      table.sort(Dist, function(a, b) return a.dist < b.dist end)
      if #Dist >= k then
        local kthDistance = Dist[k].dist
        for _, entry in ipairs(Dist) do 
          if entry.dist<=kthDistance then
            table.insert(edges, { a = point1, b = entry.point})
          end
        end
      end
    end
  end
  
  if version == "yao" then
    local sectorAngle = 2 * math.pi / k
    local totalDist = 0
    for _, point1 in ipairs(vertices) do
      totalDist = totalDist + math.sqrt(edgeLength(point1,findNearest(point1,vertices)[1]))
    end
    local lineLen = totalDist / (2*#vertices)
    for _, point in ipairs(vertices) do
      for s = 0, k - 1 do
        local angle = s * sectorAngle
        local endX = point.x + lineLen * math.cos(angle)
        local endY = point.y + lineLen * math.sin(angle)
        
        local shape = { type = "curve", closed = false }
        table.insert(shape, { type = "segment", ipe.Vector(point.x, point.y), ipe.Vector(endX, endY) })
        local obj = ipe.Path(model.attributes, { shape })
        obj:set("dashstyle", "dashed")
        obj:set("stroke", "red")
        model:creation("Yao sector", obj)
      end
    end
    for _, point1 in ipairs(vertices) do
      local nearestInSector = {}
      local sectorsDistance = {}
      for s=1,k do
        sectorsDistance[s] = math.huge
        nearestInSector[s] = {}
      end
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 then
          local point2Angle = math.atan2(point2.y-point1.y,point2.x-point1.x)
          if point2Angle < 0 then
            point2Angle = point2Angle + 2*math.pi
          end
          local sector = math.floor(point2Angle / sectorAngle) + 1
          local dist = edgeLength(point1, point2) 
          
          if dist<sectorsDistance[sector] then
            sectorsDistance[sector]=dist
            nearestInSector[sector] = { point2 } 
          elseif dist == sectorsDistance[sector] then
            table.insert(nearestInSector[sector], point2) 
          end
        end
      end
      
      for s = 1, k do
        for _, point2 in ipairs(nearestInSector[s]) do
          table.insert(edges, { a = point1, b = point2 })
        end
      end
    end
  end
  
  if version == "k" then
    local kNearest = {}
    for _, point1 in ipairs(vertices) do
      local Dist = {}
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 then
          table.insert(Dist, { point = point2, dist = edgeLength(point1, point2) })
        end
      end
      table.sort(Dist, function(a, b) return a.dist < b.dist end)
      kNearest[point1] = {}
      local kDist = Dist[math.min(k, #Dist)].dist
      for _, entry in ipairs(Dist) do 
        if entry.dist<=kDist then
          table.insert(kNearest[point1], entry.point)
        end
      end
    end
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          if inTable(point2, kNearest[point1]) and inTable(point1, kNearest[point2]) then
            table.insert(edges, { a = point1, b = point2 })
          end
        end
      end
    end
  end
 
 
  if version == "kthmutual" then
    local kNearest = {}
    for _, point1 in ipairs(vertices) do
      local Dist = {}
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 then
          table.insert(Dist, { point = point2, dist = edgeLength(point1, point2) })
        end
      end
      table.sort(Dist, function(a, b) return a.dist < b.dist end)
      kNearest[point1] = {}
      local kDist = Dist[math.min(k, #Dist)].dist
      for _, entry in ipairs(Dist) do 
        if entry.dist==kDist then
          table.insert(kNearest[point1], entry.point)
        end
      end
    end
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          if inTable(point2, kNearest[point1]) and inTable(point1, kNearest[point2]) then
            table.insert(edges, { a = point1, b = point2 })
          end
        end
      end
    end
  end
 
 
  if version == "symmetric_k" then
      local kNearest = {}
    for _, point1 in ipairs(vertices) do
      local Dist = {}
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 then
          table.insert(Dist, { point = point2, dist = edgeLength(point1, point2) })
        end
      end
      table.sort(Dist, function(a, b) return a.dist < b.dist end)
      kNearest[point1] = {}
      local kDist = Dist[math.min(k, #Dist)].dist
      for _, entry in ipairs(Dist) do 
        if entry.dist<=kDist then
          table.insert(kNearest[point1], entry.point)
        end
      end
    end
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          if inTable(point2, kNearest[point1]) ~= inTable(point1, kNearest[point2]) then
            table.insert(edges, { a = point1, b = point2 })
          end
        end
      end
    end
  end
  
 if version == "symmetric_kth" then
    local kNearest = {}
    for _, point1 in ipairs(vertices) do
      local Dist = {}
      for _, point2 in ipairs(vertices) do
        if point1 ~= point2 then
          table.insert(Dist, { point = point2, dist = edgeLength(point1, point2) })
        end
      end
      table.sort(Dist, function(a, b) return a.dist < b.dist end)
      kNearest[point1] = {}
      local kDist = Dist[math.min(k, #Dist)].dist
      for _, entry in ipairs(Dist) do 
        if entry.dist == kDist then
          table.insert(kNearest[point1], entry.point)
        end
      end
    end
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          if inTable(point2, kNearest[point1]) ~= inTable(point1, kNearest[point2]) then
            table.insert(edges, { a = point1, b = point2 })
          end
        end
      end
    end
  end
  
  
  if version == "mutual" then
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          if inTable(point1, findNearest(point2, vertices)) and inTable(point2, findNearest(point1, vertices)) then
            table.insert(edges, { a = point1, b = point2 })
          end
        end
      end
    end
  end
  

  
  if version == "epsilon" then
    for _, point in ipairs(vertices) do
      local circle = { type = "ellipse", ipe.Matrix(epsi, 0, 0, epsi, point.x, point.y) }
      local obj = ipe.Path(model.attributes, { circle })
      obj:set("dashstyle", "dashed")
      obj:set("stroke", "red")
      model:creation("Epsilon circle", obj)
    end
    epsi = epsi*epsi
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then
          if edgeLength(point1, point2) <= epsi then
            table.insert(edges, { a = point1, b = point2})
          end
        end
      end
    end
  end
  
  if version == "relative" then
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          local dist12 = edgeLength(point1, point2) 
          local relative = true
          for _, point3 in ipairs(vertices) do
            if point3 ~= point1 and point3 ~= point2 then
              if edgeLength(point1, point3)  < dist12  and edgeLength(point2, point3)  < dist12 then
                relative = false
                break
              end
            end
          end
          if relative == true then
            table.insert(edges, { a = point1, b = point2})
          end
        end
      end
    end
  end
  
  if version == "gabriel" then
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          local dist12 = edgeLength(point1, point2) 
          local isGabriel = true
          for _, point3 in ipairs(vertices) do
            if point3 ~= point1 and point3 ~= point2 then
              local dist13 = edgeLength(point1,point3)
              local dist23 = edgeLength(point2,point3)
              if dist13 + dist23 < dist12 then
                isGabriel = false
                break
              end
            end
          end
          if isGabriel == true then
            table.insert(edges, { a = point1, b = point2})
          end
        end
      end
    end
  end
  
  if version == "influence" then
    local radii = {}
    for i, point1 in ipairs(vertices) do
      local near = findNearest(point1,vertices)
      radii[point1] = math.sqrt(edgeLength(point1, near[1]))
    end
    
      -- Draw circles
    for _, point in ipairs(vertices) do
      local r = radii[point]
      local circle = { type = "ellipse", ipe.Matrix(r, 0, 0, r, point.x, point.y) }
      local obj = ipe.Path(model.attributes, { circle })
      obj:set("dashstyle", "dashed")
      obj:set("stroke", "red")
      model:creation("Influence  region", obj)
    end
    
    for i, point1 in ipairs(vertices) do
      for j, point2 in ipairs(vertices) do
        if i < j then 
          local d = math.sqrt(edgeLength(point1, point2))
          if d <= radii[point1] + radii[point2 ] then
            table.insert(edges, { a = point1, b = point2})
          end
        end
      end
    end
  end
  
  for _, t in ipairs(edges) do
    local verts = {
        ipe.Vector(t.a.x, t.a.y),
        ipe.Vector(t.b.x, t.b.y),
    }
    
    local shape = create_shape_from_vertices_open(verts, model)
    local path = ipe.Path(model.attributes, { shape })
    model:creation("Nearest Neighbor Edge", path)
  end
end


function Ur(model)
  local p = model:page()
  if not p:hasSelection() then incorrect("Please select at least 3 points", model) return end
  
  local points = {}
  -- Count the number of points
  local count = 0
  for _, obj, sel, _ in p:objects() do
    if sel then
      count = count + 1
      -- Make sure they are points
      if obj:type() ~= "reference" then
        incorrect("One or more selections are not points", model)
        return
      else
        -- Stick them in a table
        table.insert(points, obj:matrix() * obj:position())
      end
    end
  end
  
  if count < 3 then incorrect("Please select at least 3 points", model) return end
    
  
  local vertices = get_unique_selected_points(points,model)

  if #vertices == 0 then
      model:warning("No reference points selected!")
      return
  end

  local superTriangle = boundingTriangle(vertices,model)
  local triangles = {superTriangle}

  for _, point in ipairs(vertices) do
      local badTriangles = {}
      for _, triangle in ipairs(triangles) do
          if Circumcircle(triangle.a, triangle.b, triangle.c, point)>0 then
              table.insert(badTriangles, triangle)
          end
      end
      local polygon = {}
      for _,triangle in ipairs(badTriangles) do
          local edges = {{ triangle.a, triangle.b },{ triangle.b, triangle.c },{ triangle.c, triangle.a }}
          for _, edge in ipairs(edges) do
            local shared = false
            for _,tri in ipairs(badTriangles) do 
              if tri ~= triangle and checkEdges(tri, edge) == true then 
                shared = true
                break
              end
            end
            if not shared then
              table.insert(polygon, edge)
            end
          end
      end
      for _,triangle in ipairs(badTriangles) do
          for i = #triangles, 1, -1 do
            if triangles[i] == triangle then
               table.remove(triangles, i)
            end
          end
        end
      for _, edge in ipairs(polygon) do
          local newTriangle = {a = edge[1], b = edge[2], c = point}
          table.insert(triangles, newTriangle)
      end
  end
  for i = #triangles, 1, -1 do
    local t = triangles[i]
    if t.a == superTriangle.a or t.a == superTriangle.b or t.a == superTriangle.c or
       t.b == superTriangle.a or t.b == superTriangle.b or t.b == superTriangle.c or
       t.c == superTriangle.a or t.c == superTriangle.b or t.c == superTriangle.c then
        table.remove(triangles, i)
    end
  end
    

  local badEdges = {}
  for _, t in ipairs(triangles) do
    local edges = {{ t.a, t.b }, { t.b, t.c }, { t.c, t.a }}
    local maxLength, maxIndex = 0, 1
    for i, edge in ipairs(edges) do
      local len = edgeLength(edge[1], edge[2])
      if len > maxLength then maxLength = len; maxIndex = i end
    end
    table.insert(badEdges, edges[maxIndex])
  end
  
  for _, t in ipairs(triangles) do
    for _, e in ipairs({ { t.a, t.b }, { t.b, t.c }, { t.c, t.a } }) do
      local bad = false
      for _, be in ipairs(badEdges) do
        if (be[1]==e[1] and be[2]==e[2]) or (be[1]==e[2] and be[2]==e[1]) then bad=true; break end
      end
      if not bad then
        local shape = create_shape_from_vertices_open({ ipe.Vector(e[1].x, e[1].y), ipe.Vector(e[2].x, e[2].y) }, model)
        model:creation("Urquhart Edge", ipe.Path(model.attributes, { shape }))
      end
    end
  end
end

methods = {
  { label = "Nearest neighbor", run = nearest_neighbor },
  { label = "kth nearest neighbor", run = kth_nearest_neighbor },
  { label = "k-nearest neighbor", run = k_nearest_neighbor },
  { label = "Mutual nearest neighbor", run = mutual_nearest_neighbor },
  { label = "kth mutual neighbor", run = kth_mutual_neighbor },
  { label = "k-mutual neighbor", run = k_mutual_neighbor },
  { label = "Asymmetric nearest neighbor", run = symmetric_nearest_neighbor },
  { label = "Asymmetric kth-nearest neighbor", run = symmetric_kth_nearest_neighbor },
  { label = "Asymmetric k-nearest neighbor", run = symmetric_k_nearest_neighbor },
  { label = "Furthest neighbor", run = furthest_neighbor },
  { label = "Relative neighbor", run = relative_neighbor },
  { label = "Gabriel graph", run = gabriel },
  { label = "Epsilon neighbor", run = epsilon_neighbor },
  { label = "Sphere of influence", run = sphere_of_influence },
  { label = "Urquhart", run = Urquhart},
  { label = "Yao graph", run = yao_graph },
}
label = "Points"

about = [[
Point Ipelets:
- Onion peeling
- Random point triangulation
- Smallest circle
]]


-- ---------------------------------------------------------------------------
-- Smallest circle
-- ---------------------------------------------------------------------------
do
    local get_dist
    local get_center
    local create_circle
    local is_in_circle
    local convex_hull_points
    local extreme_point
    local generate_smallest_circle
    local incorrect
    local run

    function get_dist(center,point)
        return math.sqrt((center.x - point.x)^2 + (center.y - point.y)^2)
    end

    function get_center(p1,p2,p3)
        local bi1 = ipe.Bisector(p1, p2)
        local bi2 = ipe.Bisector(p2, p3)
        local center = bi1:intersects(bi2)
        return center
    end

    function create_circle(model, center, radius)
        local shape =  { type="ellipse";
            ipe.Matrix(radius, 0, 0, radius, center.x, center.y) }
        model:creation("Smallest Circle",ipe.Path(model.attributes, { shape } ))
    end

    function is_in_circle(point, center, radius)
        if (point.x - center.x)^2 + (point.y - center.y)^2 <= radius^2 + 0.000001 then -- deal with floating points overflow
            return true
        else
            return false
        end
    end

    function convex_hull_points(points)

        function orient(p, q, r) return p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y) end
        function sortByX(a,b) return a.x < b.x end

        table.sort(points, sortByX)

        local upper = {}
        table.insert(upper, points[1])
        table.insert(upper, points[2])
        for i=3, #points do
            while #upper >= 2 and orient(points[i], upper[#upper], upper[#upper-1]) <= 0 do
                table.remove(upper, #upper)
            end
            table.insert(upper, points[i])
        end

    local lower = {}
        table.insert(lower, points[#points])
        table.insert(lower, points[#points-1])
        for i = #points-2, 1, -1 do
            while #lower >= 2 and orient(points[i], lower[#lower], lower[#lower-1]) <= 0 do
                table.remove(lower, #lower)
            end
            table.insert(lower, points[i])
        end

        table.remove(upper, 1)
        table.remove(upper, #upper)

        local S = {}
        for i=1, #lower do table.insert(S, lower[i]) end
        for i=1, #upper do table.insert(S, upper[i]) end

        return S

    end

    function extreme_point(points)
        local max = -1
        local p1 = nil 
        local p2 = nil  
        for i, point1 in ipairs(points) do 
            for j, point2 in ipairs(points) do
                if i ~= j then  -- Ensure we're not comparing the same points
                    local dist = math.sqrt((point1.x - point2.x)^2 + (point1.y - point2.y)^2)
                    if dist > max then
                        max = dist
                        p1 = point1
                        p2 = point2
                    end
                end
            end
        end
        return p1, p2  -- Return the two points with the maximum distance between them
    end

    function generate_smallest_circle(model, p1, p2, points)
        local temp_radius = math.maxinteger
        local temp_center = nil


        for i=1, #points do
            if not ((points[i].x == p1.x and points[i].y == p1.y) or (points[i].x == p2.x and points[i].y == p2.y)) then
                local center = get_center(p1,p2,points[i])
                local radius = get_dist(center,p1)
                local flag = true

                for j=1, #points do
                    if not is_in_circle(points[j], center,radius) then
                        flag = false
                        break
                    end
                end

                if flag and (radius <= temp_radius) then
                    temp_radius = radius
                    temp_center = center
                end


            end
        end

        
        local radius = get_dist(p1,p2)/2
        local center2 = ipe.Vector((p1.x + p2.x)/2, (p1.y + p2.y)/2)


        if temp_center == nil then
            create_circle(model, center2, radius)
            return
        elseif radius <= temp_radius then
            for j=1, #points do
                if not is_in_circle(points[j], center2,radius) then
                    create_circle(model, temp_center, temp_radius)
                    return
                end
            end
        end

        -- failed to generate the circle of 3 points
        create_circle(model, center2, radius)

    end

    function incorrect(title, model) model:warning(title) end

    function run(model)

        local p = model:page()

        if not p:hasSelection() then incorrect("Please select at least 2 points", model) return end

        local referenceObjects = {}
        local count = 0
        for _, obj, sel, _ in p:objects() do
            if sel then
            count = count + 1
                if obj:type() ~= "reference" then
                    incorrect("One or more selections are not points", model)
                    return
                else
                    table.insert(referenceObjects, obj:matrix() * obj:position())
                end
            end
        end
        
        if count < 2 then incorrect("Please select at least 2 points", model) return end

        if count == 2 then
            local p1 = referenceObjects[1]
            local p2 = referenceObjects[2]
            local center = ipe.Vector((p1.x + p2.x)/2, (p1.y + p2.y)/2)
            local radius = get_dist(p1,p2)/2
            create_circle(model, center, radius)
            return
        end

        local edge_points = convex_hull_points(referenceObjects)
        local extreme1, extreme2 = extreme_point(edge_points)


        generate_smallest_circle(model, extreme1, extreme2, edge_points)

    end

    smallest_circle_run = run


end


-- ---------------------------------------------------------------------------
-- Onion peeling
-- ---------------------------------------------------------------------------
do
    local incorrect
    local squared_distance
    local create_compare_function
    local orientation
    local convex_hull
    local create_shape_from_vertices
    local point_on_segment
    local create_segments_from_vertices
    local not_in_table
    local onion_peeling
    local run

    function incorrect(title, model) model:warning(title) end

    -- ========================================================================================================================
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
    -- ========================================================================================================================

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

    function not_in_table(t, v)
        for i=1, #t do
            if t[i] == v then return false end
        end
        return true
    end


    local creation_objects = {}
    function onion_peeling(points, model)

        if points == nil or #points <= 1 then return end
        if #points == 2 then
            table.insert(creation_objects, create_shape_from_vertices(points, model))
            return
        end
        local hull = convex_hull(points, model)
        local shape = create_shape_from_vertices(hull, model)
        table.insert(creation_objects, shape)

        
        local non_hull = {}
        local segments = create_segments_from_vertices(hull)
        for i=1, #points do
            if not_in_table(hull, points[i]) then
                local on_boundary = false
                for j=1, #segments do
                    if point_on_segment(points[i], segments[j]) then
                        on_boundary = true
                        break
                    end
                end
                if not on_boundary then table.insert(non_hull, points[i]) end
            end
        end
        
        onion_peeling(non_hull, model)
    end

    function run(model)
        local p = model:page()
        if not p:hasSelection() then incorrect("Please select at least 1 points", model) return end

        local referenceObjects = {}
        local count = 0
        for _, obj, sel, _ in p:objects() do
            if sel then
            count = count + 1
                if obj:type() ~= "reference" then
                    incorrect("One or more selections are not points", model)
                    return
                else
                    table.insert(referenceObjects, obj:matrix() * obj:position())
                end
            end
        end
        
        if count < 1 then incorrect("Please select at least 1 points", model) return end

        
        onion_peeling(referenceObjects, model)
        
        for i=1, #creation_objects do
            local shape = creation_objects[i]
            model:creation("onion peeling", ipe.Path(model.attributes, {shape}))
        end
        creation_objects = {};
    end

    onion_peeling_run = run

end


-- ---------------------------------------------------------------------------
-- Random point triangulation
-- ---------------------------------------------------------------------------

do
    local incorrect
    local create_shape_from_vertices
    local get_polygon_vertices
    local get_pt_and_polygon_selection
    local not_in_table
    local unique_points
    local triangleArea
    local newNode
    local newLinkedList
    local ll_add
    local ll_remove
    local ll_next_loop
    local ll_size
    local angleCCW
    local isConvex
    local insideTriangle
    local triangulate
    local reorient_ccw
    local orient
    local reverse_list
    local generateRandomPoints
    local randomPointInTriangle
    local run

    function incorrect(title, model)
    model:warning(title)
    end

    function create_shape_from_vertices(v, model)
    local shape = { type = "curve", closed = true }
    for i = 1, #v - 1 do 
        table.insert(shape, { type = "segment", v[i], v[i+1] })
    end
    table.insert(shape, { type = "segment", v[#v], v[1] })
    return shape
    end


    function get_polygon_vertices(obj, model)
    local shape = obj:shape()
    local m = obj:matrix()
    local vertices = {}
    local vertex = m * shape[1][1][1]
    table.insert(vertices, vertex)
    for i = 1, #shape[1] do
        vertex = m * shape[1][i][2]
        table.insert(vertices, vertex)
    end
    return vertices
    end

    function get_pt_and_polygon_selection(model)
    local p = model:page()
    if not p:hasSelection() then 
        incorrect("Please select a polygon", model)
        return 
    end

    local pathObject = nil
    local count = 0
    for _, obj, sel, _ in p:objects() do
        if sel then
        count = count + 1
        if obj:type() == "path" then
            pathObject = obj
        end
        end
    end

    if count ~= 1 then 
        incorrect("Please select one item.", model)
        return 
    end

    local vertices = get_polygon_vertices(pathObject, model)
    return vertices
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

    function unique_points(points, model)
    local uniquePoints = {}
    for i = 1, #points do
        if not_in_table(uniquePoints, points[i]) then
        table.insert(uniquePoints, points[i])
        end
    end
    return uniquePoints
    end


    function triangleArea(A, B, C)
    return 0.5 * math.abs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y))
    end



    function newNode(value)
    return { value = value, next = nil }
    end


    function newLinkedList()
    return { head = nil, length = 0 }
    end


    function ll_add(list, value)
    local node = newNode(value)
    if not list.head then
        list.head = node
    else
        local current = list.head
        while current.next do
        current = current.next
        end
        current.next = node
    end
    list.length = list.length + 1
    end


    function ll_remove(list, value)
    local current = list.head
    local previous = nil
    while current do
        if current.value == value then
        if previous then
            previous.next = current.next
        else
            list.head = current.next
        end
        list.length = list.length - 1
        return true
        end
        previous = current
        current = current.next
    end
    return false
    end


    function ll_next_loop(list, node)
    if not list.head then return nil end
    if node.next then
        return node.next
    else
        return list.head
    end
    end


    function ll_size(list)
    return list.length
    end


    function angleCCW(a, b)
    local dot = a.x * b.x + a.y * b.y
    local det = a.x * b.y - a.y * b.x
    local angle = math.atan2(det, dot)
    if angle < 0 then
        angle = 2 * math.pi + angle
    end
    return angle
    end


    function isConvex(vertex_prev, vertex, vertex_next)
    local a = { x = vertex_prev.x - vertex.x, y = vertex_prev.y - vertex.y }
    local b = { x = vertex_next.x - vertex.x, y = vertex_next.y - vertex.y }
    local internal_angle = angleCCW(b, a)
    return internal_angle <= math.pi
    end


    function insideTriangle(a, b, c, p)
    local v0 = { x = c.x - a.x, y = c.y - a.y }
    local v1 = { x = b.x - a.x, y = b.y - a.y }
    local v2 = { x = p.x - a.x, y = p.y - a.y }
    
    local dot00 = v0.x * v0.x + v0.y * v0.y
    local dot01 = v0.x * v1.x + v0.y * v1.y
    local dot02 = v0.x * v2.x + v0.y * v2.y
    local dot11 = v1.x * v1.x + v1.y * v1.y
    local dot12 = v1.x * v2.x + v1.y * v2.y
    
    local denom = dot00 * dot11 - dot01 * dot01
    if math.abs(denom) < 1e-20 then
        return true
    end
    local invDenom = 1.0 / denom
    local u = (dot11 * dot02 - dot01 * dot12) * invDenom
    local v = (dot00 * dot12 - dot01 * dot02) * invDenom
    
    return (u >= 0) and (v >= 0) and (u + v < 1)
    end


    function triangulate(vertices, model) 
    local triangles={}
    local n = #vertices
    local indices = {}  

    local vertlist = newLinkedList()
    for i = 1, n do
        ll_add(vertlist, i)
    end

    --local index_counter = 1
    local node = vertlist.head
    while ll_size(vertlist) > 2 do
        local i = node.value
        local j = ll_next_loop(vertlist, node).value
        local k = ll_next_loop(vertlist, ll_next_loop(vertlist, node)).value

        local vert_prev = vertices[i]
        local vert_current = vertices[j]
        local vert_next = vertices[k]
        
        local is_convex = isConvex(vert_prev, vert_current, vert_next)
        local is_ear = true
        if is_convex then
        local test_node = ll_next_loop(vertlist, ll_next_loop(vertlist, ll_next_loop(vertlist, node)))
        while test_node ~= node and is_ear do
            local vert_test = vertices[test_node.value]
            is_ear = not insideTriangle(vert_prev, vert_current, vert_next, vert_test)
            test_node = ll_next_loop(vertlist, test_node)
        end
        else
        is_ear = false
        end

        -- temp
        --[[
        if is_ear then
        indices[index_counter] = {vert_prev, vert_current, vert_next}
        index_counter = index_counter + 1
        ll_remove(vertlist, ll_next_loop(vertlist, node).value)
        end
        --]]
        
        if is_ear then
        --   local triangle = { vert_prev, vert_current, vert_next }
        local triangle_objs = { ipe.Vector(vert_prev.x, vert_prev.y), ipe.Vector(vert_current.x, vert_current.y), ipe.Vector(vert_next.x, vert_next.y)}
        local Tri = create_shape_from_vertices(triangle_objs, model)
        table.insert(triangles, triangle_objs)
        ll_remove(vertlist, j)
        end

        node = ll_next_loop(vertlist, node)
    end

    return triangles
    end

    --[=[
        Given: 
            
    {vertices}
    Return:
    {vertices ordered in clockwise fashion}
    ]=]
    function reorient_ccw(vertices)
        if orient(vertices[1], vertices[2], vertices[3]) < 0 then
            return reverse_list(vertices)
        end
        return vertices
    end

    function orient(p, q, r)
        local val = p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y)
        return val
    end

    function reverse_list(lst)
        local i = 1
        local j = #lst
        while i < j do
            local temp = lst[i]
            lst[i] = lst[j]
            lst[j] = temp
            i = i + 1
            j = j -1
        end

        return lst
    end

    -- Ok this guy should give us the area of a triangle
    function triangleArea(A, B, C)
        return 0.5 * math.abs(A.x*(B.y - C.y) + B.x*(C.y - A.y) + C.x*(A.y - B.y))
    end

    -- Function to generate random points inside the convex polygon
    function generateRandomPoints(numPoints, polygon)
    local triangles = triangulate(reorient_ccw(polygon), model)
    local points = {}
    local totalMass = 0
    
    for _, tri in ipairs(triangles) do
        totalMass = totalMass+triangleArea(tri[1],tri[2],tri[3])
    end
    
    
    -- This is 100% not the best way to do this, I know, but we just keep adding weight until it's greater than choice
    for i = 1, numPoints do
        local currentMass = 0
        local choice = math.random()*totalMass
        for _, tri in ipairs(triangles) do
        currentMass=currentMass+triangleArea(tri[1],tri[2],tri[3])
        if currentMass>=choice then 
            local randomPoint = randomPointInTriangle(tri[1], tri[2], tri[3])
            table.insert(points, randomPoint)
            break
        end
        end
    end

    return points
    end


    -- Function to generate a random point inside a triangle
    function randomPointInTriangle(p1, p2, p3)
    -- Get two random numbers between 0 and 1
    local u = math.random()
    local v = math.random()

    -- Make sure they aleways add up to less than or equal to 1
    if u + v > 1 then
        u = 1 - u
        v = 1 - v
    end

    -- This uses barycentric coordinates 
    local x = (1 - u - v) * p1.x + u * p2.x + v * p3.x
    local y = (1 - u - v) * p1.y + u * p2.y + v * p3.y
    -- returns the thing
    return {x = x, y = y}
    end

    function run(model)
        -- Set the seed for random generation (optional but useful for reproducibility)
        math.randomseed()
        local amount = model:getString("How many points to place?")

        local vertices = get_pt_and_polygon_selection(model)

        local points = generateRandomPoints(amount, unique_points(vertices))
    
        for i = 1, amount, 1 do
        -- Add the point
        local pointObj = ipe.Reference(
                                        model.attributes,
                                        model.attributes.markshape, 
                                        ipe.Vector(points[i].x, points[i].y)
                                )
                model:creation("A random Point", pointObj)
        end
    end

    random_point_triangulation_run = run

end


-- ---------------------------------------------------------------------------
-- Methods
-- ---------------------------------------------------------------------------

methods = {
  { label = "Smallest circle",            run = smallest_circle_run },
  { label = "Random point triangulation", run = random_point_triangulation_run },
  { label = "Onion peeling",              run = onion_peeling_run },
}
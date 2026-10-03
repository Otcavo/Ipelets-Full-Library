label = "Polygons"

about = [[
Polygon Ipelets:
- Minkowski sum
- Floating body
- Polygon intersection
- Polygon subtraction
- Polygon union
- Macbeath region
]]

-- ---------------------------------------------------------------------------
-- Polygon union, intersection and subtraction
--
-- Union and intersection work on any number of selected polygons, with the
-- same result as applying the operation one polygon at a time:
--   union         P1 ∪ P2 ∪ ... ∪ Pn
--   intersection  P1 ∩ P2 ∩ ... ∩ Pn
--   subtraction   primary selection minus the other one (exactly 2)
--
-- Method: split every edge where it crosses another polygon's edges, keep
-- the pieces that lie on the boundary of the result, and join the kept
-- pieces into closed loops. This handles convex and non-convex polygons,
-- disjoint pieces, holes and shared edges. Results with holes are drawn as
-- one path whose outer loops run counter-clockwise and holes clockwise.
-- ---------------------------------------------------------------------------

do
    local TOL = 1e-6      -- distance below which two points are the same

    local function cross(a, b) return a.x * b.y - a.y * b.x end
    local function dot(a, b) return a.x * b.x + a.y * b.y end
    local function len(a) return math.sqrt(dot(a, a)) end

    local function signed_area(v)
        local area = 0
        for i = 1, #v do
            local p, q = v[i], v[i % #v + 1]
            area = area + cross(p, q)
        end
        return area / 2
    end

    -- Vertices of one selected object in page coordinates, counter-clockwise.
    -- Returns nil and a message if the object is not a straight-edged polygon.
    local function polygon_vertices(obj)
        if obj:type() ~= "path" then
            return nil, "Every selected object must be a polygon"
        end
        local shape = obj:shape()
        if #shape ~= 1 or shape[1].type ~= "curve" then
            return nil, "Every selected object must be a single polygon (not a circle, ellipse or multi-part path)"
        end
        local m = obj:matrix()
        local raw = {}
        for i, seg in ipairs(shape[1]) do
            if seg.type ~= "segment" then
                return nil, "Polygons must have straight edges only"
            end
            if i == 1 then raw[#raw + 1] = m * seg[1] end
            raw[#raw + 1] = m * seg[2]
        end
        local v = {}
        for _, p in ipairs(raw) do
            if #v == 0 or len(p - v[#v]) > TOL then v[#v + 1] = p end
        end
        if #v > 1 and len(v[1] - v[#v]) <= TOL then table.remove(v) end
        if #v < 3 then
            return nil, "Every polygon needs at least three vertices"
        end
        if signed_area(v) < 0 then
            local r = {}
            for i = #v, 1, -1 do r[#r + 1] = v[i] end
            v = r
        end
        return v
    end

    -- Read the selection: every selected object must be a polygon and there
    -- must be at least two. The primary selection comes first.
    local function selected_polygons(model)
        local page = model:page()
        local primary = page:primarySelection()
        local polys = {}
        for i, obj, sel, _ in page:objects() do
            if sel then
                local v, err = polygon_vertices(obj)
                if not v then model:warning(err) return nil end
                if i == primary then
                    table.insert(polys, 1, v)
                else
                    polys[#polys + 1] = v
                end
            end
        end
        if #polys < 2 then
            model:warning("Please select at least 2 polygons")
            return nil
        end
        return polys
    end

    -- Shared pool of points, so that a crossing computed from both edges
    -- becomes one point and the kept pieces join up exactly.
    local function new_pool()
        local pts = {}
        local function id(p)
            for i, q in ipairs(pts) do
                if len(p - q) <= TOL then return i end
            end
            pts[#pts + 1] = p
            return #pts
        end
        return pts, id
    end

    -- Parameters t along segment a->b where it meets segment c->d.
    local function crossing_params(a, b, c, d)
        local r, s = b - a, d - c
        local ca = c - a
        local den = cross(r, s)
        local rl, sl = len(r), len(s)
        if math.abs(den) <= 1e-12 * rl * sl then
            -- parallel: only collinear overlaps matter
            if math.abs(cross(ca, r)) > TOL * rl then return {} end
            local rr = dot(r, r)
            return { dot(c - a, r) / rr, dot(d - a, r) / rr }
        end
        local t = cross(ca, s) / den
        local u = cross(ca, r) / den
        local et, eu = TOL / rl, TOL / sl
        if t >= -et and t <= 1 + et and u >= -eu and u <= 1 + eu then
            return { t }
        end
        return {}
    end

    -- Where a piece with midpoint m and direction dir sits relative to polygon v:
    -- "inside", "outside", "on_same" (along an edge, same direction) or
    -- "on_opposite" (along an edge, opposite direction).
    local function classify(m, dir, v)
        for i = 1, #v do
            local a, b = v[i], v[i % #v + 1]
            local ab = b - a
            local abl = len(ab)
            if math.abs(cross(m - a, ab)) <= TOL * abl then
                local t = dot(m - a, ab) / (abl * abl)
                if t >= -TOL / abl and t <= 1 + TOL / abl then
                    if dot(dir, ab) > 0 then return "on_same" end
                    return "on_opposite"
                end
            end
        end
        local inside = false
        local j = #v
        for i = 1, #v do
            local pi, pj = v[i], v[j]
            if ((pi.y > m.y) ~= (pj.y > m.y)) and
               (m.x < (pj.x - pi.x) * (m.y - pi.y) / (pj.y - pi.y) + pi.x) then
                inside = not inside
            end
            j = i
        end
        return inside and "inside" or "outside"
    end

    -- Split every edge of every polygon at its crossings with the others.
    -- Returns pieces {poly = k, a = id, b = id, mid = point, dir = vector}.
    local function split_edges(polys, id, pts)
        local pieces = {}
        for k, v in ipairs(polys) do
            for i = 1, #v do
                local a, b = v[i], v[i % #v + 1]
                local ts = { 0, 1 }
                for l, w in ipairs(polys) do
                    if l ~= k then
                        for j = 1, #w do
                            for _, t in ipairs(crossing_params(a, b, w[j], w[j % #w + 1])) do
                                if t > 0 and t < 1 then ts[#ts + 1] = t end
                            end
                        end
                    end
                end
                table.sort(ts)
                local prev = id(a)
                for n = 2, #ts do
                    local q = (n == #ts) and id(b) or id(a + ts[n] * (b - a))
                    if q ~= prev then
                        local pa, pb = pts[prev], pts[q]
                        pieces[#pieces + 1] = { poly = k, a = prev, b = q,
                                                mid = 0.5 * (pa + pb), dir = pb - pa }
                        prev = q
                    end
                end
            end
        end
        return pieces
    end

    -- Would this piece be on the boundary of the union of the polygons in
    -- `group` (indices)? Shared edges are kept once; abutting edges dropped.
    local function keep_for_union(piece, polys, group)
        for _, l in ipairs(group) do
            if l ~= piece.poly then
                local c = classify(piece.mid, piece.dir, polys[l])
                if c == "inside" or c == "on_opposite" then return false end
                if c == "on_same" and l < piece.poly then return false end
            end
        end
        return true
    end

    local function keep_for_intersection(piece, polys)
        for l, w in ipairs(polys) do
            if l ~= piece.poly then
                local c = classify(piece.mid, piece.dir, w)
                if c == "outside" or c == "on_opposite" then return false end
                if c == "on_same" and l < piece.poly then return false end
            end
        end
        return true
    end

    -- Join directed pieces {a, b} (point ids) into closed loops of points.
    local function link(kept, pts)
        local from = {}
        for n, s in ipairs(kept) do
            from[s[1]] = from[s[1]] or {}
            table.insert(from[s[1]], n)
        end
        local used, loops = {}, {}
        for n = 1, #kept do
            if not used[n] then
                local loop, cur = {}, n
                while cur do
                    used[cur] = true
                    loop[#loop + 1] = pts[kept[cur][1]]
                    local nxt
                    for _, k in ipairs(from[kept[cur][2]] or {}) do
                        if not used[k] then nxt = k break end
                    end
                    cur = nxt
                end
                -- drop vertices that only split a straight edge
                local clean = {}
                for i = 1, #loop do
                    local p, q, r = loop[(i - 2) % #loop + 1], loop[i], loop[i % #loop + 1]
                    local d1, d2 = q - p, r - q
                    if not (math.abs(cross(d1, d2)) <= TOL * len(d1) and dot(d1, d2) > 0) then
                        clean[#clean + 1] = q
                    end
                end
                if #clean >= 3 then loops[#loops + 1] = clean end
            end
        end
        return loops
    end

    local function draw(model, loops, label)
        local shape = {}
        for _, loop in ipairs(loops) do
            local curve = { type = "curve", closed = true }
            for i = 1, #loop - 1 do
                curve[#curve + 1] = { type = "segment", loop[i], loop[i + 1] }
            end
            shape[#shape + 1] = curve
        end
        model:creation(label, ipe.Path(model.attributes, shape))
    end

    local function all_indices(n)
        local t = {}
        for i = 1, n do t[i] = i end
        return t
    end

    function polygon_union_run(model)
        local polys = selected_polygons(model)
        if not polys then return end
        local pts, id = new_pool()
        local kept = {}
        local group = all_indices(#polys)
        for _, piece in ipairs(split_edges(polys, id, pts)) do
            if keep_for_union(piece, polys, group) then
                kept[#kept + 1] = { piece.a, piece.b }
            end
        end
        draw(model, link(kept, pts), "Create polygon union")
    end

    function polygon_intersect_run(model)
        local polys = selected_polygons(model)
        if not polys then return end
        local pts, id = new_pool()
        local kept = {}
        for _, piece in ipairs(split_edges(polys, id, pts)) do
            if keep_for_intersection(piece, polys) then
                kept[#kept + 1] = { piece.a, piece.b }
            end
        end
        local loops = link(kept, pts)
        if #loops == 0 then
            model:warning("The selected polygons have no region in common")
            return
        end
        draw(model, loops, "Create polygon intersection")
    end

    -- The primary selection minus the other selected polygon.
    -- Subtraction is not commutative, so it takes exactly two polygons.
    function polygon_sub_run(model)
        local polys = selected_polygons(model)
        if not polys then return end
        if #polys ~= 2 then
            model:warning("Please select exactly 2 polygons")
            return
        end
        local pts, id = new_pool()
        local cutters = {}
        for i = 2, #polys do cutters[#cutters + 1] = i end
        local kept = {}
        for _, piece in ipairs(split_edges(polys, id, pts)) do
            if piece.poly == 1 then
                -- keep the primary's boundary where no cutter covers it
                local keep = true
                for _, l in ipairs(cutters) do
                    local c = classify(piece.mid, piece.dir, polys[l])
                    if c == "inside" or c == "on_same" then keep = false break end
                end
                if keep then kept[#kept + 1] = { piece.a, piece.b } end
            else
                -- cutters' outline inside the primary, reversed
                if classify(piece.mid, piece.dir, polys[1]) == "inside"
                   and keep_for_union(piece, polys, cutters) then
                    kept[#kept + 1] = { piece.b, piece.a }
                end
            end
        end
        local loops = link(kept, pts)
        if #loops == 0 then
            model:warning("Nothing is left after the subtraction")
            return
        end
        draw(model, loops, "Create polygon subtraction")
    end

end

-- ---------------------------------------------------------------------------
-- Minkowski sum
-- ---------------------------------------------------------------------------
do
    local incorrect
    local print_vertices
    local print_table
    local print_vertex
    local debug_print
    local get_polygon_vertices
    local is_convex
    local copy_table
    local get_two_polygons_selection
    local minkowski
    local orient
    local convex_hull
    local create_shape_from_vertices
    local calculate_centroid
    local shift_polygon
    local center_minkowski_sum
    local not_in_table
    local unique_points
    local run

    function incorrect(title, model) model:warning(title) end

    function print_vertices(vertices, title, model)
        local msg = title ..  ": "
        for _, vertex in ipairs(vertices) do
            msg = msg .. ": " .. string.format("Vertex: (%f, %f), ", vertex.x, vertex.y)
        end
        model:warning(msg)
    end

    function print_table(t, title, model)
        -- Print lua table
        local msg = title ..  ": "
        for k, v in pairs(t) do
            msg = msg .. k .. " = " .. v .. ", "
        end
        model:warning(msg)
    end

    function print_vertex(v, title, model)
        local msg = title
        msg = msg .. ": " .. string.format("(%f, %f), ", v.x, v.y)
        model:warning(msg)
    end

    function print(x, title, model)
        local msg = title .. ": " .. x
        model:warning(msg)
    end

    function get_polygon_vertices(obj, model)

        local shape = obj:shape()
        local polygon = obj:matrix()

        local vertices = {}

            -- Apply transformation to the first vertex to handle translation
        local vertex = polygon * shape[1][1][1]
        table.insert(vertices, vertex)

            -- Apply transformation to the rest of the vertices to handle translation
        for i=1, #shape[1] do
            vertex = polygon * shape[1][i][2]
            table.insert(vertices, vertex)
        end

        return vertices
    end

    function is_convex(vertices)
        local _, convex_hull_vectors = convex_hull(vertices)
        return #convex_hull_vectors == #vertices
    end

    function copy_table(orig_table)
        local new_table = {}
        for i=1, #orig_table do new_table[i] = orig_table[i] end
        return new_table
    end

    function get_two_polygons_selection(model)
        local p = model:page()
        
        if not p:hasSelection() then incorrect("Please select 2 convex polygons", model) return end

        local pathObject1
        local pathObject2
        local count = 0
        local flag = true

        for _, obj, sel, _ in p:objects() do
            if sel then
                count = count + 1
                if obj:type() == "path" and flag then
                    pathObject1 = obj
                    flag = not flag
                else
                    if obj:type() == "path" then pathObject2 = obj end
                end
            end
        end

        if not pathObject1 or not pathObject2 then incorrect("Please select 2 convex polygons", model) return end

        local vertices1 = unique_points(get_polygon_vertices(pathObject1, model))
        local vertices2 = unique_points(get_polygon_vertices(pathObject2, model))

        local poly1_convex = is_convex(copy_table(vertices1))
        local poly2_convex = is_convex(copy_table(vertices2))

        if poly1_convex == false or poly2_convex == false then incorrect("Polygons must be convex", model) return end
        return vertices1, vertices2
    end

    --! MINKOWSKI SUM
    -- Compute the Minkowski Sum
    -- Uses the oriented cross product to ensure convexity and consistent vertex ordering
    function minkowski(P, Q, model)
        local result = {}
        for i=1, #P do for j=1, #Q do table.insert(result, P[i] + Q[j]) end end
        return result
    end

    function orient(p, q, r) return ((q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y-q.y)) < 0 end

    -- CONVEX HULL
    --[=[
    Given:
    - vertices: () -> {Vector}
    Return:
    - shape of the convex hull of points: () -> Shape
    --]=]
    function convex_hull(points)
        table.sort(points, function(a,b)
            if a.x < b.x then
                return true
            elseif a.x == b.x then
                return a.y < b.y
            else
                return false
            end
        end)
        if #points < 3 then return end
        local hull, left_most, p, q = {}, 1, 1, 0
        while true do
            table.insert(hull, points[p])
            q = (p % #points) + 1
            for i=1, #points do
                if orient(points[p], points[i], points[q]) then q = i end
            end
            p = q
            if p == left_most then break end
        end
        return create_shape_from_vertices(hull), hull
    end


    -- SHAPE CREATION
    function create_shape_from_vertices(v, model)
        local shape = {type="curve", closed=true;}
        for i=1, #v-1 do 
            table.insert(shape, {type="segment", v[i], v[i+1]})
        end
        table.insert(shape, {type="segment", v[#v], v[1]})
        return shape
    end

    --! CENTERING FUNCTIONS
    -- Function to calculate the centroid of a polygon
    function calculate_centroid(vertices)
        local sum_x, sum_y = 0, 0
        for _, v in ipairs(vertices) do
            sum_x = sum_x + v.x
            sum_y = sum_y + v.y
        end
        return ipe.Vector(sum_x / #vertices, sum_y / #vertices)
    end

    -- Function to shift the vertices of a polygon by a given vector
    function shift_polygon(vertices, shift_vector, model)
        local shifted_vertices = {}
        for _, v in ipairs(vertices) do
            table.insert(shifted_vertices, v + shift_vector)
        end
        return shifted_vertices
    end

    -- Function to center the Minkowski sum around the two input shapes
    function center_minkowski_sum(primary, secondary, minkowski_result, model)
        local centroid_primary = calculate_centroid(primary)
        local centroid_secondary = calculate_centroid(secondary)
        local centroid_minkowski = calculate_centroid(minkowski_result)

        -- Calculate the midpoint between the two input centroids
        local midpoint = ipe.Vector((centroid_primary.x + centroid_secondary.x) / 2, 
                                    (centroid_primary.y + centroid_secondary.y) / 2)

        -- Calculate the vector required to shift the Minkowski sum's centroid to the midpoint
        -- local shift_vector = ipe.Vector(midpoint.x - centroid_minkowski.x, 
        --                                 midpoint.y - centroid_minkowski.y)
        local shift_vector = midpoint - centroid_minkowski

        -- Shift the Minkowski sum to be centered around the midpoint
        return shift_polygon(minkowski_result, shift_vector, model)
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
        -- Check for duplicate points and remove them
        local uniquePoints = {}
        for i = 1, #points do
            if not_in_table(uniquePoints, points[i]) then table.insert(uniquePoints, points[i]) end
        end
        return uniquePoints
    end

    --! Run the Ipelet
    function run(model)
        if not get_two_polygons_selection(model) then return end
        local primary, secondary = get_two_polygons_selection(model)
        
        --! Compute the Minkowski sum of the two polygons and store resulting vertices
        local result_vertices = minkowski(primary, secondary, model)
        local centered_result_vertices = center_minkowski_sum(primary, secondary, result_vertices, model)

        --! Center the Minkowski sum around the two input shapes
        local result_shape_obj, _ = convex_hull(result_vertices)
        local centered_shape_obj, _ = convex_hull(centered_result_vertices)

        model:creation("Create Minkowski Sum", ipe.Path(model.attributes, { result_shape_obj }))
        model:creation("Create Centered Minkowski Sum", ipe.Path(model.attributes, { centered_shape_obj }))
    end

    minkowski_run = run

end

-- ---------------------------------------------------------------------------
-- Macbeath region
-- ---------------------------------------------------------------------------
do
    local get_polygon_segments
    local get_polygon_vertices
    local create_segments_from_vertices
    local get_polygon_vertices_and_segments
    local apply_transform
    local compute_macbeath_vertices
    local get_intersection_points
    local is_in_polygon
    local get_overlapping_points
    local create_shape_from_vertices
    local orient
    local sortByX
    local convex_hull
    local polygon_intersection
    local incorrect
    local is_convex
    local copy_table
    local get_pt_and_polygon_selection
    local not_in_table
    local unique_points
    local run

    function get_polygon_segments(obj, model)

        local shape = obj:shape()
        local translation = obj:matrix():translation()

        local segment_matrix = shape[1]

        local segments = {}
        for _, segment in ipairs(segment_matrix) do
            table.insert(segments, ipe.Segment(segment[1]+translation, segment[2]+translation))
        end
        
        table.insert(
            segments,
            ipe.Segment(segment_matrix[#segment_matrix][2]+translation, segment_matrix[1][1]+translation)
        )

        return segments
    end

    function get_polygon_vertices(obj, model)

        local shape = obj:shape()
        local polygon = obj:matrix()

        vertices = {}

        vertex = polygon * shape[1][1][1]
        table.insert(vertices, vertex)

        for i=1, #shape[1] do
            vertex = polygon * shape[1][i][2]
            table.insert(vertices, vertex)
        end

        return vertices
    end

    function create_segments_from_vertices(vertices)
        local segments = {}
        for i=1, #vertices-1 do
            table.insert( segments, ipe.Segment(vertices[i], vertices[i+1]) )
        end

        table.insert( segments, ipe.Segment(vertices[#vertices], vertices[1]) )
        return segments
    end

    function get_polygon_vertices_and_segments(obj, model)
        local vertices = get_polygon_vertices(obj)
        vertices = unique_points(vertices)
        local segments = create_segments_from_vertices(vertices)
        return vertices, segments
    end

    function apply_transform(v, point)
        return 2*point-v
    end

    function macbeath_vertices(orig_vertices, point)
        new_vertices = {}
        for i=1, #orig_vertices do 
            table.insert(new_vertices, apply_transform(orig_vertices[i], point))
        end
        return new_vertices
    end

    function get_intersection_points(s1,s2)
        local intersections = {}
        for i=1,#s2 do
            for j=1,#s1 do
                local intersection = s2[i]:intersects(s1[j])
                if intersection then
                    table.insert(intersections, intersection)
                end
            end
        end

        return intersections
    end

    function is_in_polygon(point, polygon)
        local x, y = point.x, point.y
        local j = #polygon
        local inside = false

        for i = 1, #polygon do
            local xi, yi = polygon[i].x, polygon[i].y
            local xj, yj = polygon[j].x, polygon[j].y

            if ((yi > y) ~= (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi) then
                inside = not inside
            end
            j = i
        end

        return inside
    end

    function get_overlapping_points(v1, v2)
        local overlap = {}
        for i=1, #v1 do
            if is_in_polygon(v1[i], v2) then
                table.insert(overlap, v1[i])
            end
        end
        return overlap
    end

    function create_shape_from_vertices(v, model)
        local shape = {type="curve", closed=true;}
        for i=1, #v-1 do 
            table.insert(shape, {type="segment", v[i], v[i+1]})
        end
        table.insert(shape, {type="segment", v[#v], v[1]})
        return shape
    end

    function orient(p, q, r)
        val = p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y)
        return val
    end

    function sortByX(a,b) return a.x < b.x end

    function convex_hull(points, model)
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

        return create_shape_from_vertices(S), S

    end

    function polygon_intersection(v1, s1, v2, s2, model)
        local intersections = get_intersection_points(s1, s2)
        local overlap1 = get_overlapping_points(v1, v2)
        local overlap2 = get_overlapping_points(v2, v1)

        local region = {}
        for i=1, #intersections do table.insert(region, intersections[i]) end
        for i=1, #overlap1 do table.insert(region, overlap1[i]) end
        for i=1, #overlap2 do table.insert(region, overlap2[i]) end

        local shape, _ = convex_hull(region)
        local region_obj = ipe.Path(model.attributes, { shape })
        region_obj:set("pathmode", "strokedfilled")

        return region_obj
    end

    function incorrect(title, model) model:warning(title) end

    function is_convex(vertices)
        local _, convex_hull_vectors = convex_hull(vertices)
        return #convex_hull_vectors == #vertices
    end

    function copy_table(orig_table)
        local new_table = {}
        for i=1, #orig_table do new_table[i] = orig_table[i] end
        return new_table
    end

    function get_pt_and_polygon_selection(model)

        local p = model:page()

        if not p:hasSelection() then incorrect("Please select a convex polygon and a point", model) return end

        local referenceObject
        local pathObject
        local count = 0

        for _, obj, sel, _ in p:objects() do
        if sel then
            count = count + 1
            if obj:type() == "path" then pathObject = obj end  -- assign pathObject
            if obj:type() == "reference" then referenceObject = obj end -- assign referenceObject
            end
        end

        if not referenceObject or not pathObject then incorrect("Please select a convex polygon and a point", model) return end

        local point = referenceObject:matrix() * referenceObject:position()  -- retrieve the point position (Vector)
        local vertices, segments = get_polygon_vertices_and_segments(pathObject, model)

        local poly1_convex = is_convex(copy_table(vertices))
        if poly1_convex == false then incorrect("Polygon must be convex", model) return end
        if not is_in_polygon(point, copy_table(vertices)) then incorrect("Point must be inside the polygon", model) return end

        return point, vertices, segments
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
        -- Check for duplicate points and remove them
        local uniquePoints = {}
        for i = 1, #points do
            if (not_in_table(uniquePoints, points[i])) then
                table.insert(uniquePoints, points[i])
            end
        end
        return uniquePoints
    end

    function run(model)

        if not get_pt_and_polygon_selection(model) then return end
        local point, original_vertices, segments = get_pt_and_polygon_selection(model)
        local macbeath_vertices = macbeath_vertices(original_vertices, point)
        local macbeath_shape = create_shape_from_vertices(macbeath_vertices)
        local macbeath_obj = ipe.Path(model.attributes, { macbeath_shape })
        local macbeath_segments = get_polygon_segments(macbeath_obj)

        local macbeath_region_obj = polygon_intersection(original_vertices, segments, macbeath_vertices, macbeath_segments, model)
        local obj2 =  ipe.Reference(model.attributes,model.attributes.markshape, point)

        model:creation("Macbeath Region", ipe.Group({macbeath_obj,macbeath_region_obj, obj2}))

    end

    macbeath_run = run

end

-- ---------------------------------------------------------------------------
-- Floating bodies
--
-- For a convex polygon and delta (as a percentage of its area), cuts off a
-- cap of exactly delta% of the area in each of 360 directions. Draws either
-- the cutting lines (as one group) or the polygon through their midpoints.
--
-- For each direction the cutting line is found by bisection on the exact
-- area of the polygon clipped to one side of the line, so every direction
-- gives a line: no special cases for edges parallel to the cut or vertices
-- at the same height.
-- ---------------------------------------------------------------------------

do
    local DIRECTIONS = 360
    local TOL = 1e-9

    local function cross(a, b) return a.x * b.y - a.y * b.x end
    local function dot(a, b) return a.x * b.x + a.y * b.y end

    local function area(v)
        local a = 0
        for i = 1, #v do a = a + cross(v[i], v[i % #v + 1]) end
        return math.abs(a) / 2
    end

    -- Vertices of the first selected closed polygon, without repeats and
    -- without vertices in the middle of a straight edge.
    local function selected_polygon(model)
        local p = model:page()
        for _, obj, sel, _ in p:objects() do
            if sel and obj:type() == "path" then
                local shape = obj:shape()
                if #shape ~= 1 or shape[1].type ~= "curve" then return nil end
                local m = obj:matrix()
                local raw = {}
                for i, seg in ipairs(shape[1]) do
                    if seg.type ~= "segment" then return nil end
                    if i == 1 then raw[#raw + 1] = m * seg[1] end
                    raw[#raw + 1] = m * seg[2]
                end
                local v = {}
                for _, q in ipairs(raw) do
                    if #v == 0 or math.abs(q.x - v[#v].x) + math.abs(q.y - v[#v].y) > TOL then
                        v[#v + 1] = q
                    end
                end
                if #v > 1 and math.abs(v[1].x - v[#v].x) + math.abs(v[1].y - v[#v].y) <= TOL then
                    table.remove(v)
                end
                local changed = true
                while changed and #v >= 3 do
                    changed = false
                    for i = 1, #v do
                        local a, b, c = v[(i - 2) % #v + 1], v[i], v[i % #v + 1]
                        local ab, bc = b - a, c - b
                        if math.abs(cross(ab, bc)) <= TOL * math.sqrt(dot(ab, ab) * dot(bc, bc)) then
                            table.remove(v, i)
                            changed = true
                            break
                        end
                    end
                end
                return v
            end
        end
        return nil
    end

    local function is_convex(v)
        local sign = 0
        for i = 1, #v do
            local o = cross(v[i % #v + 1] - v[i], v[(i + 1) % #v + 1] - v[i % #v + 1])
            if o ~= 0 then
                local s = o > 0 and 1 or -1
                if sign == 0 then sign = s elseif s ~= sign then return false end
            end
        end
        return true
    end

    -- The part of polygon v where dot(dir, x) <= s (Sutherland-Hodgman).
    local function clip(v, dir, s)
        local out = {}
        for i = 1, #v do
            local p, q = v[i], v[i % #v + 1]
            local dp, dq = dot(dir, p) - s, dot(dir, q) - s
            if dp <= 0 then out[#out + 1] = p end
            if (dp < 0 and dq > 0) or (dp > 0 and dq < 0) then
                out[#out + 1] = p + (dp / (dp - dq)) * (q - p)
            end
        end
        return out
    end

    -- Endpoints of the chord where the line dot(dir, x) = s crosses v.
    local function chord(v, dir, s)
        local pts = {}
        for i = 1, #v do
            local p, q = v[i], v[i % #v + 1]
            local dp, dq = dot(dir, p) - s, dot(dir, q) - s
            if dp == 0 then pts[#pts + 1] = p end
            if (dp < 0 and dq > 0) or (dp > 0 and dq < 0) then
                pts[#pts + 1] = p + (dp / (dp - dq)) * (q - p)
            end
        end
        if #pts < 2 then return nil end
        -- keep the two points farthest apart along the line
        local t = ipe.Vector(-dir.y, dir.x)
        local lo, hi = pts[1], pts[1]
        for _, p in ipairs(pts) do
            if dot(t, p) < dot(t, lo) then lo = p end
            if dot(t, p) > dot(t, hi) then hi = p end
        end
        return lo, hi
    end

    -- Offset s of the line that cuts off `target` area on the low side.
    local function cut_offset(v, dir, target)
        local lo, hi = math.huge, -math.huge
        for _, p in ipairs(v) do
            local d = dot(dir, p)
            if d < lo then lo = d end
            if d > hi then hi = d end
        end
        for _ = 1, 60 do
            local mid = (lo + hi) / 2
            if area(clip(v, dir, mid)) < target then lo = mid else hi = mid end
        end
        return (lo + hi) / 2
    end

    local function segment_path(model, a, b)
        local shape = { type = "curve", closed = false, { type = "segment", a, b } }
        return ipe.Path(model.attributes, { shape })
    end

    local function run(model)
        local v = selected_polygon(model)
        if not v or #v < 3 then
            model:warning("Please select a convex polygon")
            return
        end
        if not is_convex(v) then
            model:warning("The polygon must be convex")
            return
        end

        local delta = tonumber(model:getString("Enter delta value (1-99, where x means x% of the total area)"))
        if delta == nil or delta < 1 or delta > 99 then
            model:warning("Invalid delta input")
            return
        end
        local showType = model:getString("Would you like the halfspace lines (Default) or a polygon of midpoints (1)")

        local target = area(v) * delta / 100
        local lines, midpoints = {}, {}
        for i = 0, DIRECTIONS - 1 do
            local angle = 2 * math.pi * i / DIRECTIONS
            local dir = ipe.Vector(math.cos(angle), math.sin(angle))
            local a, b = chord(v, dir, cut_offset(v, dir, target))
            if a then
                lines[#lines + 1] = segment_path(model, a, b)
                midpoints[#midpoints + 1] = 0.5 * (a + b)
            end
        end

        if showType == "1" then
            local shape = { type = "curve", closed = true }
            for i = 1, #midpoints - 1 do
                shape[#shape + 1] = { type = "segment", midpoints[i], midpoints[i + 1] }
            end
            model:creation("Floating body midpoints", ipe.Path(model.attributes, { shape }))
        else
            model:creation("Floating body lines", ipe.Group(lines))
        end
    end

    floating_bodies_run = run

end


-- ---------------------------------------------------------------------------
-- Methods
-- ---------------------------------------------------------------------------

methods = {
  { label = "Polygon union",        run = polygon_union_run },
  { label = "Polygon subtraction",  run = polygon_sub_run },
  { label = "Polygon intersection", run = polygon_intersect_run },
  { label = "Minkowski sum",        run = minkowski_run },
  { label = "Macbeath region",      run = macbeath_run },
  { label = "Floating bodies",      run = floating_bodies_run },
}
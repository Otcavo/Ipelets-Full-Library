label = "Polygons"

about = [[
Polygon Ipelets:
- Minkowski sum
- Floating body
- Polygon intersection
- Polygon subtraction
- Polygon union
- Polar body
- Macbeath region
]]

-- ---------------------------------------------------------------------------
-- Polygon intersection 
-- ---------------------------------------------------------------------------

do
    local incorrect
    local get_selection_data
    local not_in_table
    local collect_vertices
    local is_in_polygon
    local distance_squared
    local make_clockwise
    local orientation
    local reorder_polygon
    local get_polygon_segments
    local get_lower_polygon
    local get_intersections
    local traverse_rest_points
    local get_polygon_vertices
    local create_shape_from_vertices
    local run


    function incorrect(model)
    model:warning("One or more selections are not polygons")
    end

    function get_selection_data(model)
        local page = model:page()
        local polygons = {}

        for i = 1, #page do
            local obj = page[i]
            if page:select(i) then
                if obj:type() ~= "path" then
                    incorrect(model)
                    return
                end
                table.insert(polygons, obj)
            end
        end

        if #polygons ~= 2 then
            model:warning("Please select 2 polygons")
        end

        return polygons[1], polygons[2]
    end

    -- COLLECT VERTICES
    function not_in_table(vertices, vertex_comp)
        local flag = true
        for _, vertex in ipairs(vertices) do
            if vertex == vertex_comp then
                flag = false
            end
        end
        return flag
    end

    function collect_vertices(obj)
        local vertices = {}

        local shape = obj:shape()

        for _, subpath in ipairs(shape) do
            for _, segment in ipairs(subpath) do
                if not_in_table(vertices, segment[1]) then
                    table.insert(vertices, segment[1])
                end
                if not_in_table(vertices, segment[2]) then
                    table.insert(vertices, segment[2])
                end
            end
        end

        return vertices
    end

    -- Check whether a vertex is in a polyon
    -- Adapted from the C Code on this website: https://alienryderflex.com/polygon/
    -- For now just sticking to the simplest code(though ineffecient)
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

    -- Closed set calculator
    -- Translated code to lua from https://www.geeksforgeeks.org/find-simple-closed-path-for-a-given-set-of-points/
    function distance_squared(point_a, point_b)
        return (point_a.x - point_b.x) * (point_a.x - point_b.x) + (point_a.y - point_b.y) * (point_a.y - point_b.y)
    end

    function make_clockwise(poly1, model)
        local reference = poly1[1]
        local should_reverse = false
        if orientation(poly1[1], poly1[2], poly1[3]) == 2 then
            should_reverse = true
        end
        if should_reverse then
            local i = 1
            local j = #poly1
            while i < j do
                local temp = poly1[i]
                poly1[i] = poly1[j]
                poly1[j] = temp
                i = i + 1
                j = j - 1
            end
        end
        return poly1
    end

    function orientation(point_a, point_b, point_c, model)
        local val = (point_b.y - point_a.y) * (point_c.x - point_b.x) - (point_b.x - point_a.x) * (point_c.y - point_b.y)
        if (val > 0) then
            return 1
        elseif (val < 0) then
            return 2
        else
            return 0
        end
    end

    function reorder_polygon(P)
        local pos = 1
        for i = 2, #P do
            if P[i].y < P[pos].y or (P[i].y == P[pos].y and P[i].x < P[pos].x) then
                pos = i
            end
        end
        local reordered = {}
        for i = pos, #P do
            table.insert(reordered, P[i])
        end
        for i = 1, pos - 1 do
            table.insert(reordered, P[i])
        end
        return reordered
    end

    function get_polygon_segments(obj, model)
        local shape = obj:shape()

        local segment_matrix = shape[1]

        local segments = {}
        for _, segment in ipairs(segment_matrix) do
            table.insert(segments, ipe.Segment(segment[1], segment[2]))
        end

        table.insert(
            segments,
            ipe.Segment(segment_matrix[#segment_matrix][2], segment_matrix[1][1])
        )

        return segments
    end

    function get_lower_polygon(poly1, poly2, model)
        if poly1[1].y < poly2[1].y then
            return poly1, poly2
        elseif poly1[1].y == poly2[1].y and poly1[1].x < poly2[1].x then
            return poly1, poly2
        else
            return poly2, poly1
        end
    end

    function get_intersections(poly1, poly2, model)
        local union_points = {}

        local visited_points = {}
        poly1 = make_clockwise(poly1, model)
        poly2 = make_clockwise(poly2, model)
        poly1 = reorder_polygon(poly1)
        poly2 = reorder_polygon(poly2)

        poly1, poly2 = get_lower_polygon(poly1, poly2, model)
        for i = 1, #poly1 do
            table.insert(union_points, poly1[i])
            local s1 = ipe.Segment(poly1[i], poly1[(i % #poly1) + 1])
            for j = 1, #poly2 do
                local s2 = ipe.Segment(poly2[j], poly2[(j % #poly2) + 1])
                local intersection = s1:intersects(s2)
                if intersection then
                    if not (is_in_polygon(poly2[(j % #poly2) + 1], poly1, model)) then
                        table.insert(union_points, intersection)
                        table.insert(visited_points, poly2[(j % #poly2) + 1])
                        union_points = traverse_rest_points(poly1, poly2, i + 1, (j % #poly2) + 1, union_points,
                            visited_points, model)

                        return union_points
                    end
                end
            end
        end
        return union_points
    end

    function traverse_rest_points(poly1, poly2, is, jn, union_points, visited_points, model)
        local ci = is
        local cj = jn
        local traversing_2 = true
        local has_intersected = false
        while ci ~= 1 do
            if traversing_2 then
                has_intersected = false
                local s1 = ipe.Segment(poly2[cj], poly2[(cj % #poly2) + 1])
                table.insert(union_points, poly2[cj])
                for i = 1, #poly1 do
                    -- print_vertex(poly1[i], "start", model)
                    -- print_vertex(poly1[(i % #poly1) + 1], "stop", model)
                    local s2 = ipe.Segment(poly1[i], poly1[(i % #poly1) + 1])
                    local intersection = s1:intersects(s2)
                    if intersection then
                        local x = math.floor(intersection.x)
                        local y = math.floor(intersection.y)
                        intersection = ipe.Vector(x, y)
                        if not_in_table(visited_points, poly1[(i % #poly1) + 1]) then
                            if not (is_in_polygon(poly1[(i % #poly1) + 1], poly2, model)) then
                                table.insert(union_points, intersection)
                                table.insert(visited_points, poly1[(i % #poly1) + 1])

                                traversing_2 = false
                                has_intersected = true
                                ci = (i % #poly1) + 1
                                break
                            end
                        else
                        end
                    end
                end
                if not has_intersected then
                    cj = cj % # poly2 + 1
                end
            else
                has_intersected = false
                table.insert(union_points, poly1[ci])
                local s1 = ipe.Segment(poly1[ci], poly1[(ci % #poly1) + 1])
                for j = 1, #poly2 do
                    local s2 = ipe.Segment(poly2[j], poly2[(j % #poly2) + 1])
                    local intersection = s1:intersects(s2)
                    if intersection then
                        local x = math.floor(intersection.x)
                        local y = math.floor(intersection.y)
                        intersection = ipe.Vector(x, y)
                        if not_in_table(visited_points, poly2[(j % #poly2) + 1]) then
                            if not (is_in_polygon(poly2[(j % #poly2) + 1], poly1, model)) then
                                table.insert(union_points, intersection)
                                table.insert(visited_points, poly2[(j % #poly2) + 1])
                                traversing_2 = true
                                cj = (j % #poly2) + 1
                                has_intersected = true
                                break
                            end
                        end
                    end
                end
                if not has_intersected then
                    ci = ci % #poly1 + 1
                end
            end
        end
        return union_points
    end

    function get_polygon_vertices(obj, model)
        local shape = obj:shape()
        local polygon = obj:matrix()

        vertices = {}

        vertex = polygon * shape[1][1][1]
        table.insert(vertices, vertex)

        for i = 1, #shape[1] do
            vertex = polygon * shape[1][i][2]
            table.insert(vertices, vertex)
        end

        return vertices
    end

    function create_shape_from_vertices(v, model)
        local shape = { type = "curve", closed = true, }
        for i = 1, #v - 1 do
            table.insert(shape, { type = "segment", v[i], v[i + 1] })
        end
        table.insert(shape, { type = "segment", v[#v], v[1] })
        return shape
    end

    function run(model)
        -- Obtain the first page of the Ipe document.
        -- Typically, work with thes objects (like polygons) on this page.
        --local page = model.doc[1]
        local page = model:page()
        local obj1, obj2 = get_selection_data(model)
        local obj1_vertices = get_polygon_vertices(obj1, model)
        local obj2_vertices = get_polygon_vertices(obj2, model)
        local shape = create_shape_from_vertices(get_intersections(obj1_vertices, obj2_vertices, model))
        local result_obj = ipe.Path(model.attributes, { shape }) -- Generate the original result shape
        model:creation("Create Polygon Union", result_obj)
    end

    polygon_union_run = run

end

-- ---------------------------------------------------------------------------
-- Polygon subtraction 
-- ---------------------------------------------------------------------------

do
    local incorrect
    local get_selection_data
    local not_in_table
    local collect_vertices
    local unique_points
    local create_segments_from_vertices
    local get_polygon_vertices_and_segments
    local get_intersection_points
    local distance_squared
    local is_between
    local cmp_pt_by_dist
    local get_midpoint
    local insert_intersection
    local orient
    local reverse_list
    local reorient_cw
    local is_in_polygon
    local find_cross_index
    local add_flags
    local process_vertices
    local get_outside_unused_point
    local perform_subtraction
    local incorrect_with_title
    local convex_hull
    local draw_shape
    local is_convex
    local run

    function incorrect(model)
        model:warning("One or more selections are not polygons")
    end

    function get_selection_data(model)
        local page = model:page()
        local primary_obj, secondary_obj
        local j = 0

        for i = 1, #page do
            if primary_obj ~= nil and secondary_obj ~= nil then
                break
            end

            local obj = page[i]
            if page:select(i) then
                if obj:type() ~= "path" then
                    incorrect(model)
                    return
                end

                if page:primarySelection() == i then
                    primary_obj = obj
                else
                    secondary_obj = obj
                end
                j = j + 1
            end
        end

        if j ~= 2 then
            model:warning("Please select 2 polygons")
            return nil, nil 
        end

        return primary_obj, secondary_obj
    end

    --[=[
        Given:
            - {vertices}
            - Vertex
        Return:
            True if vertex doesnt exist in vertices
            False otherweise 
    ]=]
    function not_in_table(vertices, vertex_comp)
        local flag = true
        for _, vertex in ipairs(vertices) do
            if vertex == vertex_comp then
                flag = false
            end
        end
        return flag
    end
    --[=[
        Given:
            - Ipelet Object
        Return:
            - List of vertices of Ipelet Object
    ]=]
    function collect_vertices(obj)
        local vertices = {}

        local shape = obj:shape()
        local m = obj:matrix()

        for _, subpath in ipairs(shape) do
            for _, segment in ipairs(subpath) do
                if not_in_table(vertices, m*segment[1]) then
                    table.insert(vertices, m*segment[1])
                end
                if not_in_table(vertices, m*segment[2]) then
                    table.insert(vertices, m*segment[2])
                end
            end
        end

        return vertices
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

    function create_segments_from_vertices(vertices)
        local segments = {}
        for i=1, #vertices-1 do
            table.insert( segments, ipe.Segment(vertices[i], vertices[i+1]) )
        end

        table.insert( segments, ipe.Segment(vertices[#vertices], vertices[1]) )
        return segments
    end

    function get_polygon_vertices_and_segments(obj, model)
        local vertices = collect_vertices(obj)
        vertices = unique_points(vertices)
        local segments = create_segments_from_vertices(vertices)
        return vertices, segments
    end

    --[=[
    Given:
    - vertices, segments of polygon A: () -> {Vector}, () -> {Segment} 
    - vertices, segments of polygon B: () -> {Vector}, () -> {Segment} 
    Return:
    - table of interection points: () -> {Vector}
    --]=]
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

    -- Closed set calculator
    -- Translated code to lua from https://www.geeksforgeeks.org/find-simple-closed-path-for-a-given-set-of-points/
    function distance_squared(point_a, point_b)
        return (point_a.x - point_b.x) * (point_a.x - point_b.x) + (point_a.y - point_b.y) * (point_a.y - point_b.y)
    end

    --[=[
        Taken from https://stackoverflow.com/questions/328107/how-can-you-determine-a-point-is-between-two-other-points-on-a-line-segment
        Given:
            - Vertex a, b, c
        Return:
            - true if vertex c exists between the line created between vertex a and b
            - false otherwise
    ]=]

    function is_between(a, b, c)
        
        local crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y)
        local epsilon = 1e-10

        if math.abs(crossproduct) > epsilon then
            return false
        end

        local dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y)
        if dotproduct < 0 then
            return false
        end

        local squaredlengthba = distance_squared(a,b)
        if dotproduct > squaredlengthba then
            return false
        end
        
        return true

    end

    --[=[
        Given:
            - Vertex c
        Return:
            - Comparator function comparing the distances of both vertices from c
    ]=]
    function cmp_pt_by_dist(c)
        return function (a, b)
            return distance_squared(a,c) < distance_squared(b, c)
        end
    end

    function get_midpoint(a, b)
        local x0 = 0 
        local y0 = 0
        x0 = (a.x + b.x)/2
        y0 = (a.y + b.y)/2
        return ipe.Vector(x0, y0)
    end

    --[=[
        Given:
            - Intersection points between two vertices 
            - {vertices from single polygon}
            - isPrimary check if it is primary polygon in order to insert extra points
                between intsection points on the same side
            - other_poly - for adding extra points, only add extra points if outside of other poly
        Return:
            {vertices in eithr CW ordering with intersection points in correct ordering}

    ]=]
    function insert_intersection(vertices, intersections, is_primary, other_poly)
        local new_v = {}
        for i=1, #vertices do
            local p_a = vertices[i]
            local p_b = vertices[(i % #vertices) + 1]
            table.insert(new_v, p_a)
            local seg_intersections = {}
            for _, i_point in ipairs(intersections) do 
                if is_between(p_a, p_b, i_point) then
                    table.insert(seg_intersections, i_point)
                end
            end

            -- In case we get multiple intersections on a segment, we need to figure out the ordering of inserting intersections by distance
            if #seg_intersections ~= 0 then
                table.sort(seg_intersections, cmp_pt_by_dist(p_a))
            end

            for i = 1, #seg_intersections do 
                table.insert(new_v, seg_intersections[i])
                if is_primary then
                    if i < #seg_intersections then
                        local new_pt = get_midpoint(seg_intersections[i], seg_intersections[i+1])
                        if not is_in_polygon(new_pt, other_poly) then
                            table.insert(new_v, new_pt)
                        end
                    end
                end
            end
        end

        return new_v
    end 

    -- val > 0 => CCW
    -- val < 0 => CW
    -- val == 0 => collinear
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

    --[=[
        Given: 
            - {vertices}
        Return:
            - {vertices ordered in clockwise fashion}
    ]=]
    function reorient_cw(vertices)
        if orient(vertices[1], vertices[2], vertices[3]) > 0 then
            return reverse_list(vertices)
        else
            return vertices
        end


    end

    --[=[
    Given:
    - point: () -> Vector
    - vertices of a polygon: () -> {Vector}
    Return:
    - returns true if point is inside the polygon, false otherwise
    - if the point is on the edge of a polygon, then false is returned
    - () -> Bool
    --]=]
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

    --[=[
        Given: 
            - Point p
            - List of vetices of given polygon 
        Return:
            - Index of vertex in polygon_v matching point p 
            - Otherwise -1 
    --]=]
    function find_cross_index(p, polygon_v)
        for i, v in ipairs(polygon_v) do 
            if v.x == p.x and v.y == p.y then 
                return i
            end
        end

        return -1
    end

    --[=[
        Reprocess the list of vertices with the following flags for each vertex:
        If the vertex has been processed
        The vertex reference
        If the vertex is outside of the alternate polygon
        The location of the point index in the other polygon, -1 otherwise
        Given:
            - {vertices to processed}
            - {vertices of alternate polygon}
        Return:
            - table{table{["vertex"], ["processed"], ["outside"], ["cross"]}}
    ]=]
    function add_flags(vertices, other_vertices)
        new_v = {}
        for _, v in ipairs(vertices) do
            table.insert(new_v, {["processed"] = false, ["vertex"] = v,
            ["outside"] = not is_in_polygon(v, other_vertices), ["cross"] = find_cross_index(v, other_vertices)})
        end

        return new_v
    end

    --[=[
        Processes the vertices by reorienting the vertices in counter clockwise, and adding in the intersection points
        Adds in the following flags to the vertices:
            - If it has been processed,
            - The index of the intersection point, in the other polygon vertex list 
            - If point is outside of the other polygon
            - 
        Return:
            - table{table{["vertex"], ["processed"], ["outside"], ["cross"]}}
    --]=]
    function process_vertices(primary_v, secondary_v, intersection)
        -- Reorient vertices to clockwise directions, then insert intersections in order 
        local primary_v = reorient_cw(primary_v)
        local secondary_v = reorient_cw(secondary_v)
        primary_v = insert_intersection(primary_v, intersection, true, secondary_v)
        secondary_v = insert_intersection(secondary_v, intersection, false, primary_v)
        -- Add flags have to be done after both vertices are processed with intersections
        local temp = add_flags(primary_v, secondary_v)
        secondary_v = add_flags(secondary_v, primary_v)
        primary_v = temp

        -- print_processed(primary_v)
        -- print("--------------------------")
        -- print("secondary")
        -- print_processed(secondary_v)

        return primary_v, secondary_v
    end

    function  get_outside_unused_point(vertices) 
        for i, v in ipairs(vertices) do 
            if v["outside"] and not v["processed"] and v["cross"] == -1 then
                return i 
            end
        end

        return -1
    end

    --[=[ 
        Perform polygon subtraction
        Algorithm pulled from: https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-SA-97135.pdf
        Given: 
            - Primary - primary vertices fully processed with following flags: vertex, processed, outside, crossed
            - Secondary - secondary vertices fully processed with following flags: vertex, processed, outside, crossed
            - table{table{["vertex"], ["processed"], ["outside"], ["cross"]}}
        Return:
            - list of resulting shape vertices from subtraction
            - table{list of table {list of vertices}}
    --]=] 
    function perform_subtraction(primary, secondary, model)
        local offset = 1
        local index = get_outside_unused_point(primary)
        local poly_operands = {[1] = primary, [2] = secondary}
        local curr_poly = 0
        local res_polys = {}
        local poly = {}

        while index ~= -1 do
            local v = poly_operands[(curr_poly % #poly_operands) + 1][index]
            if not v["processed"] then
                table.insert(poly, v["vertex"])
                v["processed"] = true
            
                if v["cross"] ~= -1 then
                    index = v["cross"]
                    offset = offset * -1
                    curr_poly = curr_poly + 1
                    poly_operands[(curr_poly % #poly_operands) + 1][index]["processed"] = true
                end
                index = index + offset
                if index == 0 then
                    index = #poly_operands[(curr_poly % #poly_operands) + 1]
                elseif index == (#poly_operands[(curr_poly % #poly_operands) + 1] + 1) then
                    index = 1
                end

            else
                table.insert(res_polys, poly)
                index = get_outside_unused_point(primary)
                poly = {}
                curr_poly = 0
                offset = 1
            end
        end 

        return res_polys

    end

    function incorrect_with_title(title, model) model:warning(title) end

    function convex_hull(points)

        function sortByX(a,b) return a.x < b.x end
        function orient_ch(p, q, r) return p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y) end
        
        local pts = {}
        for i=1, #points do pts[i] = points[i] end
        table.sort(pts, sortByX)
        
        local upper = {}
        table.insert(upper, pts[1])
        table.insert(upper, pts[2])
        for i=3, #pts do
            while #upper >= 2 and orient_ch(pts[i], upper[#upper], upper[#upper-1]) <= 0 do
                table.remove(upper, #upper)
            end
            table.insert(upper, pts[i])
        end

    local lower = {}
        table.insert(lower, pts[#pts])
        table.insert(lower, pts[#pts-1])
        for i = #pts-2, 1, -1 do
            while #lower >= 2 and orient_ch(pts[i], lower[#lower], lower[#lower-1]) <= 0 do
                table.remove(lower, #lower)
            end
            table.insert(lower, pts[i])
        end

        table.remove(upper, 1)
        table.remove(upper, #upper)
        
        local S = {}
        for i=1, #lower do table.insert(S, lower[i]) end
        for i=1, #upper do table.insert(S, upper[i]) end

        return S

    end

    --[=[
        Function for drawing a new instance of the same shape wihtout overriding 
        orignal shape's properties
        Given: table{vertices}
        Return: Shape object using given vertices
    --]=]
    function draw_shape(vertices, model) 
        local result_shape = { type = "curve", closed = true, }
        
        for i = 1, #vertices - 1 do
            table.insert(result_shape, { type = "segment", vertices[i], vertices[i + 1]})
        end

        local result_obj = ipe.Path(model.attributes, { result_shape })
        result_obj:set("pathmode", "stroked")
        result_obj:set("stroke", "red")
        return result_obj
    end

    function is_convex(v)
        local convex_hull_vectors = convex_hull(v)
        return #convex_hull_vectors == #v
    end

    function run(model)
        local page = model:page()
        local primary_obj, secondary_obj = get_selection_data(model)
        if primary_obj == nil or secondary_obj == nil then
            return
        end
        
        local p_v, p_s = get_polygon_vertices_and_segments(primary_obj, model)
        local s_v, s_s = get_polygon_vertices_and_segments(secondary_obj, model)

        if not is_convex(p_v) or not is_convex(s_v) then
            incorrect_with_title("Polygons are not convex. Polygon subtraction might not work as expected.", model)
        end

        local intersections = get_intersection_points(p_s, s_s)

        if #intersections == 0 then
            -- If completely enclosed, draw shape with hole
            if is_in_polygon(s_v[1], p_v ) then
                res_obj_lst = {draw_shape(p_v, model),draw_shape(s_v, model)}
                model:creation("Create polygon subtraction", ipe.Group(res_obj_lst))
            end

            return 
        end

        p_v, s_v = process_vertices(p_v, s_v, intersections)
        local res_polys = perform_subtraction(p_v, s_v, model)
        local objs = {}
        for _, s in ipairs(res_polys) do 
            local result_shape = { type = "curve", closed = true, }
            if #s >= 3 then 
                for i = 1, #s - 1 do
                    table.insert(result_shape, { type = "segment", s[i], s[i + 1]})
                end

                local result_obj = ipe.Path(model.attributes, { result_shape })
                table.insert(objs, result_obj)
            end
        end

        model:creation("Create polygon subtraction", ipe.Group(objs))

    end

    polygon_sub_run = run

end

-- ---------------------------------------------------------------------------
-- Polygon intersection 
-- ---------------------------------------------------------------------------
do
    local incorrect
    local get_polygon_vertices
    local create_segments_from_vertices
    local not_in_table
    local unique_points
    local get_polygon_vertices_and_segments
    local is_convex
    local copy_table
    local get_two_polygons_selection
    local get_intersection_points
    local is_in_polygon
    local get_overlapping_points
    local create_shape_from_vertices
    local orient
    local convex_hull
    local polygon_intersection
    local run

    function incorrect(title, model) model:warning(title) end

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

    function get_polygon_vertices_and_segments(obj, model)
        local vertices = get_polygon_vertices(obj)
        vertices = unique_points(vertices)
        local segments = create_segments_from_vertices(vertices)
        return vertices, segments
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

        local vertices1, segments1 = get_polygon_vertices_and_segments(pathObject1, model)
        local vertices2, segments2 = get_polygon_vertices_and_segments(pathObject2, model)
        
        local poly1_convex = is_convex(copy_table(vertices1))
        local poly2_convex = is_convex(copy_table(vertices2))
        if poly1_convex == false or poly2_convex == false then incorrect("Polygons must be convex", model) return end

        return vertices1, segments1, vertices2, segments2
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

    function get_overlapping_points(v1, v2, model)
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

    function convex_hull(points, model)

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

        return create_shape_from_vertices(S), S

    end

    function polygon_intersection(v1, s1, v2, s2, model)
        local intersections = get_intersection_points(s1, s2)
        if #intersections == 0 then
            return nil
        end
        local overlap1 = get_overlapping_points(v1, v2, model)
        local overlap2 = get_overlapping_points(v2, v1, model)

        local region = {}
        for i=1, #intersections do table.insert(region, intersections[i]) end
        for i=1, #overlap1 do table.insert(region, overlap1[i]) end
        for i=1, #overlap2 do table.insert(region, overlap2[i]) end

        local shape, _ = convex_hull(region)
        local region_obj = ipe.Path(model.attributes, { shape })
        return region_obj
    end

    --! Run the Ipelet
    function run(model)
        if not get_two_polygons_selection(model) then return end
        local v1,s1,v2,s2 = get_two_polygons_selection(model)
        local obj =	polygon_intersection(v1,s1,v2,s2, model)
        if obj == nil then
            return
        end
        model:creation("polygon intersection", obj)
    end

    polygon_intersect_run = run

end

-- ---------------------------------------------------------------------------
-- Polar body
-- ---------------------------------------------------------------------------

do
    local incorrect
    local convex_hull
    local is_convex
    local copy_table
    local get_original_vertices
    local vertex_dual
    local dual_transform
    local intersect
    local get_intersection_points
    local create_polar_body
    local shift_to_origin
    local shift_back
    local not_in_table
    local unique_points
    local run

    function incorrect(title, model) model:warning(title) end

    function convex_hull(points, model)

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

    function is_convex(vertices)
        local convex_hull_vectors = convex_hull(vertices)
        return #convex_hull_vectors == #vertices
    end

    function copy_table(orig_table)
        local new_table = {}
        for i=1, #orig_table do new_table[i] = orig_table[i] end
        return new_table
    end

    function get_original_vertices(model)
        local p = model:page()

        if not p:hasSelection() then incorrect("Please select a convex polygon", model) return end

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

        if not pathObject then incorrect("Please select a convex polygon", model) return end

        local shape = pathObject:shape()
        local polygon = pathObject:matrix()

        local orig_vertices = {}

        local vertex = polygon * shape[1][1][1]
        table.insert(orig_vertices, vertex)

        for i=1, #shape[1] do
            vertex = polygon * shape[1][i][2]
            table.insert(orig_vertices, vertex)
        end

        orig_vertices = unique_points(orig_vertices)

        if not is_convex(copy_table(orig_vertices)) then incorrect("Selected polygon is not convex", model) return end

        return orig_vertices, referenceObject
    end

    function vertex_dual(v)
        local a = v.x
        local b = v.y
        local p1
        local p2
        if a == 0 then
            p1 = ipe.Vector(0, 1/b)
            p2 = ipe.Vector(1, 1/b)
        elseif b == 0 then
            p1 = ipe.Vector(0, 1/a)
            p2 = ipe.Vector(1, 1/a)
        else
            p1 = ipe.Vector(1/a, 0)
            p2 = ipe.Vector(0, 1/b)
        end

        return ipe.LineThrough(p1, p2)
    end

    function dual_transform(v,model)
        lines = {}
        for i=1, #v do table.insert(lines, vertex_dual(v[i])) end
        return lines
    end

    function intersect(l1,l2, model)
        return l1:intersects(l2)
    end


    function get_intersection_points(l,model)
        polar_vertices = {}
        for i=1, #l-1 do table.insert(polar_vertices, intersect(l[i], l[i+1])) end
        table.insert(polar_vertices, intersect(l[#l], l[1], model))
        return polar_vertices
    end

    function create_polar_body(v, model)
        local shape = {type="curve", closed=true;}
        for i=1, #v-1 do table.insert(shape, {type="segment", v[i], v[i+1]}) end
        table.insert(shape, {type="segment", v[#v], v[1]})
        local obj = ipe.Path(model.attributes, { shape })
        model:creation("Polar Dual", obj)
    end

    function shift_to_origin(v)

        -- centroid calculation
        local x = 0
        local y = 0
        for _, vertex in ipairs(v) do
            x = x + vertex.x
            y = y + vertex.y
        end

        x = x / #v
        y = y / #v
        
        local shifted_vertices = {}
        for _, vertex in ipairs(v) do
            table.insert(shifted_vertices, ipe.Vector(vertex.x-x, vertex.y-y))
        end
        
        return shifted_vertices, x, y
    end

    -- Centers the polar body within the original polygon
    -- also applies a scaling factor to make the body more visible
    function shift_back(v, x, y)
        
        local shifted_vertices = {}
        for _, vertex in ipairs(v) do
            table.insert(shifted_vertices, ipe.Vector((2048*vertex.x)+x, (2048*vertex.y)+y))
        end
        
        return shifted_vertices
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

    --! Run the Ipelet
    function run(model)
        
        local orig_vertices, origin_obj = get_original_vertices(model)
        if not orig_vertices then return end

        local origin
        if origin_obj then
            origin = origin_obj:matrix() * origin_obj:position()
        end

        local orig_vertices, offset_x, offset_y = shift_to_origin(orig_vertices)
        local lines = dual_transform(orig_vertices, model)
        local polar_vertices = get_intersection_points(lines)

        if origin then
            polar_vertices = shift_back(polar_vertices, origin.x, origin.y)
            create_polar_body(polar_vertices, model)
            local obj =  ipe.Reference(model.attributes,model.attributes.markshape, origin)
            model:creation("Polar Dual Origin", obj)
        else
            polar_vertices = shift_back(polar_vertices, offset_x, offset_y)
            create_polar_body(polar_vertices, model)
            local obj = ipe.Reference(model.attributes,model.attributes.markshape, ipe.Vector(offset_x, offset_y))
            model:creation("Polar Dual Origin", obj)
        end
    end

    polar_body_run = run

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
-- ---------------------------------------------------------------------------

do
    local get_pts_and_segment_selection
    local reorderTable
    local isConvex
    local calculate_polygon_area
    local compare
    local get_closest
    local get_intersection_points
    local euclidean_distance
    local get_segments_sorted
    local get_magnitude
    local get_slope
    local get_angle
    local get_angle_between_segments
    local get_area
    local get_line
    local get_inverse
    local get_perp_line
    local get_dist_between_lines
    local get_dist_for_points
    local calc_area
    local reverseTable
    local find_h
    local removeDuplicates
    local print_vector
    local run
    function get_pts_and_segment_selection(model)
        local p = model:page()
        if not p:hasSelection() then
            return
        end

        local tab = {}

        for _, obj, sel, _ in p:objects() do
            local transform = obj:matrix()
            if sel then
                if obj:type() == "path" then 
                    local shape = obj:shape()
                    for _, subpath in ipairs(shape) do
                        if subpath.type == "curve" then
                            local max = 0
                            for j, vertices in ipairs(subpath) do
                                table.insert(tab, transform*vertices[1])
                                max = max + 1
                            end
                            table.insert(tab, transform*subpath[max][2])
                        end
                    end
                    -- Only finds the first selected path object
                    break
                end
            end
        end

        return tab
    end

    -- Reorders a table of sequential points starting at a provided value
    -- Courtesy of GPT
    function reorderTable(t, startElement)
        -- Find the index of the start element
        local startIndex = nil
        for i, v in ipairs(t) do
            if v == startElement then
                startIndex = i
                break
            end
        end

        -- If the element is not found, return the original table
        if startIndex == nil then
            print("Element not found in the table")
            return t
        end

        -- Create a new table to store the reordered result
        local reordered = {}

        -- Add the elements from the start index to the end
        for i = startIndex, #t do
            table.insert(reordered, t[i])
        end

        -- Add the elements from the beginning to the start index - 1
        for i = 1, startIndex - 1 do
            table.insert(reordered, t[i])
        end

        return reordered
    end

    -- Givena table of adjacent points of a polygon, returns true if it is convex
    function isConvex(vertices)
        function orient(p, q, r) return p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y) end
        local side = nil
        local i = 1
        while i < #vertices do
            local temp = orient(vertices[i], vertices[(i % #vertices)+1], vertices[((i+1)%#vertices)+1])
            if side == nil then
                if temp > 0 then
                    side = 1
                end
                if temp < 0 then
                    side = -1
                end
            else
                if side*temp < 0 then
                    return false
                end
            end
            i = i + 1

        end
        return true
    end

    -- Calculate the area of a polygon given its vertices (Shoelace algorithm)
    -- Code generated via Copilot
    function calculate_polygon_area(vertices)
        local area = 0
        local num_vertices = #vertices

        for i = 1, num_vertices do
            local j = (i % num_vertices) + 1
            area = area + (vertices[i].x * vertices[j].y) - (vertices[j].x * vertices[i].y)
        end

        return math.abs(area) / 2
    end

    -- Function to compare points based on dot product, then x value, then y
    -- (Allows a certain chain to always be on a certain side of another)
    function compare(a, b)
        if a[1] < b[1] then
            return true
        end

        if a[1] > b[1] then
            return false
        end

        if a[2].x < b[2].x then
            return true
        end

        if a[2].x > b[2].x then
            return false
        end

        if a[2].y < b[2].y then
            return true
        end

        return false
    end

    -- Orders points based on dot product
    function get_closest(points, dir, model)
        -- Stores the dot product to sort by
        local temp = {}
        for _, pt in ipairs(points) do
            table.insert(temp, {dir^pt, pt})
        end

        --table.sort(temp, compare) Just trying it, shouldn't matter though
        table.sort(temp, function(a,b) return a[1]<b[1] end)

        -- Convert it back into a table of vertices
        local newPoints = {}
        for _, pt in ipairs(temp) do
            table.insert(newPoints, pt[2])
        end

        return newPoints
    end


    --[=[
    Given:
    - vertices, segments of polygon A: () -> {Vector}, () -> {Segment} 
    - vertices, segments of polygon B: () -> {Vector}, () -> {Segment} 
    Return:
    - table of interection points: () -> {Vector}
    --]=]
    function get_intersection_points(s1,s2, model)
        local intersections = {}
        
        local end_point = s1[#s1]
        for i=1,#s1 do
            local intersection = ipe.Segment(s1[i], end_point):intersects(s2)
            if intersection then
                table.insert(intersections, intersection)
            end
            end_point = s1[i]
        end
        return intersections
    end



    -- euclidean distance between two points
    function euclidean_distance(p1, p2)
        return math.sqrt(math.pow((p1.x-p2.x), 2) + math.pow((p1.y-p2.y), 2))
    end


    -- Use the first and last point in sorted to construct the chains
    function get_segments_sorted(coordinates, sorted, model)
        local chain1 = {}
        local chain2 = {}
        local first = sorted[1]
        local last = sorted[#sorted]
        local index

        for i, x in ipairs(coordinates) do
            if first == x then
                index = i
                break
            end
        end

        if index == nil then
            model:warning("Issue finding segment in get_segments_sorted: ")
        end

        local i = index - 1
        local chainFlag = true

        local cur

        local flag = true
        while flag or i ~= index-1 do
            flag = false
            cur = coordinates[i+1]
            if cur == last then
                chainFlag = false
            end

            if chainFlag then
                table.insert(chain1, cur)
            else
                table.insert(chain2, 1, cur)
            end

            i = (i + 1) % #coordinates
        end

        -- In order to make sure they have access to the segments, we add the first and last node to both
        -- It is easy enough to check to see if it is the first or the last node. 
        -- Requires a check for othogonality
        table.insert(chain1, last)
        table.insert(chain2, 1, first)
        return {chain1, chain2}

    end

    -- Used for get_angle
    function get_magnitude(p1)
        return math.sqrt(math.pow(p1.x, 2) + math.pow(p1.y, 2))
    end

    -- get the difference between two vertices
    function get_slope(a,b)
        return ipe.Vector(a.x-b.x, a.y-b.y)
    end

    -- Given two angle vectors, get the angle between them (Used in get_angle_between_segments)
    function get_angle(a, b)
        return math.acos(a^b/(get_magnitude(a) * get_magnitude(b)))
    end

    -- Given a segment and direction vector, get the angle between
    function get_angle_between_segments(a1, a2, dir)
        return get_angle(get_slope(a1, a2), dir)
    end

    -- Find the area using a given h value
    function get_area(gamma, h, theta1, theta2, model)
        if h < 0.001 then
            model:warning("Tiny h")
            return 0
        end
        local temp = gamma * h + h^2/2 * (math.tan(theta1) + math.tan(theta2))
        
        if temp < 0 then
            model:warning("Negative area")
            model:warning("Gamma = " .. gamma .. "H = " .. h .. ", theta1 = " .. theta1 .. " Theta2 = " .. theta2)
            model:warning(temp)
        end
        return gamma * h + h^2/2 * (math.tan(theta1) + math.tan(theta2))
    end

    -- Make a line at a given point in the direction dir
    function get_line(point, dir, model)
        local temp_point = ipe.Vector(point.x + dir.x*20, point.y + dir.y*20)
        if point.x == nil then
            model:warning("Breaky1")
        end
        
        local temp = ipe.LineThrough(point, temp_point)
        return temp
    end

    -- Returns the inverse of the given vertex direction
    function get_inverse(dir, model)
        return ipe.Vector(-1 * dir.y, 1 * dir.x)
    end

    -- Make a line perpendicular to the given vertex direction at point (used for the lines at each vertex)
    function get_perp_line(point, dir, model)
        local temp = get_line(point, get_inverse(dir, model), model)
    return temp
    end


    -- Pass in 2 line segments and a line for the direction vector, returns the min distance between the two vectors
    -- While you could solve the direction yourself as it is orthogonal to the segments, to reduce unneeded calculations
    -- I just pass it in
    -- Will break if the lines are parallel to dir, but not if they overlap


    function get_dist_between_lines(seg1, seg2, dir, model)
        local p1 = seg1:intersects(get_line(ipe.Vector(0,0), dir, model))
        local p2 = seg2:intersects(get_line(ipe.Vector(0,0), dir, model))

        return euclidean_distance(p1, p2)
    end


    function get_dist_for_points(p1, p2, dir, model)
        --model:warning("1" .. dir.x .. dir.y)
        return get_dist_between_lines(get_perp_line(p1, dir, model), get_perp_line(p2, dir, model), dir, model)
    end




    -- Issue lies in thetas, with figuring out whether to make it positive or negative
    -- Given a location (1) and a point for direction (2), get the area of that step
    function calc_area(a1, a2, b1, b2, dir, gamma, h, model)
        local seg1 = ipe.Segment(a1, a2)
        local seg2 = ipe.Segment(b1, b2)

        local temp1 = seg1:line():dir()
        local temp2 = seg2:line():dir()


        local theta1 = get_angle(temp1, dir)
        local theta2 = get_angle(temp2, dir)

        -- Needs to adjust the angle depending on if it should increase or decrease the area
        if temp1.x*dir.y - temp1.y*dir.x < 0 then
            theta1 = -theta1
        end
        if temp2.x*dir.y - temp2.y*dir.x > 0 then
            theta2 = -theta2
        end

        return math.abs(get_area(gamma, h, theta1, theta2, model))
    end

    function reverseTable(t)
        local n = #t
        for i = 1, math.floor(n / 2) do
            t[i], t[n - i + 1] = t[n - i + 1], t[i]
        end
    end

    -- Issue lies in thetas, with figuring out whether to make it positive or negative
    -- Given a location (1) and a point for direction (2), get the area of that step
    function find_h(a1, a2, b1, b2, dir, gamma, area, model)
        if area < 0.001 then
            return 0
        end


        local seg1 = ipe.Segment(a1, a2)
        local seg2 = ipe.Segment(b1, b2)

        local temp1 = seg1:line():dir()
        local temp2 = seg2:line():dir()


        local theta1 = get_angle(temp1, dir)
        local theta2 = get_angle(temp2, dir)


        -- How do i make this work? Seems like which chain is which doesn't tend to be consistent, so it
        -- Needs to rely on orientation or be modified before the function to be consistent
        if temp1.x*dir.y - temp1.y*dir.x < 0 then
            theta1 = -theta1
        end
        if temp2.x*dir.y - temp2.y*dir.x > 0 then
            theta2 = -theta2
        end

        if math.abs(math.tan(theta1) + math.tan(theta2)) < 0.001 then
            if math.abs(gamma) < 0.001 then
                return 0
            end
            return area/gamma
        end

        local t1 = (-gamma + math.sqrt(gamma*gamma+2*(math.tan(theta1)+math.tan(theta2))*area))/(math.tan(theta1)+math.tan(theta2))
        local t2 = (-gamma - math.sqrt(gamma*gamma+2*(math.tan(theta1)+math.tan(theta2))*area))/(math.tan(theta1)+math.tan(theta2))

        if t1 < 0 then
            return t2
        else
            if t2 < 0 then
                return t1
            else
                return math.min(t1,t2)
            end
        end
    end

    -- Used for intersections
    function removeDuplicates(tbl)
        local seen = {}
        local result = {}
        
        for _, value in ipairs(tbl) do
            if not seen[value] then
                seen[value] = true
                table.insert(result, value)
            end
        end
        
        return result
    end

    --Prints out a vector
    function print_vector(vector, name, model)
        model:warning(name.." = (" .. vector.x .. ", " .. vector.y .. ")")
    end


    function run(model)

        -- Stores the unmodified coordinates in ncoordinates
        local ncoordinates = get_pts_and_segment_selection(model)
        if ncoordinates == nil or #ncoordinates == 0 then
            model:warning("No/Not enough coordinates found, exiting")
            return
        end

        -- The polygon must be convex, if not, we exit
        if not isConvex(ncoordinates) then
            model:warning("The provided shape is not convex, exiting ipelet (If three adjacent points are on the same line, it breaks)")
            return
        end

        -- Takes in the desired amount of area, doesn't accept 0 or 100 since those would do nothing
        local delta = model:getString("Enter delta value (1-99, where x means x% of the total area)")
        delta = tonumber(delta)
        if delta == nil or delta < 1 or delta > 99 then
            model:warning("Invalid delta input")
            return
        end


        -- We don't check this one since we have a default value
        local showType = model:getString("Would you like the halfspace lines (Default) or a polygon of midpoints (1)")


        local midpoints = {}
        local target_area = calculate_polygon_area(ncoordinates) * delta/100


        -- Iterates over 1 degree, 2 degrees ... 359 degrees
        for i = 0, 359 do
            -- Creates a directional vector
            local dir = ipe.Vector(math.cos(i/180*math.pi), math.sin(i/180*math.pi))

            -- Sorts the coordinates based on dot product value to the direction vector, then puts the 1st node first in the table
            local points = get_closest(ncoordinates, dir, model)
            local coordinates = reorderTable(ncoordinates, points[1])


            -- Reorders it so that chain 1 is always clockwise of dir
            if (coordinates[1].x-coordinates[2].x)*dir.y - (coordinates[1].y-coordinates[2].y)*dir.x > (coordinates[1].x-coordinates[#coordinates].x)*dir.y - (coordinates[1].y-coordinates[#coordinates].y)*dir.x then
                reverseTable(coordinates)
                coordinates = reorderTable(coordinates, points[1])
            end

            -- Gets the chains, stored as {chain1, chain2}
            local chains = get_segments_sorted(coordinates, points,  model)

            -- Indexes for chain 1 and chain 2
            local c1 = 1
            local c2 = 1

            local chain1 = chains[1]
            local chain2 = chains[2]

            -- IMPORTANT
            -- Gamma updates at the end of each segment for the following segment
            local gamma = 0
            local total_area = 0

            -- They both start on the initial point, so we just skip that useless step (Breaks index - 1 otherwise)
            if chain1[2]^dir < chain2[2]^dir then
                c1 = c1 + 1
            else
                c2 = c2 + 1
            end

            -- While neither have gone too far (This really doesn't matter, just acts as a failsafe to prevent)
            -- Infinite looping if something breaks
            while c1 <= #chain1 and c2 <= #chain2 do

                --Vectors that show the direction of a segment from a1 or b1
                local a2
                local b2

                -- Gets the length of the next segment
                local h = get_dist_for_points(points[c1+c2-2], points[c1+c2-1], dir, model)

                -- If they are just about on the same line, we don't need to account for it. Issue with == due to floating point error
                if h > 0.00001 then

                    -- We need to know which point is first, so we have two cases
                    -- We modify the second point using temp slope instead of just taking the other value in order to flip it around
                    -- This ensures angle calculation works
                    if chain1[c1]^dir < chain2[c2]^dir then   
                        a2 = chain1[c1+1]
                        b2 = chain2[c2]
                        local temp_slope = get_slope(chain2[c2-1], chain2[c2])
                        b2 = ipe.Vector(b2.x + temp_slope.x, b2.y + temp_slope.y)
                    else
                        b2 = chain2[c2+1]
                        a2 = chain1[c1]
                        local temp_slope = get_slope(chain1[c1-1], chain1[c1])
                        a2 = ipe.Vector(a2.x + temp_slope.x, a2.y + temp_slope.y)
                    end

                    -- Calculates the max area of the segment
                    local temp_area = calc_area(chain1[c1], a2, chain2[c2], b2, dir, gamma, h, model)

                    -- If it overshoots, we need to find the h that works
                    if temp_area + total_area > target_area then

                        local h = find_h(chain1[c1], a2, chain2[c2], b2, dir, gamma, target_area-total_area, model)
                        
                        -- Takes the first point and advances along dir by h, so we can find the orthogonal line to it
                        local last_point
                        if chain1[c1]^dir < chain2[c2]^dir then
                            last_point = ipe.Vector(chain1[c1].x + dir.x*h, chain1[c1].y + dir.y*h)
                        else
                            last_point = ipe.Vector(chain2[c2].x + dir.x*h, chain2[c2].y + dir.y*h)
                        end

                        -- Finds the intersection points in the main shape
                        local temp_perp_line = get_perp_line(last_point, dir, model)
                        local intersect_points = get_intersection_points(coordinates, temp_perp_line, model)


                        -- This one is possible when a line goes through a vertex
                        if #intersect_points > 2 then
                            intersect_points = removeDuplicates(intersect_points)
                        end

                        -- Shouldn't trigger, but this prevents a crash
                        if #intersect_points < 2 then
                            model:warning("Didn't intersect. You shouldn't see this message")
                            break
                        end

                        -- If half space lines, draw them
                        if showType ~= "1" then
                            local start = intersect_points[1]
                            local finish = intersect_points[2]

                            -- Create the path between the two vectors
                            local segment = {type="segment", start, finish}
                            local shape = { type="curve", closed=false, segment}
                            local pathObj = ipe.Path(model.attributes, { shape })
                            
                            -- Draw the path
                            model:creation("create basic path", pathObj)
                        end

                        -- Adds the midpoint to a table
                        table.insert(midpoints, ipe.Vector((intersect_points[1].x + intersect_points[2].x)/2,(intersect_points[1].y + intersect_points[2].y)/2))
                        
                        break

                    end

                    -- Update gamma. I wish it was simpler, but if the next point is on the same line as the first point, we need
                    -- To update the gamma there instead of at the other chain
                    if chain1[c1]^dir < chain2[c2]^dir then   
                        local temp_point
                        if chain1[c1+1]^dir < chain2[c2]^dir then
                            temp_point = get_perp_line(chain1[c1+1], dir, model):intersects(ipe.Segment(chain2[c2], chain2[c2-1]):line())
                            gamma = euclidean_distance(temp_point, chain1[c1+1])
                        else
                            temp_point = get_perp_line(chain2[c2], dir, model):intersects(ipe.Segment(chain1[c1], chain1[c1+1]):line())
                            gamma = euclidean_distance(temp_point, chain2[c2])
                        end
                    else
                        local temp_point
                        if chain2[c2+1]^dir < chain1[c1]^dir then
                            temp_point = get_perp_line(chain2[c2+1], dir, model):intersects(ipe.Segment(chain1[c1], chain1[c1-1]):line())
                            gamma = euclidean_distance(temp_point, chain2[c2+1])
                        else
                            temp_point = get_perp_line(chain1[c1], dir, model):intersects(ipe.Segment(chain2[c2], chain2[c2+1]):line())
                            gamma = euclidean_distance(temp_point, chain1[c1])
                        end

                    end

                    total_area = total_area + temp_area
                else
                    -- If h is basically 0, we can just find the distance between those two points
                    gamma = euclidean_distance(chain1[c1], chain2[c2])
                end

                


                -- Advancing on the chain (Basically, whichever segment is next we update)
                if c1 == #chain1 then
                    c2 = c2 + 1
                else
                    if c2 == #chain2 then
                        c1 = c1 + 1
                    else
                        -- Floating point issues
                        if math.abs(chain1[c1]^dir - chain2[c2]^dir) < 0.001 then
                            if chain1[c1+1]^dir < chain2[c2+1]^dir then
                                c1 = c1 + 1
                            else 
                                c2 = c2+1
                            end
                        else
                            if chain1[c1]^dir < chain2[c2]^dir then
                                c1 = c1 + 1
                            else 
                                c2 = c2+1
                            end
                        end
                    end
                end
                
            end
        end

        -- Midpoint polygon
        if showType == "1" then
            local closest_points = midpoints
            local start = closest_points[#closest_points]
            local finish
            local lines = {}
            for _, point in ipairs(closest_points) do
                finish = ipe.Vector(point.x, point.y)
                local segment = {type="segment", start, finish}
                local shape = { type="curve", closed=false, segment}
                table.insert(lines, shape)
                start = finish
            end
            local pathObj = ipe.Path(model.attributes, lines)
            model:creation("create basic path", pathObj)
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
  { label = "Polar body",           run = polar_body_run },
  { label = "Minkowski sum",        run = minkowski_run },
  { label = "Macbeath region",      run = macbeath_run },
  { label = "Floating bodies",      run = floating_bodies_run },
}
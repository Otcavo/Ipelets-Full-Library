label = "Funk and Hilbert"

about = [[
Funk and Hilbert Ipelets:
- Funk ball
- Hilbert ball
- MST with Hilbert Distance
- MST with Minimum Funk Distance
]]

-- ---------------------------------------------------------------------------
-- Funk balls
-- ---------------------------------------------------------------------------

do
    local incorrect
    local is_convex
    local is_in_polygon
    local get_polygon_vertices
    local create_segments_from_vertices
    local unique_points
    local get_polygon_vertices_and_segments
    local copy_table
    local get_pt_and_polygon_selection
    local create_ray
    local create_rays
    local intersect
    local equal
    local get_spokes
    local get_spokes_path_objs
    local compute_K_constants
    local get_point_vertical
    local get_point_vertical_reverse
    local get_point_forward
    local get_point_reverse
    local get_points_on_spokes_forward
    local get_points_on_spokes_reverse
    local create_shape_from_vertices
    local orient
    local convex_hull
    local not_in_table

    function incorrect(title, model) model:warning(title) end

    function is_convex(vertices)
        local _, convex_hull_vectors = convex_hull(vertices)
        return #convex_hull_vectors == #vertices
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
        local segments_start_finish = {}
        for i=1, #vertices-1 do
            table.insert( segments, ipe.Segment(vertices[i], vertices[i+1]) )
            table.insert( segments_start_finish, {vertices[i],vertices[i+1]} )
        end

        table.insert( segments, ipe.Segment(vertices[#vertices], vertices[1]) )
        table.insert( segments_start_finish, {vertices[#vertices], vertices[1]} )
        return segments, segments_start_finish
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
        local segments, segments_start_finish = create_segments_from_vertices(vertices)
        return vertices, segments, segments_start_finish
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

        local point = referenceObject:matrix() * referenceObject:position()
        local vertices, segments, segments_start_finish = get_polygon_vertices_and_segments(pathObject, model)

        local poly1_convex = is_convex(copy_table(vertices))
        if poly1_convex == false then incorrect("Polygon must be convex", model) return end
        if not is_in_polygon(point, copy_table(vertices)) then incorrect("Point must be inside the polygon", model) return end

        return point, vertices, segments, segments_start_finish
    end

    function create_ray(v,c, model)
        return ipe.LineThrough(v, c)
    end

    function create_rays(v, c, model)
        local rays = {}
        for i=1, #v do table.insert( rays, { v[i],  create_ray(v[i], c, model) } ) end
        return rays
    end

    function intersect(s,r)
        return r:intersects(s)
    end

    function equal(v1,v2)
        if v1.x == v2.x and v1.y == v2.y then
            return true
        else
            return false
        end
    end

    function get_spokes(rays, segments, model)

    local spokes = {}
    local vertex_intersect = {}
    for j=1, #segments do
        for i=1, #rays do
            -- local line = ipe.LineThrough(segments[j][1],segments[j][2])
            if not equal(segments[j][1], rays[i][1]) and not equal(segments[j][2], rays[i][1]) then
                local intersection = ipe.Segment(segments[j][1],segments[j][2]):intersects(rays[i][2])
                if intersection ~= nil then
                    table.insert( vertex_intersect, {rays[i][1].x, rays[i][1].y, intersection.x, intersection.y} )
                    table.insert( spokes, {type="segment", rays[i][1], intersection} )
                end
            end
        end
    end
    return spokes, vertex_intersect
    end

    function get_spokes_path_objs(spokes, model)
        local spoke_obj_list = {}
        for _, spoke in ipairs(spokes) do
            local shape = {type="curve", closed=true; spoke}
            local obj = ipe.Path(model.attributes, { shape })
            table.insert(spoke_obj_list, obj)
        end
        return spoke_obj_list
    end

    function compute_K_constants(a,d)
        local top_fraction = (a.x-d.x)*(a.x-d.x) + (a.y-d.y)*(a.y-d.y)
        local bottom_fraction = math.exp(2)
        local K1 = top_fraction / bottom_fraction
        local K2 = (a.y-d.y) / (a.x-d.x)
        local K3 = (math.exp(1)*math.exp(1)) * ((a.x-d.x)*(a.x-d.x) + (a.y-d.y)*(a.y-d.y))
        return K1, K2, K3
    end

    function get_point_vertical(a,d)
        local K1,_ = compute_K_constants(a,d)
        local A = 1
        local B = -2*d.y
        local C = (a.x*a.x) - 2*d.x*a.x + (d.x*d.x) + (d.y*d.y) - K1
        local x
        local y
        local discriminant = (B*B)-4*A*C
        
        if d.y > a.y then
            x = a.x
            y = (-B - math.sqrt(discriminant)) / (2*A)
        else
            x = a.x
            y = (-B + math.sqrt(discriminant)) / (2*A)
        end
        
        return ipe.Vector(x,y)
    end

    function get_point_vertical_reverse(a,d)
        local _, _, K3 = compute_K_constants(a,d)
        local A = 1
        local B = -2*d.y
        local C = (d.y*d.y) - K3
        local x
        local y
        local discriminant = (B*B)-4*A*C
        
        if a.y > d.y then
            x = a.x
            y = (-B + math.sqrt(discriminant)) / (2*A)
        else
            x = a.x
            y = (-B - math.sqrt(discriminant)) / (2*A)
        end
        
        return ipe.Vector(x,y)
    end

    function get_point_forward(a,d)
        local K1,K2, _ = compute_K_constants(a,d)
        local A = 1 + (K2*K2)
        local B = -2*d.x-2*(K2*K2)*a.x+2*K2*a.y-2*d.y*K2
        local C = (d.x*d.x) + (K2*K2) * (a.x*a.x) + (a.y*a.y) - 2*K2*a.x*a.y+2*d.y*K2*a.x-2*d.y*a.y+(d.y*d.y)-K1

        local discriminant = (B*B)-4*A*C
        local x_1 = (-B + math.sqrt(discriminant)) / (2*A)
        local y_1 = K2*(x_1-a.x)+a.y

        local x_2 = (-B - math.sqrt(discriminant)) / (2*A)
        local y_2 = K2*(x_2-a.x)+a.y

        return ipe.Vector(x_1,y_1), ipe.Vector(x_2, y_2)
    end

    function get_point_reverse(a,d)
        local _,K2,K3 = compute_K_constants(a,d)
        local A = 1 + (K2*K2)
        local B = -2 * d.x - 2 * (K2*K2) * d.x
        local C = (d.x*d.x) + (K2*K2) * (d.x*d.x) - K3

        local discriminant = (B*B)-4*A*C
        local x_1 = (-B + math.sqrt(discriminant)) / (2*A)
        local y_1 = K2*(x_1-a.x)+a.y

        local x_2 = (-B - math.sqrt(discriminant)) / (2*A)
        local y_2 = K2*(x_2-a.x)+a.y

        return ipe.Vector(x_1,y_1), ipe.Vector(x_2, y_2)

    end

    function get_points_on_spokes_forward(vertex_intersect, center)
        
        local points_on_spokes = {}

        for i=1, #vertex_intersect do

            local vertex = ipe.Vector(vertex_intersect[i][1], vertex_intersect[i][2])

            if center.x == vertex.x then
                local point = get_point_vertical(center, vertex)
                table.insert(points_on_spokes, point)
            else
                local point1, point2 = get_point_forward(center, vertex)
                
                if vertex.x > center.x  then
                    table.insert(points_on_spokes, point2)
                else
                    table.insert(points_on_spokes, point1)
                end
            end
        end
        return points_on_spokes
    end


    function get_points_on_spokes_reverse(vertex_intersect, center, model)
        
        local points_on_spokes = {}

        for i=1, #vertex_intersect do

            local intersection = ipe.Vector(vertex_intersect[i][1], vertex_intersect[i][2])
            if center.x == intersection.x then
                local point = get_point_vertical_reverse(center, intersection)
                table.insert(points_on_spokes, point)
            else
                local point1, point2 = get_point_reverse(center, intersection, model)
                if center.x < intersection.x then
                    table.insert(points_on_spokes, point2)
                else
                    table.insert(points_on_spokes, point1)
                end
            end

        end
        return points_on_spokes
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

    function run_forward_spokes(model)
        if not get_pt_and_polygon_selection(model) then return end
        local center, vertices, _, segments_start_finish = get_pt_and_polygon_selection(model)
        local v_rays = create_rays(vertices, center, model)
        local spokes, vertex_intersect = get_spokes(v_rays, segments_start_finish, model)
        
        local spoke_obj_list = get_spokes_path_objs(spokes, model)
        local points_on_spokes = get_points_on_spokes_forward(vertex_intersect, center)
        local shape, _ = convex_hull(points_on_spokes)
        table.insert(spoke_obj_list, ipe.Path(model.attributes, { shape }))
        model:creation("points on spokes", ipe.Group(spoke_obj_list) )
    end

    function run_forward_without_spokes(model)
        if not get_pt_and_polygon_selection(model) then return end
        local center, vertices, _, segments_start_finish = get_pt_and_polygon_selection(model)
        local v_rays = create_rays(vertices, center, model)
        local spokes, vertex_intersect = get_spokes(v_rays, segments_start_finish, model)
        
        local points_on_spokes = get_points_on_spokes_forward(vertex_intersect, center)
        local shape, _ = convex_hull(points_on_spokes)
        
        model:creation("points on spokes", ipe.Path(model.attributes, { shape }))
    end

    function run_reverse_spokes(model)
        if not get_pt_and_polygon_selection(model) then return end
        local center, vertices, _, segments_start_finish = get_pt_and_polygon_selection(model)
        local v_rays = create_rays(vertices, center, model)
        local spokes, vertex_intersect = get_spokes(v_rays, segments_start_finish, model)
        
        local spoke_obj_list = get_spokes_path_objs(spokes, model)
        local points_on_spokes = get_points_on_spokes_reverse(vertex_intersect, center)
        local shape, _ = convex_hull(points_on_spokes)
        table.insert(spoke_obj_list, ipe.Path(model.attributes, { shape }))
        model:creation("points on spokes", ipe.Group(spoke_obj_list) )
    end

    function run_reverse_without_spokes(model)
        if not get_pt_and_polygon_selection(model) then return end
        local center, vertices, _, segments_start_finish = get_pt_and_polygon_selection(model)
        local v_rays = create_rays(vertices, center, model)
        local spokes, vertex_intersect = get_spokes(v_rays, segments_start_finish, model)
        
        local points_on_spokes = get_points_on_spokes_reverse(vertex_intersect, center)
        local shape, _ = convex_hull(points_on_spokes)
        
        model:creation("points on spokes", ipe.Path(model.attributes, { shape }))
    end
end

-- ---------------------------------------------------------------------------
-- Hilbert balls
-- ---------------------------------------------------------------------------

do
    local incorrect
    local is_convex
    local is_in_polygon
    local get_polygon_vertices
    local create_segments_from_vertices
    local unique_points
    local get_polygon_vertices_and_segments
    local copy_table
    local get_pt_and_polygon_selection
    local create_ray
    local create_rays
    local intersect
    local equal
    local get_spokes
    local get_spokes_path_objs
    local compute_K_constants
    local get_point_vertical
    local get_point
    local get_points_on_spokes
    local create_shape_from_vertices
    local orient
    local convex_hull
    local not_in_table

    function incorrect(title, model) model:warning(title) end

    function not_in_table(vectors, vector_comp)
        local flag = true
        for _, vertex in ipairs(vectors) do
            if vertex == vector_comp then
                flag = false
            end
        end
        return flag
    end

    function is_convex(vertices)
        local _, convex_hull_vectors = convex_hull(vertices)
        return #convex_hull_vectors == #vertices
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
        local segments_start_finish = {}
        for i=1, #vertices-1 do
            table.insert( segments, ipe.Segment(vertices[i], vertices[i+1]) )
            table.insert( segments_start_finish, {vertices[i],vertices[i+1]} )
        end

        table.insert( segments, ipe.Segment(vertices[#vertices], vertices[1]) )
        table.insert( segments_start_finish, {vertices[#vertices], vertices[1]} )
        return segments, segments_start_finish
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
        local segments, segments_start_finish = create_segments_from_vertices(vertices)
        return vertices, segments, segments_start_finish
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

        local point = referenceObject:matrix() * referenceObject:position()
        local vertices, segments, segments_start_finish = get_polygon_vertices_and_segments(pathObject, model)

        local poly1_convex = is_convex(copy_table(vertices))
        if poly1_convex == false then incorrect("Polygon must be convex", model) return end
        if not is_in_polygon(point, copy_table(vertices)) then incorrect("Point must be inside the polygon", model) return end

        return point, vertices, segments, segments_start_finish
    end

    function create_ray(v,c, model)
        return ipe.LineThrough(v, c)
    end

    function create_rays(v, c, model)
        local rays = {}
        for i=1, #v do table.insert( rays, { v[i],  create_ray(v[i], c, model) } ) end
        return rays
    end

    function intersect(s,r)
        return r:intersects(s)
    end

    function equal(v1,v2)
        if v1.x == v2.x and v1.y == v2.y then
            return true
        else
            return false
        end
    end

    function get_spokes(rays, segments, model)

    local spokes = {}
    local vertex_intersect = {}
    for j=1, #segments do
        for i=1, #rays do
            if not equal(segments[j][1], rays[i][1]) and not equal(segments[j][2], rays[i][1]) then
                local intersection = ipe.Segment(segments[j][1],segments[j][2]):intersects(rays[i][2])
                if intersection ~= nil then
                    table.insert( vertex_intersect, {rays[i][1].x, rays[i][1].y, intersection.x, intersection.y} )
                    table.insert( spokes, {type="segment", rays[i][1], intersection} )
                end
            end
        end
    end
    return spokes, vertex_intersect
    end

    function get_spokes_path_objs(spokes, model)
        local spoke_obj_list = {}
        for _, spoke in ipairs(spokes) do
            local shape = {type="curve", closed=true; spoke}
            local obj = ipe.Path(model.attributes, { shape })
            table.insert(spoke_obj_list, obj)
        end
        return spoke_obj_list
    end

    function compute_K_constants(a,b,d,R, model)
        local top_fraction = (b.x-d.x)*(b.x-d.x) + (b.y-d.y)*(b.y-d.y)
        local bottom_fraction = (b.x-a.x)*(b.x-a.x) + (b.y-a.y)*(b.y-a.y)
        local square_root_part = math.sqrt(top_fraction/bottom_fraction)
        local exponent_part = 2*(2*R - math.log(square_root_part))
        local K1 = math.exp(exponent_part)
        local K2 = (a.y-d.y) / (a.x-d.x)
        return K1, K2
    end

    function get_point_vertical(a,b,d,R, model) 
        local K1,_ = compute_K_constants(a,b,d,R, model)
        local A = 1-K1
        local B = -2*a.y+2*d.y*K1
        local C_1 = (a.x*a.x) - 2*(a.x*a.x) - K1*(a.x*a.x) + 2*K1*d.x*a.x
        local C_2 = K1*((d.x*d.x)+(d.y*d.y))-(a.x*a.x)-(a.y*a.y)
        local C = C_1 - C_2
        local x
        local y
        local discriminant = (B*B)-4*A*C
        
        if d.y > a.y then
            x = a.x
            y = (-B + math.sqrt(discriminant)) / (2*A)
        else
            x = a.x
            y = (-B - math.sqrt(discriminant)) / (2*A)
        end
        
        return ipe.Vector(x,y)
    end

    function get_point(a,b,d,R, model)
        local K1,K2 = compute_K_constants(a,b,d,R, model)
        local A = 1 + (K2*K2) - K1 - K1*(K2*K2)
        local B = -2*a.x-2*a.x*(K2*K2)+2*a.y*K2-2*a.y*K2+2*K1*d.x+2*K1*a.x*(K2*K2)-2*K1*a.y*K2+2*d.y*K1*K2
        local C_1 = (a.x*a.x)*(K2*K2)-2*a.x*a.y*K2+(a.y*a.y)-2*a.y*(-K2*a.x+a.y)-K1*((a.x*a.x)*(K2*K2)-2*a.x*a.y*K2+(a.y*a.y))+2*d.y*K1*(-K2*a.x+a.y)
        local C_2 = K1*((d.x*d.x)+(d.y*d.y))-(a.x*a.x)-(a.y*a.y)
        local C = C_1 - C_2

        local discriminant = (B*B)-4*A*C
        local x_1 = (-B + math.sqrt(discriminant)) / (2*A)
        local y_1 = K2*(x_1-a.x)+a.y

        local x_2 = (-B - math.sqrt(discriminant)) / (2*A)
        local y_2 = K2*(x_2-a.x)+a.y

        return ipe.Vector(x_1,y_1), ipe.Vector(x_2, y_2)
    end
    
    function get_points_on_spokes(vertex_intersect, center, radius, polygon, model)
        
        local points_on_spokes = {}

        for i=1, #vertex_intersect do

            local vertex = ipe.Vector(vertex_intersect[i][1], vertex_intersect[i][2])
            local intersect = ipe.Vector(vertex_intersect[i][3], vertex_intersect[i][4])

            if intersect.x == vertex.x then
                local point1 = get_point_vertical(vertex, center, intersect, radius, model)
                local point2 = get_point_vertical(intersect, center, vertex, radius, model)
                table.insert(points_on_spokes, point1)
                table.insert(points_on_spokes, point2)
            else
                local point1, point2 = get_point(vertex, center, intersect, radius, model)
                local point3, point4 = get_point(intersect, center, vertex, radius, model)
                
                if vertex.x < intersect.x  then
                    table.insert(points_on_spokes, point1)
                else
                    table.insert(points_on_spokes, point2)
                end
                if intersect.x < vertex.x then
                    table.insert(points_on_spokes, point3)
                else
                    table.insert(points_on_spokes, point4)
                end
            end
            

        end
        return points_on_spokes
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

    function run_spokes(model)
        if not get_pt_and_polygon_selection(model) then return end
        local center, vertices, _, segments_start_finish = get_pt_and_polygon_selection(model)
        local v_rays = create_rays(vertices, center, model)
        local spokes, vertex_intersect = get_spokes(v_rays, segments_start_finish, model)
        
        local spoke_obj_list = get_spokes_path_objs(spokes, model)
        local points_on_spokes = get_points_on_spokes(vertex_intersect, center, 1, vertices, model)
        local shape, _ = convex_hull(points_on_spokes)
        table.insert(spoke_obj_list, ipe.Path(model.attributes, { shape }))

        model:creation("points on spokes", ipe.Group(spoke_obj_list) )
    end

    function run_without_spokes(model)
        if not get_pt_and_polygon_selection(model) then return end
        local center, vertices, _, segments_start_finish = get_pt_and_polygon_selection(model)
        local v_rays = create_rays(vertices, center, model)
        local _, vertex_intersect = get_spokes(v_rays, segments_start_finish, model)
        
        local points_on_spokes = get_points_on_spokes(vertex_intersect, center, 1, vertices, model)
        local shape, _ = convex_hull(points_on_spokes)

        model:creation("points on spokes", ipe.Path(model.attributes, { shape }) )
    end
end

-- ---------------------------------------------------------------------------
-- MSTs
-- ---------------------------------------------------------------------------
do
    local print_vertices
    local print_table
    local print_vertex
    -- local print
    local findIntersectionsWithShape
    local minKey
    local primMST
    local orient
    local convex_hull
    local is_convex
    local get_polygon_vertices
    local copy_table
    local getSelectedPointsAndShape
    local getSelectedPoints
    local isPointInShape

    --! PRINT FUNCTIONS
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


    --! DISTANCE FUNCTIONS

    --! Euclidean Distance (2 points)
    function euclideanDistance(model, shape_vertices, p, q)
        return (p - q):len()
    end

    --! Helper function to find the intersections of a ray through a line segment with a polygon
    function findIntersectionsWithShape(ray, shape_vertices)
        local intersectionPoints = {}
        
        -- Go through each edge of the shape and check for intersection with the ray
        for i = 1, #shape_vertices do
            local next_index = (i % #shape_vertices) + 1
            local edge = {shape_vertices[i], shape_vertices[next_index]}
            local intersection = ipe.Segment(edge[1], edge[2]):intersects(ray)
            
            -- If there is an intersection, add it to the list of intersection points
            if intersection then
                table.insert(intersectionPoints, intersection)
            end
        end
        
        -- Sort intersection points by distance from the starting point of the ray
        table.sort(intersectionPoints, function(p, q)
            -- Calculate the distance of point 'p' from the ray's starting point
            local distP = (p - ray:point()):len()
            -- Calculate the distance of point 'q' from the ray's starting point
            local distQ = (q - ray:point()):len()
            -- Compare the distances; if distP is less than distQ, then 'p' should come before 'q' in the sorted list
            return distP < distQ
        end)
        
        return intersectionPoints -- Return the sorted list of intersection points
    end

    --! HILBERT DISTANCE (shape, 4 points)
    -- d_{H_K}(B, C) = \frac{1}{2}\ln[\frac{\lVert D-B \rVert}{\lVert D-C \rVert} \cdot \frac{\lVert A-C \rVert}{\lVert A-B \rVert}]
    function hilbertDistance(model, shape_vertices, B, C)
        -- Create the rays that extend in both directions from B to C and from C to B
        local BC = ipe.LineThrough(B, C)
        
        -- Find the intersection points with the shape
        local intersections = findIntersectionsWithShape(BC, shape_vertices)
        -- A is the first intersection, D is the last
        local A, D = intersections[1], intersections[#intersections]
        
        -- Compute the Euclidean distances needed for the Hilbert distance formula
        local DB = (D - B):len()
        local DC = (D - C):len()
        local AC = (A - C):len()
        local AB = (A - B):len()
        
        -- Calculate the Hilbert distance using the formula
        local hilbertDist1 = 0.5 * math.log((DB / DC) * (AC / AB))
        local hilbertDist2 = 0.5 * math.log((AB / AC) * (DC / DB))

        if hilbertDist1 > 0 then return hilbertDist1 else return hilbertDist2 end
    end


    --! REVERSE FUNK DISTANCE (shape, 3 points)
    -- d_{F_\Omega} (B, C) = ln(\frac{\lVert A-B \rVert}{\lVert A-C \rVert})
    function reverseFunkDistance(model, shape_vertices, B, C)
        local BC = ipe.LineThrough(B, C) -- Ray from B through C

        -- D is the point on the shape reached by shooting a ray from B through C
        local intersections = findIntersectionsWithShape(BC, shape_vertices)
        local A = intersections[1] -- Farthest intersection point

        -- Calculate Euclidean distances
        local AB = (A - B):len()
        local AC = (A - C):len()

        -- Calculate Forward Funk distance
        local reverseFunkDist1 = math.log(AB / AC)
        local reverseFunkDist2 = math.log(AC / AB)

        if reverseFunkDist1 > 0 then return reverseFunkDist1 else return reverseFunkDist2 end
    end

    --! FORWARD FUNK DISTANCE (shape, 3 points)
    -- d_{F_\Omega} (B, C) = ln(\frac{\lVert D-C \rVert}{\lVert D-B \rVert})
    function forwardFunkDistance(model, shape_vertices, B, C)
        local BC = ipe.LineThrough(B, C) -- Ray from C through B

        -- A is the point on the shape reached by shooting a ray from C through B
        local intersections = findIntersectionsWithShape(BC, shape_vertices)
        local D = intersections[#intersections] -- Farthest intersection point

        -- Calculate Euclidean distances
        local DC = (D - C):len()
        local DB = (D - B):len()

        -- Calculate Reverse Funk distance
        local forwardFunkDist1 = math.log(DC / DB)
        local forwardFunkDist2 = math.log(DB / DC)

        if forwardFunkDist1 > 0 then return forwardFunkDist1 else return forwardFunkDist2 end
    end

    --! MINIMUM FUNK DISTANCE
    -- d_{F_{{min}_\Omega}} (B, C) = {min}(d_{F_\Omega} (B, C), d_{F_\Omega} (C, B))
    -- necessary since Funk is inherently directional
    function minimumFunkDistance(model, shape_vertices, B, C)
        local forwardFunkDist = forwardFunkDistance(model, shape_vertices, B, C)
        local reverseFunkDist = reverseFunkDistance(model, shape_vertices, B, C)

        -- print forward and reverse funk distances

        return math.min(forwardFunkDist, reverseFunkDist)
    end


    --! PRIM'S ALGORITHM
    -- Modular distance function
    function distance(func, model, shape_vertices, ...)
        local dist = func(model, shape_vertices, ...)
        return dist
    end

    -- Utility function to find the vertex with minimum key value from the set of vertices not yet included in MST
    function minKey(key, mstSet)
        local min = math.huge
        local min_index = -1

        -- Search not yet included vertices
        for v, k in ipairs(key) do
            if k < min and not mstSet[v] then
                min = k
                min_index = v
            end
        end

        return min_index
    end

    -- Function to construct the Minimum Spanning Tree (MST) using Prim's algorithm
    -- vertices: An array of points within the specified shape
    -- shape: The shape within which the points lie, used for distance calculations
    function primMST(points, shape_vertices, distance_func, model)
        local V = #points -- Number of vertices
        local parent = {} -- Array to store the MST
        local key = {} -- Key values used to pick minimum weight edge
        local mstSet = {} -- To represent set of vertices not yet included in MST
        local computedEdges = {} -- To store computed edges

        -- Initialize all keys as infinite and mstSet[] as false
        for i = 1, V do
            key[i] = math.huge
            mstSet[i] = false
        end

        -- Initialize all edge pairs as not computed
        for i = 1, V do
            key[i] = math.huge
            mstSet[i] = false
            computedEdges[i] = {}
            for j = 1, V do
                computedEdges[i][j] = false
            end
        end

        -- First vertex in the MST will always be the first point
        key[1] = 0
        parent[1] = -1 -- First node is always root of MST

        for count = 1, V do
            -- Pick the minimum key vertex from the set of vertices not yet included in MST
            local u = minKey(key, mstSet)

            -- Add the picked vertex to the MST set
            mstSet[u] = true

            -- Update key value and parent index of the adjacent vertices of the picked vertex
            for v = 1, V do
                -- Only update the key if vertices[v] is not in mstSet, there is an edge from u to v,
                -- and weight of edge from u to v is smaller than key[v]
                -- Distance calculation can be modular, based on the shape's requirements

                if u ~= v and not computedEdges[u][v] and not computedEdges[v][u] then -- checks if the vertices are the same

                    -- Ensure points[u] is top-leftmost and points[v] is bottom-rightmost
                    local pointU = points[u]
                    local pointV = points[v]
                    if pointV.x < pointU.x or (pointV.x == pointU.x and pointV.y < pointU.y) then
                        -- Swap points if pointV is top-leftmost
                        pointU, pointV = pointV, pointU
                    end

                    local weight = distance(distance_func, model, shape_vertices, pointU, pointV)
                    -- print(weight, "Weight", model)
                    
                    -- weight > 0 ensures that the distance is not 0, i.e. the vertices are not the same
                    if weight > 0 and not mstSet[v] and key[v] > weight then
                        key[v] = weight
                        parent[v] = u
                        computedEdges[u][v] = true
                        computedEdges[v][u] = true
                    end
                end
            end
        end

        -- parent array now represents the MST
        return parent
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

        return S
    end


    --! IPE FUNCTIONS
    function is_convex(vertices)
        local convex_hull_vectors = convex_hull(vertices)
        return #convex_hull_vectors == #vertices
    end

    --[=[
    Given:
    - path obj: () -> Path
    Return:
    - table of vertices: () -> {Vector}
    --]=]
    function get_polygon_vertices(obj, model)

        local shape = obj:shape()
        local polygon = obj:matrix()

        vertices = {}

        -- Apply transformation to the first vertex to handle translation
        vertex = polygon * shape[1][1][1]
        table.insert(vertices, vertex)

        -- Apply transformation to the rest of the vertices to handle translation
        for i=1, #shape[1] do
            vertex = polygon * shape[1][i][2]
            table.insert(vertices, vertex)
        end

        return vertices
    end

    function copy_table(orig_table)
        local new_table = {}
        for i=1, #orig_table do new_table[i] = orig_table[i] end
        return new_table
    end

    -- Function to collect vertices and a shape from Ipe selection
    function getSelectedPointsAndShape(model)
        local points = {}
        local shapeObj
        local p = model:page()
        if not p:hasSelection() then
            model:warning("noselection") -- explain and quit if no selection
            return
        end

        for _, obj, sel, _ in p:objects() do
            if sel then
                if obj:type() == "reference" then
                    -- Apply transformation to the point's position to handle translation
                    local transformedPoint = obj:matrix() * obj:position()
                    table.insert(points, transformedPoint)
                elseif obj:type() == "path" then
                    shapeObj = obj
                end
            end
        end

        if not shapeObj then
            model:warning("Please select a shape")
            return
        end

        if #points < 2 then
            model:warning("Please select at least two points")
            return
        end

        local shape_vertices = get_polygon_vertices(shapeObj, model)
        local shape_convex = is_convex(copy_table(shape_vertices))
        if shape_convex == false then model:warning("Polygon must be convex") return end

        return points, shape_vertices
    end

    -- Function to collect vertices from Ipe selection
    function getSelectedPoints(model)
        local points = {}
        local p = model:page()
        if not p:hasSelection() then
            model:warning("noselection") -- explain and quit if no selection
            return
        end

        for _, obj, sel, _ in p:objects() do
            if sel then
                if obj:type() == "reference" then
                    -- Apply transformation to the point's position to handle translation
                    local transformedPoint = obj:matrix() * obj:position()
                    table.insert(points, transformedPoint)
                end
            end
        end

        if #points < 2 then
            model:warning("Please select at least two points")
            return
        end

        return points
    end

    -- Function to determine if a point is inside a polygon
    --[=[
    Given:
    - point: () -> Vector
    - vertices of a polygon: () -> {Vector}
    Return:
    - returns true if point is inside the polygon, false otherwise
    - if the point is on the edge of a polygon, then false is returned
    - () -> Bool
    --]=]
    function isPointInShape(point, shape_vertices)
        local x, y = point.x, point.y
        local j = #shape_vertices
        local inside = false

        for i = 1, #shape_vertices do
            local xi, yi = shape_vertices[i].x, shape_vertices[i].y
            local xj, yj = shape_vertices[j].x, shape_vertices[j].y

            if ((yi > y) ~= (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi) then
                inside = not inside
            end
            j = i
        end

        return inside
    end

    --! RUN FUNCTIONS FOR GIVEN DISTANCE FUNCTION
    function runWithDistanceFunc(model, distFunc)
        if not getSelectedPoints(model) then return end
        local points = getSelectedPoints(model)

        -- print_vertices(points, "Selected points", model)


        -- Compute the MST
        local mst = primMST(points, {}, distFunc, model)

        -- Print the MST
        -- print_table(mst, "MST", model)
        
        -- Visualize the MST
        local paths = {}
        for i, parentIndex in ipairs(mst) do
            if parentIndex ~= -1 then
                -- Create a segment for each edge in the MST
                local segment = {type="segment", points[i], points[parentIndex]}
                local shape = {type="curve", closed=false, segment}  -- Create an open curve with the single segment
                local path = ipe.Path(model.attributes, {shape})
                table.insert(paths, path)
            end
        end
        
        model:creation("Create MST edges", ipe.Group(paths))
    end
    function runWithDistanceShape(model, distFunc)
        if not getSelectedPointsAndShape(model) then return end
        local points, shape_vertices = getSelectedPointsAndShape(model)

        -- print_vertices(points, "Selected points", model)
        -- print_vertices(shape_vertices, "Shape vertices", model)

        -- Check if all points are within the shape
        for _, point in ipairs(points) do
            if not isPointInShape(point, shape_vertices) then
                model:warning("All points must be within the selected shape")
                return
            end
        end


        -- Compute the MST
        local mst = primMST(points, shape_vertices, distFunc, model)

        -- Print the MST
        -- print_table(mst, "MST", model)
        
        -- Visualize the MST
        local paths = {}
        for i, parentIndex in ipairs(mst) do
            if parentIndex ~= -1 then
                -- Create a segment for each edge in the MST
                local segment = {type="segment", points[i], points[parentIndex]}
                local shape = {type="curve", closed=false, segment}  -- Create an open curve with the single segment
                local path = ipe.Path(model.attributes, {shape})
                table.insert(paths, path)
            end
        end
        
        model:creation("Create MST edges", ipe.Group(paths))
    end

end

methods = {
    { label="Forward Funk ball with spokes", run = run_forward_spokes },
    { label="Forward Funk ball without spokes", run = run_forward_without_spokes },
    { label="Reverse Funk ball with spokes", run = run_reverse_spokes },
    { label="Reverse Funk ball without spokes", run = run_reverse_without_spokes },
    { label="Hilbert ball with spokes", run = run_spokes },
    { label="Hilbert ball without spokes", run = run_without_spokes },
    { label = "MST with Hilbert Distance", run = function(model) runWithDistanceShape(model, hilbertDistance) end },
    { label = "MST with Minimum Funk Distance", run = function(model) runWithDistanceShape(model, minimumFunkDistance) end },
}


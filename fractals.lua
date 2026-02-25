label = "Fractals"

about = [[
Fractal Ipelets:
- Sierpinski triangle
- Sierpinski carpet
- Dragon curve
]]

-- ---------------------------------------------------------------------------
-- Helpers 
-- ---------------------------------------------------------------------------

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

function create_shape_from_vertices(v, model)
    local shape = { type = "curve", closed = true, }
    for i = 1, #v - 1 do
        table.insert(shape, { type = "segment", v[i], v[i + 1] })
    end
    table.insert(shape, { type = "segment", v[#v], v[1] })
    return shape
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

function get_selected_poly_vertices(model)
    local shape
    local page = model:page()

    if not page:hasSelection() then
        return
    end

    for _, obj, sel, _ in page:objects() do
        if sel and obj:type() == "path" then
            shape = obj
        end
    end

    return get_polygon_vertices(shape, model)
end

function unique_points(points, model)
    -- Check for duplicate points and remove them
    local uniquePoints = {}
    if points == nil then 
        incorrect("No points selected", model)
        return 
    end

    for i = 1, #points do
        if (not_in_table(uniquePoints, points[i])) then
                    table.insert(uniquePoints, points[i])
                end
    end
    return uniquePoints
end

function merge(t1, t2)
    local merged = {}

    for _, v in ipairs(t1) do
        merged[#merged + 1] = v
    end
    for _, v in ipairs(t2) do
        merged[#merged + 1] = v
    end

    return merged
end

function add_vectors(a, b) return ipe.Vector((a.x + b.x), (a.y + b.y)) end

function midpoint(a, b) return ipe.Vector((a.x + b.x) / 2, (a.y+ b.y) / 2) end

function sub_vectors(a, b) return ipe.Vector((a.x - b.x), (a.y - b.y)) end

function scale_vector(v, scale) return ipe.Vector(v.x * scale, v.y * scale) end

function is_L(v)
    if #v ~= 3 then
        return false
    else

        -- gets lengths of sides
        local s1 = (v[1].x - v[2].x) ^ 2 + (v[1].y - v[2].y) ^ 2
        local s2 = (v[2].x - v[3].x) ^ 2 + (v[2].y - v[3].y) ^ 2

        -- dot product (for right angle check)
        local dp = (v[1].x - v[2].x)*(v[2].x - v[3].x) + (v[1].y - v[2].y)*(v[2].y - v[3].y)

        return s1 == s2 and dp == 0
    end
end


function is_square(v, model)
    if #v ~= 4 then
        return false
    else
        local s1 = (v[1].x - v[2].x) ^ 2 + (v[1].y - v[2].y) ^ 2
        local s2 = (v[2].x - v[3].x) ^ 2 + (v[2].y - v[3].y) ^ 2
        local s3 = (v[3].x - v[4].x) ^ 2 + (v[3].y - v[4].y) ^ 2
        local s4 = (v[1].x - v[4].x) ^ 2 + (v[1].y - v[4].y) ^ 2

        local dp = (v[1].x - v[2].x)*(v[2].x - v[3].x) + (v[1].y - v[2].y)*(v[2].y - v[3].y)

        return s1 == s2 and s2 == s3 and s3 == s4 and dp == 0
    end
end

-- ---------------------------------------------------------------------------
-- Sierpinski's Triangle
-- ---------------------------------------------------------------------------

function sierpinski_triangle(v, iterations, model)
    if iterations > 0 then
        local p12 = midpoint(v[1], v[2])
        local p13 = midpoint(v[1], v[3])
        local p23 = midpoint(v[2], v[3])
        
        local shape = create_shape_from_vertices({ p12, p13, p23 }, model)

        model:creation(1, ipe.Path(model.attributes, {shape}))
        
        sierpinski_triangle({ v[1], p12, p13 }, iterations - 1, model)
        sierpinski_triangle({ v[2], p12, p23 }, iterations - 1, model)
        sierpinski_triangle({ v[3], p13, p23 }, iterations - 1, model)
    end 
end

function sierpinski_triangle_run(model)
    local v = get_selected_poly_vertices(model)
    local vu = unique_points(v, model)

    local depth = tonumber(model:getString("Enter Depth"))

    if vu == nil then 
        incorrect("waiter! waiter! more vertices, please!", model)
        return
    end

    if #vu ~= 3 then
        incorrect("waiter! waiter! i need a triangle!", model)
        return
    end

    sierpinski_triangle(vu, depth, model)
end

-- ---------------------------------------------------------------------------
-- Sierpinski's Carpet
-- ---------------------------------------------------------------------------

function sierpinski_carpet(v, iterations, model)

    if iterations > 0 then
        -- distance vectors for splitting v[1] to v[2] to thirds
        local p1to2d = ipe.Vector((v[2].x - v[1].x) / 3, (v[2].y - v[1].y) / 3)
        -- distance vectors for splitting v[2] to v[3] to thirds
        local p2to3d = ipe.Vector((v[3].x - v[2].x) / 3, (v[3].y - v[2].y) / 3)
        -- distance vectors for splitting v[3] to v[4] to thirds
        local p3to4d = ipe.Vector((v[4].x - v[3].x) / 3, (v[4].y - v[3].y) / 3)
        -- distance vectors for splitting v[4] to v[1] to thirds
        local p4to1d = ipe.Vector((v[1].x - v[4].x) / 3, (v[1].y - v[4].y) / 3)
        -- Vectors above get added to the corners of the square to get each mid point

        local m1 = add_vectors(add_vectors(v[1], p1to2d), p2to3d)
        local m2 = add_vectors(add_vectors(v[2], p2to3d), p3to4d)
        local m3 = add_vectors(add_vectors(v[3], p3to4d), p4to1d)
        local m4 = add_vectors(add_vectors(v[4], p4to1d), p1to2d)

        -- points a third of the distance between corners
        local p1to2t = add_vectors(v[1], p1to2d)
        local p2to1t = add_vectors(v[2], p3to4d)

        local p2to3t = add_vectors(v[2], p2to3d)
        local p3to2t = add_vectors(v[3], p4to1d)

        local p3to4t = add_vectors(v[3], p3to4d)
        local p4to3t = add_vectors(v[4], p1to2d)

        local p4to1t = add_vectors(v[4], p4to1d)
        local p1to4t = add_vectors(v[1], p2to3d)

        local shape = create_shape_from_vertices({ m1, m2, m3, m4 }, model)
        local obj = ipe.Path(model.attributes, { shape })
        -- obj:set("fill", "black")
        -- obj:set("pathmode", "filled")
        model:creation(1, obj)
        
        sierpinski_carpet({v[1], p1to2t, m1, p1to4t}, iterations - 1, model)
        sierpinski_carpet({p1to2t, p2to1t, m2, m1}, iterations - 1, model)
        sierpinski_carpet({p2to1t, v[2], p2to3t, m2}, iterations - 1, model)
        sierpinski_carpet({m2, p2to3t, p3to2t, m3}, iterations - 1, model)
        sierpinski_carpet({m3, p3to2t, v[3], p3to4t}, iterations - 1, model)
        sierpinski_carpet({m4, m3, p3to4t, p4to3t}, iterations - 1, model)
        sierpinski_carpet({p4to1t, m4, p4to3t, v[4]}, iterations - 1, model)
        sierpinski_carpet({p1to4t, m1, m4, p4to1t}, iterations - 1, model)
    end
end

function sierpinski_carpet_run(model)
    local depth = tonumber(model:getString("Enter Depth"))
    local v = get_selected_poly_vertices(model)
    local vu = unique_points(v, model)

    if vu == nil then
        incorrect("waiter! waiter! more vertices, please!", model)
        return
    end

    if not is_square(vu, model) then
        incorrect("waiter! waiter! i need a square!", model)
        return
    end

    sierpinski_carpet(vu, depth, model)
end


-- ---------------------------------------------------------------------------
-- Dragon Curve
-- ---------------------------------------------------------------------------

function dragon(v, iterations)
    if iterations > 0 then
        -- directions (to manipulate later)
        local d1to2 = sub_vectors(v[2], v[1])
        local d2to3 = sub_vectors(v[3], v[2])

        -- sin/cos of 45 degrees (they are the same for 45 degrees)
        local sqrt2o2 = (2 ^ (1 / 2)) / 2
        
        -- iteration on "left" side of curve (extending outward), 45 degree counterclockwise rotation
        local d1to2rot = ipe.Vector(d1to2.x * sqrt2o2 - d1to2.y * sqrt2o2, d1to2.x * sqrt2o2 + d1to2.y * sqrt2o2)
        -- making vector shorter to create proper right angle isosceles triangle edge
        local new_v1 = scale_vector(d1to2rot, sqrt2o2)

        -- new point for next iteration
        local btw1to2 = add_vectors(v[1], new_v1)
        local t1 = dragon({ v[1], btw1to2, v[2] }, iterations - 1)
        
        -- iteration on "top" side of curve (extending inward), 45 degree clockwise rotation
        local d2to3rot = ipe.Vector(d2to3.x * sqrt2o2 - d2to3.y * (-sqrt2o2), d2to3.x * (-sqrt2o2) + d2to3.y * sqrt2o2)
        -- making vector shorter to create proper right angle isosceles triangle edge
        local new_v2 = scale_vector(d2to3rot, sqrt2o2)

        -- new point for next iteration
        local btw2to3 = add_vectors(v[2], new_v2)
        local t2 = dragon({ v[2], btw2to3, v[3] }, iterations - 1)

        -- the vertices of the sides get merged to one big table.
        -- this is then drawn!
        return merge(t1, t2)
    else
        -- vertices at lowest level simply returned.
        -- the table of vertices then gets passed to ipe to
        -- be drawn.
        return v
    end
end

function dragon_run(model)
    local v = get_selected_poly_vertices(model)
    local vu = unique_points(v, model)

    if vu == nil then
        incorrect("waiter! waiter! more vertices, please!", model)
        return
    end

    if not is_L(vu) then
        incorrect("Make sure the L is right-angled and has equal sides.", model)
        return
    end

    local out = model:getString("Enter iterations. Anything above 15-17\nwill take a while and may slow your computer.\nYou need to delete the original L.")

    if string.match(out, "^%d+$") then
        local dr = dragon(vu, tonumber(out))

        local shape = create_shape_from_vertices(dr, model)
        local obj = ipe.Path(model.attributes, { shape })
        model:creation(1, obj)
    else
        incorrect("waiter! waiter! i need a number!", model)
        return
    end
end

-- ---------------------------------------------------------------------------
-- Methods
-- ---------------------------------------------------------------------------

methods = {
  { label = "Sierpinski triangle", run = sierpinski_triangle_run },
  { label = "Sierpinski carpet",   run = sierpinski_carpet_run   },
  { label = "Dragon curve",        run = dragon_run              },
}
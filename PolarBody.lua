----------------------------------------------------------------------
-- PolarBody ipelet
----------------------------------------------------------------------
--[[

	Computes the polar body of a convex polygon with respect to a circle
	whose center lies inside it. Select the polygon and the circle, in
	either order.

	This file is intended to be used with the extensible drawing editor Ipe.
	Copyright (c) 2026 Auguste Gezalyan, Ryan Parker, and David M. Mount

	This software is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with Ipe; if not, you can find it at
			"http://www.gnu.org/copyleft/gpl.html",
	or write to the Free Software Foundation, Inc., 675 Mass Ave,
	Cambridge, MA 02139, USA.

--]]

label = "Polar body"

about = [[
Computes the polar body of a convex polygon relative to a circle whose
center lies inside the polygon. Select the polygon and the circle, in
either order. Each edge of the polygon becomes a vertex of the polar body.
]]

do
    local TOL = 1e-9

    local function incorrect(model, message)
        model:warning(message)
    end

    -- Positive, zero or negative as p, q, r turn counterclockwise, are
    -- collinear, or turn clockwise.
    local function orient(p, q, r)
        return p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y)
    end

    local function wrap_index(i, n)
        return ((i - 1) % n + n) % n + 1
    end

    -- Remove repeated vertices and vertices in the middle of a straight edge.
    local function clean_vertices(vertices)
        local v = {}
        for _, p in ipairs(vertices) do
            if #v == 0 or (p - v[#v]):len() > TOL then table.insert(v, p) end
        end
        if #v > 1 and (v[1] - v[#v]):len() <= TOL then table.remove(v) end
        local changed = true
        while changed and #v >= 3 do
            changed = false
            for i = 1, #v do
                local p, q, r = v[wrap_index(i - 1, #v)], v[i], v[wrap_index(i + 1, #v)]
                if math.abs(orient(p, q, r)) <= TOL * (q - p):len() * (r - q):len() then
                    table.remove(v, i)
                    changed = true
                    break
                end
            end
        end
        return v
    end

    -- True if every turn of the closed chain goes the same way.
    local function is_convex(vertices)
        local sign = 0
        for i = 1, #vertices do
            local o = orient(vertices[i], vertices[wrap_index(i + 1, #vertices)],
                             vertices[wrap_index(i + 2, #vertices)])
            if o ~= 0 then
                local s = o > 0 and 1 or -1
                if sign == 0 then sign = s elseif s ~= sign then return false end
            end
        end
        return true
    end

    -- True if the center is strictly inside, i.e. every edge sees it on the
    -- same side.
    local function contains_center(vertices, center)
        local this_orient = 1
        if orient(vertices[1], vertices[2], center) < 0 then
            this_orient = -1
        end
        for i = 1, #vertices do
            if this_orient * orient(vertices[i], vertices[wrap_index(i + 1, #vertices)], center) <= 0 then
                return false
            end
        end
        return true
    end

    -- Dual vertex of each edge: the point y with (y - c).(x - c) = r^2 for
    -- both endpoints x of the edge. Returns nil if it cannot be computed.
    local function get_dual_vertices(vertices, center, radius)
        local dual_verts = {}
        local rad2 = radius * radius
        for i, vertex in ipairs(vertices) do
            local i2 = wrap_index(i + 1, #vertices)
            local a = (vertex - center) * (1 / rad2)
            local b = (vertices[i2] - center) * (1 / rad2)
            local delta = b - a
            local denom = a.x * delta.y - a.y * delta.x
            if denom == 0 then
                return nil
            end
            table.insert(dual_verts, ipe.Vector(delta.y / denom, -delta.x / denom) + center)
        end
        return dual_verts
    end

    -- Vertices of a closed straight-edged polygon in page coordinates, or nil.
    local function polygon_from(obj)
        local shape = obj:shape()
        if #shape ~= 1 or shape[1].type ~= "curve" or not shape[1].closed or #shape[1] < 2 then
            return nil
        end
        local m = obj:matrix()
        local vertices = {}
        local subpath = shape[1]
        for i = 1, #subpath do
            if subpath[i].type ~= "segment" then return nil end
            table.insert(vertices, m * subpath[i][1])
        end
        table.insert(vertices, m * subpath[#subpath][2])
        return vertices
    end

    -- Center and radius of a circle, or nil (ellipses are rejected).
    local function circle_from(obj)
        local shape = obj:shape()
        if #shape ~= 1 or shape[1].type ~= "ellipse" then return nil end
        local m = obj:matrix()
        local ellipse = shape[1][1] -- maps the unit circle onto the ellipse
        local center = m * ellipse:translation()
        local r1 = (m * (ellipse * ipe.Vector(1, 0)) - center):len()
        local r2 = (m * (ellipse * ipe.Vector(0, 1)) - center):len()
        if math.abs(r1 - r2) > 1e-6 * math.max(r1, r2) then return nil end
        return center, r1
    end

    -- Exactly two selected objects: a closed polygon and a circle, any order.
    local function get_inputs(model)
        local p = model:page()
        local objs = {}
        for _, obj, sel, _ in p:objects() do
            if sel then table.insert(objs, obj) end
        end
        if #objs ~= 2 then
            incorrect(model, "Select a convex polygon and a circle")
            return
        end
        local vertices, center, radius
        for _, obj in ipairs(objs) do
            if obj:type() == "path" then
                local c, r = circle_from(obj)
                if c then
                    center, radius = c, r
                else
                    vertices = polygon_from(obj) or vertices
                end
            end
        end
        if not vertices or not center then
            incorrect(model, "Select a convex polygon (closed, straight edges) and a circle")
            return
        end
        return vertices, center, radius
    end

    local function create_ipe_polygon(vertices, model)
        local shape = { type = "curve", closed = true }
        for i = 1, #vertices - 1 do
            table.insert(shape, { type = "segment", vertices[i], vertices[i + 1] })
        end
        model:creation("Polar body", ipe.Path(model.attributes, { shape }))
    end

    function run(model)
        local vertices, center, radius = get_inputs(model)
        if not vertices then return end

        vertices = clean_vertices(vertices)
        if #vertices < 3 then
            incorrect(model, "The polygon needs at least three corners")
            return
        end
        if not is_convex(vertices) then
            incorrect(model, "The polygon must be convex")
            return
        end
        if not contains_center(vertices, center) then
            incorrect(model, "The circle's center must lie inside the polygon")
            return
        end

        local dual_verts = get_dual_vertices(vertices, center, radius)
        if not dual_verts then
            incorrect(model, "Could not compute the polar body")
            return
        end
        create_ipe_polygon(dual_verts, model)
    end


end

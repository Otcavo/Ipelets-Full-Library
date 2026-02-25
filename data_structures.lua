label = "Data structures"

about = [[
Data structure Ipelets:
- Quad tree
- Trapezoidal map
- Randomized trapezoidal map
- Triangulate
]]

-- ---------------------------------------------------------------------------
-- Helpers 
-- ---------------------------------------------------------------------------

function incorrect(title, model) model:warning(title) end

function display(title, message, model) 
    local s = title
    local d = ipeui.Dialog(model.ui:win(), "Output")
    d:add("label1", "label", {label=s}, 1, 1, 1, 2)
    d:add("input", "input", {}, 2, 1, 1, 2)
    d:addButton("ok", "&Ok", "accept")
    d:setStretch("column", 2, 1)
    d:setStretch("column", 1, 1)
    d:set("input", message)
    d:execute()
end

function dump(o)
    if _G.type(o) == 'table' then
       local s = '{ '
       for k,v in pairs(o) do
          if _G.type(k) ~= 'number' then k = '"'..k..'"' end
          s = s .. '['..k..'] = ' .. dump(v) .. ','
         do return s .. '} ' end
       end
       
    else
       return tostring(o)
    end
 end


 function get_polygon_segments(obj, model)

	local shape = obj:shape()
	local transform = obj:matrix()

	local segment_matrix = shape[1]

	local segments = {}
	for _, segment in ipairs(segment_matrix) do
		table.insert(segments, ipe.Segment(transform * segment[1], transform * segment[2]))
	end
	 
	table.insert(
		segments,
		ipe.Segment(transform * segment_matrix[#segment_matrix][2], transform * segment_matrix[1][1])
	)

	return segments
end

function in_box(bbox, startPoint, endPoint)
    local x_min, x_max, y_min, y_max = bbox[1], bbox[2], bbox[3], bbox[4]
    if ((x_min <= endPoint.x and endPoint.x <= x_max) and
        (x_min <= startPoint.x and startPoint.x <= x_max) and
        (y_min <= endPoint.y and endPoint.y <= y_max) and
        (y_min <= startPoint.y and startPoint.y <= y_max)) then
            return true
        else
            return false
        end
end

function get_pt_and_polygon_selection(model)
	local p = model:page()

    if not p:hasSelection() then
        incorrect(dump("Nothing Selected"), model)
        return
	end

	local count = 0

    local path_objects = {}

	for _, obj, sel, _ in p:objects() do
        if sel then
            count = count + 1
            if obj:type() == "path" then 
                table.insert(path_objects, obj)
            end
        end
	end

    local segments_table = {}
    local bounding_box = {}

    for i = 1, #path_objects do
        local segments = get_polygon_segments(path_objects[i], model)
        if #segments == 4 then
            table.insert(bounding_box, segments)
        else
            table.insert(segments_table, segments)
        end
    end

    -- Store the points of selected Line Segments and Calculate
    -- the max and min (extremities) of all coordinates

    local output_table = {{}, {}, {}}

    local x_min, y_min, x_max, y_max = math.huge, math.huge, -1 * math.huge, -1 * math.huge


    if #bounding_box ~= 0 then
        x_min = bounding_box[1][1]:endpoints().x 
        y_min = bounding_box[1][2]:endpoints().y 
        x_max = bounding_box[1][3]:endpoints().x 
        y_max = bounding_box[1][1]:endpoints().y 
    end


    for i = 1, #segments_table do
        local startPoint, endPoint = segments_table[i][1]:endpoints()
        if #bounding_box == 0 then
            x_min = math.min(endPoint.x, startPoint.x, x_min)
            y_min = math.min(endPoint.y, startPoint.y, y_min)
            x_max = math.max(endPoint.x, startPoint.x, x_max)
            y_max = math.max(endPoint.y, startPoint.y, y_max)
        end


        if (startPoint.x > endPoint.x) then -- Ensures Line Segments are organized left to right
            if (#bounding_box == 0 or in_box({x_min, x_max, y_min, y_max}, startPoint, endPoint)) then
                table.insert(output_table[1], {{endPoint.x, startPoint.x}, {endPoint.y, startPoint.y}})
            else
                display("The Following Segment was Ignored - Please Ensure it is fully contained in the bounding box", 
                        tableToString({{endPoint.x, startPoint.x}, {endPoint.y, startPoint.y}}), model)
            end
        else
            if (#bounding_box == 0 or in_box({x_min, x_max, y_min, y_max}, startPoint, endPoint)) then
                table.insert(output_table[1], {{startPoint.x, endPoint.x}, {startPoint.y, endPoint.y}})
            else
                display("The Following Segment was Ignored - Please Ensure it is fully contained in the bounding box", 
                        tableToString({{startPoint.x, endPoint.x}, {startPoint.y, endPoint.y}}), model)
            end
        end

    end

    if  #output_table[1] ~= #segments_table then
        incorrect(dump("Some Points were Ignored - Please draw the Bounding Box after modifying (Translating, Shearing...) segments"), model)
    end

    local scale = 20

    if #bounding_box == 0 then -- bounding box not given
        table.insert(output_table[2], {{x_min - scale, x_max + scale}, {y_min - scale, y_max + scale}})
        table.insert(output_table[3], false)
    else -- bounding box given
        table.insert(output_table[2], {{x_min, x_max}, {y_min, y_max}})
        table.insert(output_table[3], true)
    end
        
	return output_table
end

function create_boundary(x_min, x_max, y_min, y_max, scale, model)

    -- Draws a Boundary around the highlighted line sections

    local start = ipe.Vector(x_min - scale , y_min - scale) 
    local finish = ipe.Vector(x_max + scale, y_min - scale) 

    local segment = {type="segment", start, finish}
    local shape = { type="curve", closed=false, segment}
    local pathObj = ipe.Path(model.attributes, { shape })

    model:creation("create basic path", pathObj) 

    local start = ipe.Vector(x_min - scale, y_min - scale) 
    local finish = ipe.Vector(x_min - scale, y_max + scale)

    local segment = {type="segment", start, finish}
    local shape = { type="curve", closed=false, segment}
    local pathObj = ipe.Path(model.attributes, { shape })

    model:creation("create basic path", pathObj) 

    local start = ipe.Vector(x_min -  scale, y_max + scale) 
    local finish = ipe.Vector(x_max + scale, y_max + scale) 

    local segment = {type="segment", start, finish}
    local shape = { type="curve", closed=false, segment}
    local pathObj = ipe.Path(model.attributes, { shape })

    model:creation("create basic path", pathObj) 

    local start = ipe.Vector(x_max + scale, y_min - scale) 
    local finish = ipe.Vector(x_max + scale, y_max + scale)

    local segment = {type="segment", start, finish}
    local shape = { type="curve", closed=false, segment}
    local pathObj = ipe.Path(model.attributes, { shape })

    model:creation("create basic path", pathObj) 

end

function tableToString(tbl, indent)
    indent = ""
    local str = "{"
    for i, v in ipairs(tbl) do
        if _G.type(v) == 'table' then
            str = str .. tableToString(v, indent .. "  ")
        else
            str = str .. tostring(v)
        end
        if i < #tbl then
            str = str .. ", "
        end
    end
    str = str .. indent .. "}"
    return str
end

local PointRegionQuadtreeNode = {}

function PointRegionQuadtreeNode.new(boundary, capacity)
    local instance = {
        boundary = boundary,
        capacity = capacity,
        points = {},
        divided = false,
    }

    -- instance methodsx
    instance.subdivide = PointRegionQuadtreeNode.subdivide
    instance.insert = PointRegionQuadtreeNode.insert
    instance.belongs = PointRegionQuadtreeNode.belongs
    instance.draw = PointRegionQuadtreeNode.draw
    instance.to_string = PointRegionQuadtreeNode.to_string

    return instance
end

-- returns string representation of entire point-region quadtree
function PointRegionQuadtreeNode.to_string(self, depth)
    local indent = string.rep("         ", depth)
    local result = ""

    if #self.points > 0 then
        result = result .. "Points: {"
        for i, point in ipairs(self.points) do
            result = result .. "(" .. point.x .. ", " .. point.y .. ")"
            if i < #self.points then
                result = result .. ", "
            end
        end
        result = result .. "}\n"
    elseif depth ~= 0 then
        result = result .. "\n"
    end

    -- Recursively print child nodes if subdivided
    if self.divided then
        result = result .. indent .. "Northwest: " .. self.northwest:to_string(depth + 1)
        result = result .. indent .. "Northeast: " .. self.northeast:to_string(depth + 1)
        result = result .. indent .. "Southwest: " .. self.southwest:to_string(depth + 1)
        result = result .. indent .. "Southeast: " .. self.southeast:to_string(depth + 1)
    end

    return result
end

-- subdivision of node in equal quadrants when its number of points exceeds capacity
function PointRegionQuadtreeNode.subdivide(self)
    local midX = (self.boundary.min_x + self.boundary.max_x) / 2
    local midY = (self.boundary.min_y + self.boundary.max_y) / 2

    self.northwest = PointRegionQuadtreeNode.new({min_x = self.boundary.min_x, min_y = midY, max_x = midX, max_y = self.boundary.max_y}, self.capacity)
    self.northeast = PointRegionQuadtreeNode.new({min_x = midX, min_y = midY, max_x = self.boundary.max_x, max_y = self.boundary.max_y}, self.capacity)
    self.southwest = PointRegionQuadtreeNode.new({min_x = self.boundary.min_x, min_y = self.boundary.min_y, max_x = midX, max_y = midY}, self.capacity)
    self.southeast = PointRegionQuadtreeNode.new({min_x = midX, min_y = self.boundary.min_y, max_x = self.boundary.max_x, max_y = midY}, self.capacity)

    self.divided = true
end


-- recursively inserts new point into quadtree
function PointRegionQuadtreeNode.insert(self, point, root_boundary)
    
    -- checks if the point even belongs to this node's range 
    if not self:belongs(point, root_boundary) then
        return false
    end

    -- base case of recursive calls: point simply gets added to node
    if #self.points < self.capacity and not self.divided then
        table.insert(self.points, point)
        return true
    
    -- recursively calls insert to insert point into appropriate node (which may or may need subdivisions)
    else

        if not self.divided then
            self:subdivide()

            -- re-inserts each point from the subdivided node to the appropriate quadrant
            for i = #self.points, 1, -1 do

                local removed_point = table.remove(self.points, i)

                if self.northwest:belongs(removed_point, root_boundary) then
                    self.northwest:insert(removed_point, root_boundary)
                elseif self.northeast:belongs(removed_point, root_boundary) then
                    self.northeast:insert(removed_point, root_boundary)
                elseif self.southwest:belongs(removed_point, root_boundary) then
                    self.southwest:insert(removed_point, root_boundary)
                else
                    self.southeast:insert(removed_point, root_boundary)
                end
            end
        end

        if self.northwest:belongs(point, root_boundary) then
            return self.northwest:insert(point, root_boundary)
        elseif self.northeast:belongs(point, root_boundary) then
            return self.northeast:insert(point, root_boundary)
        elseif self.southwest:belongs(point, root_boundary) then
            return self.southwest:insert(point, root_boundary)
        else
            return self.southeast:insert(point, root_boundary)
        end
    end
end


-- returns whether point belongs within the boundary of the current object 
function PointRegionQuadtreeNode.belongs(self, point, root_boundary)
    return point.x >= self.boundary.min_x 
    and point.y >= self.boundary.min_y
    and (point.y < self.boundary.max_y or (point.y == self.boundary.max_y and point.y == root_boundary.max_y)) 
    and (point.x < self.boundary.max_x or (point.x == self.boundary.max_x and point.x == root_boundary.max_x)) 
end


-- recursively draws entire point-region quadtree when called on root
function PointRegionQuadtreeNode.draw(self, model)
    local min_x, min_y, max_x, max_y = self.boundary.min_x, self.boundary.min_y, self.boundary.max_x, self.boundary.max_y

    local box = ipe.Path(model.attributes, {{
        type = "curve",
        closed = true,
        {type = "segment", ipe.Vector(min_x, min_y), ipe.Vector(max_x, min_y)},
        {type = "segment", ipe.Vector(max_x, min_y), ipe.Vector(max_x, max_y)},
        {type = "segment", ipe.Vector(max_x, max_y), ipe.Vector(min_x, max_y)},
        {type = "segment", ipe.Vector(min_x, max_y), ipe.Vector(min_x, min_y)}
    }})

    model:creation("Box around node", box)

    if self.divided then
        self.northwest:draw(model)
        self.northeast:draw(model)
        self.southwest:draw(model)
        self.southeast:draw(model)
    end
end

-- ---------------------------------------------------------------------------
-- Quad tree 
-- ---------------------------------------------------------------------------

-- point quadtree class
local PointQuadtreeNode = {}

function PointQuadtreeNode.new(point, boundary)
    local instance = {
        point = point,
        boundary = boundary,
        northwest = nil,
        northeast = nil,
        southwest = nil,
        southeast = nil
    }

    -- instance methods
    instance.insert = PointQuadtreeNode.insert
    instance.belongs = PointQuadtreeNode.belongs
    instance.draw = PointQuadtreeNode.draw
    instance.to_string = PointQuadtreeNode.to_string

    return instance
end

-- returns string representation of entire point quadtree
function PointQuadtreeNode.to_string(self, depth)
    local indent = string.rep("         ", depth)
    local result = ""

    if depth == 0 then 
        result = "Root Node: " 
    end

    result = result .. "(" .. self.point.x .. ", " .. self.point.y .. ")\n"
    indent = string.rep("         ", depth + 1)

    result = result .. indent .. "Northwest: "
    if self.northwest then 
        result = result .. self.northwest:to_string(depth + 1) 
    else
        result = result .. "Null\n"
    end
    
    result = result .. indent .. "Northeast: "
    if self.northeast then 
        result = result .. self.northeast:to_string(depth + 1) 
    else
        result = result .. "Null\n"
    end
    
    result = result .. indent .. "Southwest: "
    if self.southwest then 
        result = result .. self.southwest:to_string(depth + 1)  
    else
        result = result .. "Null\n"
    end
    
    result = result .. indent .. "Southeast: "
    if self.southeast then 
        result = result .. self.southeast:to_string(depth + 1) 
    else
        result = result .. "Null\n"
    end

    return result
end

-- recursively inserts new point into quadtree
function PointQuadtreeNode.insert(self, point, root_boundary)

    -- checks if the point even belongs to this node's range 
    if not self:belongs(point, root_boundary) then
        return false
    end

    -- base case of recursive calls: point simply gets added to node
    if not self.point then
        self.point = point
        return true
        
    -- recursively calls insert to insert point into appropriate node (which may or may need subdivisions)
    else

        -- subdivision necessary
        if point.x < self.point.x then
            if point.y < self.point.y then
                if not self.southwest then
                    self.southwest = PointQuadtreeNode.new(nil, {min_x = self.boundary.min_x, min_y = self.boundary.min_y, max_x = self.point.x, max_y = self.point.y})
                end
                return self.southwest:insert(point, root_boundary)
            else
                if not self.northwest then
                    self.northwest = PointQuadtreeNode.new(nil, {min_x = self.boundary.min_x, min_y = self.point.y, max_x = self.point.x, max_y = self.boundary.max_y})
                end
                return self.northwest:insert(point, root_boundary)
            end
        else
            if point.y < self.point.y then
                if not self.southeast then
                    self.southeast = PointQuadtreeNode.new(nil, {min_x = self.point.x, min_y = self.boundary.min_y, max_x = self.boundary.max_x, max_y = self.point.y})
                end
                return self.southeast:insert(point, root_boundary)
            else
                if not self.northeast then
                    self.northeast = PointQuadtreeNode.new(nil, {min_x = self.point.x, min_y = self.point.y, max_x = self.boundary.max_x, max_y = self.boundary.max_y})
                end
                return self.northeast:insert(point, root_boundary)
            end
        end
    end
end


-- returns whether point belongs within the boundary of the current object 
function PointQuadtreeNode.belongs(self, point, root_boundary)
    return point.x >= self.boundary.min_x 
    and point.y >= self.boundary.min_y
    and (point.y < self.boundary.max_y or (point.y == self.boundary.max_y and point.y == root_boundary.max_y)) 
    and (point.x < self.boundary.max_x or (point.x == self.boundary.max_x and point.x == root_boundary.max_x)) 
end


-- recursively draws entire point quadtree when called on root
function PointQuadtreeNode.draw(self, model)
    if not self.point then return end

    local horizontal_line = ipe.Path(model.attributes, {{
        type = "curve",
        closed = false,
        {type = "segment", ipe.Vector(self.boundary.min_x, self.point.y), ipe.Vector(self.boundary.max_x, self.point.y)}
    }})
    model:creation("Horizontal line through point", horizontal_line)

    local vertical_line = ipe.Path(model.attributes, {{
        type = "curve",
        closed = false,
        {type = "segment", ipe.Vector(self.point.x, self.boundary.min_y), ipe.Vector(self.point.x, self.boundary.max_y)}
    }})
    model:creation("Vertical line through point", vertical_line)

    if self.northwest then self.northwest:draw(model) end
    if self.northeast then self.northeast:draw(model) end
    if self.southwest then self.southwest:draw(model) end
    if self.southeast then self.southeast:draw(model) end
end




local function get_unique_selected_points(model)
    local page = model:page()
    local points = {}

    -- goes through the selected objects on the page and adds the points to the points table
    for i, obj, sel, _ in page:objects() do
        if sel then
            if obj:type() == "reference" then
                local point = obj:position()
                local transform = obj:matrix()

                for _, existing_point in ipairs(points) do

                    while existing_point.x == point.x do
                        if math.random() >= 0.5 then
                            point = ipe.Vector(point.x + 0.001, point.y)
                        else
                            point = ipe.Vector(point.x - 0.001, point.y)
                        end
                    end

                    while existing_point.y == point.y do
                        if math.random() >= 0.5 then
                            point = ipe.Vector(point.x, point.y + 0.001)
                        else
                            point = ipe.Vector(point.x, point.y - 0.001)
                        end
                    end
                end

                local dx, dy = point.x - obj:position().x, point.y - obj:position().y
                local new_matrix = obj:matrix() * ipe.Matrix(1, 0, 0, 1, dx, dy)
                obj:setMatrix(new_matrix)

                table.insert(points, transform * point)
            end
        end
    end

    return points
end



-- gets coordinates of top left and bottom right vertices of the box 
local function get_box_vertices(points)

    local min_x, min_y = math.huge, math.huge
    local max_x, max_y = -math.huge, -math.huge

    for _, point in ipairs(points) do
        if point.x < min_x then min_x = point.x end
        if point.x > max_x then max_x = point.x end
        if point.y < min_y then min_y = point.y end
        if point.y > max_y then max_y = point.y end
    end

    return {min_x = min_x, min_y = min_y, max_x = max_x, max_y = max_y}
end


-- draws the box on the page
local function draw_box(model, box)

    local box = ipe.Path(model.attributes, {{
        type = "curve",
        closed = true,
        {type = "segment", ipe.Vector(box.min_x, box.min_y), ipe.Vector(box.max_x, box.min_y)},
        {type = "segment", ipe.Vector(box.max_x, box.min_y), ipe.Vector(box.max_x, box.max_y)},
        {type = "segment", ipe.Vector(box.max_x, box.max_y), ipe.Vector(box.min_x, box.max_y)},
        {type = "segment", ipe.Vector(box.min_x, box.max_y), ipe.Vector(box.min_x, box.min_y)}
    }})
    
    model:creation("Box", box)
end


local function create_box_and_get_points(model)

    local points = get_unique_selected_points(model)

    -- no points were selected
    if #points < 1 then
        model:warning("Please select at least one point!")
        return
    end

    -- get coordinates of top left and bottom right vertices of the box 
    local box_vertices = get_box_vertices(points)

    -- draw the bounding box around selected points
    draw_box(model, box_vertices)

    return points
end

local function create_point_region_quadtree_run(model)

    -- ends if no points have been selected
    local unique_points = get_unique_selected_points(model)

    if #unique_points < 1 then
        model:warning("Please select at least one point!")
        return
    end

    -- getting max node capacity from user (has to be an integer greater than or equal to 1)

    local s = "Please enter an integer greater than or equal to 1.\nThis will be the maximum capacity of each node in the point-region quadtree."
    local d = ipeui.Dialog(model.ui:win(), "Input Validation")
    d:add("label1", "label", {label=s}, 1, 1, 1, 2)
    d:add("label2", "label", {label="Input:"}, 2, 1)
    d:add("input", "input", {}, 2, 2)
    d:addButton("ok", "&Ok", "accept")
    d:addButton("cancel", "&Cancel", "reject")
    d:setStretch("column", 2, 1)

    local num = -1
  
    while true do
      if not d:execute() then return end
      num = tonumber(d:get("input"))
      if num and num >= 1 and math.floor(num) == num then
        break
      else
        ipeui.messageBox(model.ui:win(), "warning", "Invalid Input", "Please enter an integer greater than or equal to 1!")
      end
    end

    -- below only runs if user has inputted valid max node capacity value

    if num and num >= 1 and math.floor(num) == num then

        -- gets all selected points and draws bounding box
        local points = create_box_and_get_points(model)

        local boundary = get_box_vertices(points)
        local quadtree = PointRegionQuadtreeNode.new(boundary, num)

        -- insert all points into the quadtree
        for _, point in ipairs(points) do
            quadtree:insert(point, boundary)
        end

        quadtree:draw(model)
        
        model:creation("", ipe.Text(model.attributes, "Right Click then Edit Text!", ipe.Vector(boundary.min_x, boundary.max_y + 50), 200))
        model:creation("", ipe.Text(model.attributes, quadtree:to_string(0), ipe.Vector(boundary.min_x, boundary.max_y + 25), 200))

    end
end


local function create_point_quadtree_run(model)

    -- gets all selected points and draws bounding box
    local points = create_box_and_get_points(model)

    if not points then return end

    local boundary = get_box_vertices(points)
    local quadtree = PointQuadtreeNode.new(nil, boundary)

    -- insert all points into the quadtree
    for _, point in ipairs(points) do
        quadtree:insert(point, boundary)
    end

    quadtree:draw(model)

    model:creation("", ipe.Text(model.attributes, "Right Click then Edit Text!", ipe.Vector(boundary.min_x, boundary.max_y + 50), 200))
    model:creation("", ipe.Text(model.attributes, quadtree:to_string(0), ipe.Vector(boundary.min_x, boundary.max_y + 25), 200))

end

-- ---------------------------------------------------------------------------
-- Trapezoidal map 
-- ---------------------------------------------------------------------------

function trapezoidal_run(model)
    local everything = get_pt_and_polygon_selection(model) 

    local inpt = everything[1]
    local bbox = everything[2][1]
    local given = everything[3][1]


    local x_min = math.max(0, bbox[1][1])
    local x_max = bbox[1][2]
    local y_min = math.max(0, bbox[2][1])
    local y_max = bbox[2][2]

    if given == false then create_boundary(x_min, x_max, y_min, y_max, 0, model) end

    local outp = {}
    for i = 1, #inpt do
         

        local a1, a2 = inpt[i][1][1], inpt[i][1][2]
        local b1, b2 = inpt[i][2][1], inpt[i][2][2]

        local arr_outp = {{a1, b1}, {a2, b2}}

        table.insert(arr_outp, {a1, y_min})
        table.insert(arr_outp, {a1, y_max})
        table.insert(arr_outp, {a2, y_min})
        table.insert(arr_outp, {a2, y_max})

        table.insert(outp, arr_outp)
    end

    
    for segment_index = 1, #outp do
        local segments = outp[segment_index]
        local left, right = segments[1], segments[2]
        for i = 3, #segments do
            for j = 1, #inpt do

                local a1, b1, a2, b2

                if (i == 3 or i == 4) then
                    a1, b1, a2, b2 = left[1], left[2], left[1], segments[i][2]
                else
                    a1, b1, a2, b2 = right[1], right[2], right[1], segments[i][2]   
                end

                local x1, x2 = inpt[j][1][1], inpt[j][1][2]
                local y1, y2 = inpt[j][2][1], inpt[j][2][2]
                
                if (not ((x2 == a1 and y2 == b1) or (a1 == x1 and b1 == y1)) and (x1 <= a1 and a1 <= x2)) then

                    local function f(x)
                        local m = (y2 - y1) / (x2 - x1)
                        local b = y1 - m * x1
                        return m * x + b
                    end
        
                    local func = f(a1)

                    if ((b2 <= func and func <= b1) or (b1 <= func and func <= b2)) then
                        segments[i][2] = func
                    end
                end
            end
        end
    end
    
    for segment_index = 1, #outp do
        local segments = outp[segment_index]
        local left, right = segments[1], segments[2]

        local a1, b1, a2, b2

        for i = 3, #segments do
            if (i == 3 or i == 4) then
                a1, b1, a2, b2 = left[1], left[2], left[1], segments[i][2]
            else
                a1, b1, a2, b2 = right[1], right[2], right[1], segments[i][2]
            end
        
            local start = ipe.Vector(a1,b1)
            local finish = ipe.Vector(a2,b2)

            local segment = {type="segment", start, finish}
            local shape = { type="curve", closed=false, segment}
            local pathObj = ipe.Path(model.attributes, { shape })
            pathObj:set("stroke", "red")

            model:creation("create basic path", pathObj)
        end
    end

    display("Array of Segments in the form [Left Endpoint, Right Endpoint, Lower Left Point, Upper Left Point, Lower Right Point, Upper Right Point]",tableToString(outp, ""), model)
end


-- ---------------------------------------------------------------------------
-- Triangulate
-- ---------------------------------------------------------------------------

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

function get_pt_and_polygon_selection_triangulate(model)
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
      model:creation("Triangle", ipe.Path(model.attributes, { Tri }))
      ll_remove(vertlist, j)
    end

    node = ll_next_loop(vertlist, node)
  end

  return indices
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

function triangulate_run(model)
  math.randomseed()
  local vertices = get_pt_and_polygon_selection_triangulate(model)
  
  if vertices then
    triangulate(reorient_ccw(vertices), model)
  end
end

methods = {
    { label = "Point-region quadtree",      run = create_point_region_quadtree_run },
    { label = "Point quadtree",             run = create_point_quadtree_run },
    { label = "Trapezoidal map (segments)", run = trapezoidal_run },
    { label = "Triangulate",                run = triangulate_run }
}
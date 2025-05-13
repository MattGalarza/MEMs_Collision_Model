using Plots
gr()

# Create a new plot with higher quality settings
plt = plot(
    aspect_ratio=:equal,
    xlim=(-180, 180),
    ylim=(-100, 100),
    legend=false,
    grid=false,
    showaxis=false,
    ticks=false,
    background=:white,
    framestyle=:none,
    size=(1200, 600),
    dpi=300
)

# Define trapezoid dimensions - FLIPPED orientation: wide at top, narrow at bottom
height = 120
width_top = 70        # Wider at top
width_bottom = 30     # Narrower at bottom
spacing = 100         # Spacing between electrodes

# Vertical offset for electrodes - SWAPPED positions
middle_offset = -20   # Move middle (blue) electrode UP to connect to top (blue) rectangle
side_offset = 20      # Move side (gray) electrodes DOWN to connect to bottom (gray) rectangle

# Define rectangular body dimensions
body_width = 250
body_height = 30
# Calculate exact positions to connect at edges
top_rect_y = -height/2 + middle_offset - body_height/2   # Top edge of middle trapezoid
bottom_rect_y = height/2 + side_offset + body_height/2   # Bottom edge of side trapezoids

# Function to draw a trapezoid
function draw_trapezoid(x_center, y_offset, width_top, width_bottom, height, color)
    xs = Float64[
        x_center - width_top/2,    # Top left
        x_center + width_top/2,    # Top right
        x_center + width_bottom/2, # Bottom right
        x_center - width_bottom/2, # Bottom left
        x_center - width_top/2     # Close the shape
    ]
    
    ys = Float64[
        -height/2 + y_offset,    # Top
        -height/2 + y_offset,    # Top
        height/2 + y_offset,     # Bottom
        height/2 + y_offset,     # Bottom
        -height/2 + y_offset     # Close the shape
    ]
    
    plot!(plt, xs, ys, seriestype=:shape, fillcolor=color, linecolor=:black, linewidth=2)
end

# Function to draw a rectangle
function draw_rectangle(x_center, y_center, width, height, color)
    xs = Float64[
        x_center - width/2,  # Left
        x_center + width/2,  # Right
        x_center + width/2,  # Right
        x_center - width/2,  # Left
        x_center - width/2   # Close the shape
    ]
    
    ys = Float64[
        y_center - height/2, # Bottom
        y_center - height/2, # Bottom
        y_center + height/2, # Top
        y_center + height/2, # Top
        y_center - height/2  # Close the shape
    ]
    
    plot!(plt, xs, ys, seriestype=:shape, fillcolor=color, linecolor=:black, linewidth=2)
end

# Draw rectangular bodies with SWAPPED colors
draw_rectangle(0, top_rect_y, body_width, body_height, :steelblue)     # Top body now BLUE
draw_rectangle(0, bottom_rect_y, body_width, body_height, :lightgray)  # Bottom body now GRAY

# Draw trapezoid electrodes with corrected offsets
draw_trapezoid(-spacing, side_offset, width_top, width_bottom, height, :lightgray)  # Left fixed
draw_trapezoid(0, middle_offset, width_top, width_bottom, height, :steelblue)       # Middle mobile
draw_trapezoid(spacing, side_offset, width_top, width_bottom, height, :lightgray)   # Right fixed


display(plt)



using CairoMakie

# Device dimensions (in microns for realistic scaling)
shuttle_length = 400
shuttle_width = 300
frame_width = 50
device_width = 1400
device_height = 1000
electrode_length = 450
electrode_base_width = 30
electrode_tip_width = 9
electrode_spacing = 25
num_electrodes_per_side = 7
spring_length = 250
spring_width = 20
spring_spacing = 60

# Create figure
fig = Figure(resolution = (1000, 800), backgroundcolor = :lightgrey)

ax = Axis(fig[1, 1], 
    title = "MEMS Device Diagram",
    xlimits = (-(device_width / 2 + 200), device_width / 2 + 200),
    ylimits = (-(device_height / 2 + 200), device_height / 2 + 200),
    aspect = 1,
    backgroundcolor = :white,
    xticks = [],
    yticks = [],
)

# Draw the main shuttle mass
shuttle_rect = FRect(-shuttle_length / 2, -shuttle_width / 2, shuttle_length, shuttle_width)
rect!(ax, shuttle_rect, color = :white, strokewidth = 2, strokecolor = :black)

# Draw the frame
frame_rect = FRect(-device_width / 2, -device_height / 2, device_width, device_height)
rect!(ax, frame_rect, color = :lightgrey, strokewidth = 2, strokecolor = :black)

# Draw the electrodes (top and bottom)
function draw_electrodes(x_offset, y_start, direction)
    for i in 0:num_electrodes_per_side - 1
        y_pos = y_start + i * (electrode_length + electrode_spacing)
        poly!(ax, [
            (x_offset, y_pos),
            (x_offset + electrode_tip_width, y_pos),
            (x_offset + electrode_base_width, y_pos + electrode_length),
            (x_offset, y_pos + electrode_length)
        ], color = :cadetblue, strokewidth = 2, strokecolor = :black)
    end
end

# Top electrodes
draw_electrodes(-device_width / 2, -device_height / 2 + frame_width + 50, 1)
draw_electrodes(device_width / 2 - electrode_base_width, -device_height / 2 + frame_width + 50, -1)

# Bottom electrodes
draw_electrodes(-device_width / 2, device_height / 2 - frame_width - 50 - electrode_length, 1)
draw_electrodes(device_width / 2 - electrode_base_width, device_height / 2 - frame_width - 50 - electrode_length, -1)

# Draw the springs (left and right sides)
function draw_springs(x_offset)
    for i in 0:3
        y_pos = -spring_length / 2 + i * (spring_length + spring_spacing)
        rect!(ax, FRect(x_offset, y_pos, spring_width, spring_length), color = :white, strokewidth = 2, strokecolor = :black)
    end
end

draw_springs(-device_width / 2 + frame_width)
draw_springs(device_width / 2 - frame_width - spring_width)

# Draw the X-Y axis indicator
arrow!(ax, Point2f0(0, 0), Point2f0(100, 0), linewidth = 2, color = :red)
arrow!(ax, Point2f0(0, 0), Point2f0(0, 100), linewidth = 2, color = :red)
text!(ax, "X", position = (120, 0), align = (:center, :center), fontsize = 14, color = :red)
text!(ax, "Y", position = (0, 120), align = (:center, :center), fontsize = 14, color = :red)

fig

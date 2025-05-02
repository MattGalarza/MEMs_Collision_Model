using GLMakie

# Animation parameters - scaled proportions
const ASPECT_RATIO = 16/9
const SHUTTLE_RATIO = 0.42  # SHUTTLE_WIDTH as proportion of total width
const SHUTTLE_HEIGHT_RATIO = 0.2  # SHUTTLE_HEIGHT as proportion of total width
const TRAPEZOID_HEIGHT_RATIO = 0.35  # TRAPEZOID_HEIGHT as proportion of total width
const TOP_WIDTH_RATIO = 0.1   # Wide at top
const BASE_WIDTH_RATIO = 0.03  # Narrow at base
const MAX_TILT_ANGLE = 25  # degrees

# Fixed electrode position
const FIXED_ELECTRODE_X_RATIO = 0.46  # proportion of total width
const FIXED_ELECTRODE_Y_RATIO = 0.73  # proportion of total height

# MODIFIED: Reduced displacement parameters
const SHUTTLE_START_X_RATIO = 0.43  # Starting position (reduced)
const SHUTTLE_MOVEMENT_AMOUNT = 0.08  # Movement amount (reduced)

function animate_mems(; width = 900)
    # Calculate dimensions based on width
    height = width / ASPECT_RATIO
    
    # Calculate actual dimensions based on ratios
    shuttle_width = width * SHUTTLE_RATIO
    shuttle_height = width * SHUTTLE_HEIGHT_RATIO
    trapezoid_height = width * TRAPEZOID_HEIGHT_RATIO
    top_width = width * TOP_WIDTH_RATIO
    base_width = width * BASE_WIDTH_RATIO
    
    # Fixed electrode position
    fixed_electrode_x = width * FIXED_ELECTRODE_X_RATIO
    fixed_electrode_y = height * FIXED_ELECTRODE_Y_RATIO
    
    # Create the figure and axis
    fig = Figure(resolution = (width, height))
    ax = Axis(fig[1, 1])
    hidedecorations!(ax)
    hidespines!(ax)
    limits!(ax, 0, width, 0, height)
    
    # Animation time - observable
    time = Observable(0.0)
    
    # Function to calculate positions based on animation time - MODIFIED for smaller displacement
    function calculate_positions(t)
        # Phase 1: Shuttle approaches from right (0 to 0.3)
        # Phase 2: Collision and deformation (0.3 to 0.7)
        # Phase 3: Return to original position (0.7 to 1.0)
        
        if t < 0.3
            # Approaching from right to left
            normalized_t = t / 0.3
            shuttle_x = width * SHUTTLE_START_X_RATIO - normalized_t * width * SHUTTLE_MOVEMENT_AMOUNT
            tilt_angle = 0.0  # No tilt before collision
        elseif t < 0.7
            # Collision and deformation (continue moving left)
            normalized_t = (t - 0.3) / 0.4
            shuttle_x = width * (SHUTTLE_START_X_RATIO - SHUTTLE_MOVEMENT_AMOUNT) - normalized_t * width * (SHUTTLE_MOVEMENT_AMOUNT/2)
            
            # Tilt angle based on collision distance - increases with collision depth
            tilt_angle = min(MAX_TILT_ANGLE, normalized_t * MAX_TILT_ANGLE)
        else
            # Return (move back to right)
            normalized_t = (t - 0.7) / 0.3
            shuttle_x = width * (SHUTTLE_START_X_RATIO - SHUTTLE_MOVEMENT_AMOUNT*1.5) + normalized_t * width * (SHUTTLE_MOVEMENT_AMOUNT*1.5)
            tilt_angle = 0.0  # Return to no tilt
        end
        
        return shuttle_x, tilt_angle
    end
    
    # Function to calculate trapezoid points
    function calculate_trapezoid(base_x, base_y, height, narrow_width, wide_width, 
                                is_pointing_up, tilt_angle=0.0, tilt_direction="right")
        # Points before rotation
        if is_pointing_up
            # Upward-facing trapezoid (narrow at base, wide at top)
            points = [
                Point2f(base_x - wide_width/2, base_y - height),   # top left
                Point2f(base_x + wide_width/2, base_y - height),   # top right
                Point2f(base_x + narrow_width/2, base_y),          # bottom right
                Point2f(base_x - narrow_width/2, base_y)           # bottom left
            ]
        else
            # Downward-facing trapezoid (narrow at top, wide at bottom)
            points = [
                Point2f(base_x - narrow_width/2, base_y),          # top left
                Point2f(base_x + narrow_width/2, base_y),          # top right
                Point2f(base_x + wide_width/2, base_y + height),   # bottom right
                Point2f(base_x - wide_width/2, base_y + height)    # bottom left
            ]
        end
        
        
        # Apply tilt rotation if needed
        if tilt_angle != 0
            # Determine center point for rotation (bottom center of trapezoid)
            center_x = base_x
            center_y = base_y
            
            # Convert tilt angle to radians (positive for left, negative for right)
            tilt_radians = (tilt_direction == "left" ? 1 : -1) * tilt_angle * π / 180
            
            # Rotation function around the center point
            function rotate(point)
                # Translate to origin
                translated_x = point[1] - center_x
                translated_y = point[2] - center_y
                
                # Rotate
                rotated_x = translated_x * cos(tilt_radians) - translated_y * sin(tilt_radians)
                rotated_y = translated_x * sin(tilt_radians) + translated_y * cos(tilt_radians)
                
                # Translate back
                return Point2f(rotated_x + center_x, rotated_y + center_y)
            end
            
            # Apply rotation to all points
            return map(rotate, points)
        else
            # No rotation needed
            return points
        end
    end
    
    # Function to update the scene for each time value
    function draw_scene(t)
        # Clear the axis for redrawing
        empty!(ax)
        
        # Get positions for current time
        shuttle_x, tilt_angle = calculate_positions(t)
        
        # Set y positions
        shuttle_y = height * 0.78
        
        # Calculate shuttle center
        shuttle_center_x = shuttle_x + shuttle_width/2
        
        # Mobile electrode base position aligned with x₁ dot
        mobile_electrode_y = shuttle_y
        mobile_electrode_base_x = shuttle_center_x
        
        # Calculate both trapezoids with corrected orientations
        # Fixed electrode: pointing up
        fixed_trapezoid = calculate_trapezoid(fixed_electrode_x, fixed_electrode_y, 
                                            trapezoid_height, base_width, top_width, true)
        
        # Mobile electrode: pointing up, with tilt to the left 
        mobile_trapezoid = calculate_trapezoid(mobile_electrode_base_x, mobile_electrode_y, 
                                            trapezoid_height, base_width, top_width, true, 
                                            tilt_angle, "left")
        
        # Text mode indicator
        mode_text = tilt_angle > 0 ? "Rotation" : "Translation"
        
        # Calculate x₁ and x₂ positions for the red dots
        x1_position = shuttle_center_x
        x2_position = (mobile_trapezoid[1][1] + mobile_trapezoid[2][1])/2
        
        # Calculate distance between electrodes for detection of contact
        distance = abs(fixed_electrode_x - mobile_electrode_base_x)
        is_approaching = distance > width * 0.019
        is_colliding = distance <= width * 0.019
        
        # Font and line sizes
        font_size_large = round(Int, width/25)
        font_size_medium = round(Int, width/40)
        font_size_small = round(Int, width/50)
        stroke_width = width/250
        marker_size = width/40
        
        # Mode indicator box - MOVED TO SIDE INSTEAD OF TOP
        box_width = width * 0.15 
        box_height = height * 0.08
        mode_x = width * 0.05
        mode_y = height * 0.5  # Moved to the middle of the left side
        
        rect_points = [
            Point2f(mode_x, mode_y),
            Point2f(mode_x + box_width, mode_y),
            Point2f(mode_x + box_width, mode_y - box_height),
            Point2f(mode_x, mode_y - box_height)
        ]
        poly!(ax, rect_points, color = :white, strokewidth = stroke_width*1.5, strokecolor = :black)
        text!(ax, mode_text, position = (mode_x + box_width/2, mode_y - box_height/2), 
            align = (:center, :center), fontsize = font_size_medium, color = :black, font = :italic)
        
        # Fixed trapezoid
        poly!(ax, fixed_trapezoid, color = :white, strokewidth = stroke_width, strokecolor = :black)
        
        # Mobile trapezoid
        poly!(ax, mobile_trapezoid, color = :white, strokewidth = stroke_width, strokecolor = :black)
        
        # Shuttle
        shuttle_points = [
            Point2f(shuttle_x, shuttle_y),
            Point2f(shuttle_x + shuttle_width, shuttle_y),
            Point2f(shuttle_x + shuttle_width, shuttle_y + shuttle_height),
            Point2f(shuttle_x, shuttle_y + shuttle_height)
        ]
        poly!(ax, shuttle_points, color = :white, strokewidth = stroke_width, strokecolor = :black)
        
        # SPRING REMOVED - Not drawing spring or anchor to remove from top left
        
        # Position indicators - x₁ and x₂ red dots
        scatter!(ax, [Point2f(x1_position, shuttle_y)], 
                color = :red, markersize = marker_size*1.5)
        text!(ax, "x₁(t)", position = (x1_position, shuttle_y + shuttle_height * 0.3), 
             align = (:center, :center), fontsize = font_size_medium)
        
        scatter!(ax, [Point2f(x2_position, mobile_trapezoid[1][2])], 
                color = :red, markersize = marker_size*1.5)
        text!(ax, "x₂(t)", position = (x2_position, mobile_trapezoid[1][2] - trapezoid_height * 0.1), 
             align = (:center, :center), fontsize = font_size_medium)
        
        # Reference vertical line
        lines!(ax, [Point2f(shuttle_center_x, height * 0.05), 
                   Point2f(shuttle_center_x, height)], 
               color = :red, linestyle = :dash, linewidth = stroke_width/3)
        
        # Force arrow during collision
        if is_colliding && tilt_angle > 0
            force_length = width * 0.06
            arrow_start = Point2f(mobile_trapezoid[1][1], mobile_trapezoid[1][2] + trapezoid_height * 0.13)
            arrow_end = Point2f(mobile_trapezoid[1][1] - force_length, mobile_trapezoid[1][2] + trapezoid_height * 0.13)
            arrow_vector = arrow_end - arrow_start
            
            # Draw arrow line
            lines!(ax, [arrow_start, arrow_end], color = :red, linewidth = stroke_width)
            
            # Draw arrowhead manually
            arrowhead_size = width * 0.015
            arrowhead_points = [
                Point2f(arrow_end[1] + arrowhead_size, arrow_end[2] - arrowhead_size/2),
                Point2f(arrow_end),
                Point2f(arrow_end[1] + arrowhead_size, arrow_end[2] + arrowhead_size/2)
            ]
            poly!(ax, arrowhead_points, color = :red, strokewidth = 0)
            
            text!(ax, "ke", position = (arrow_start[1] - force_length * 0.6, arrow_start[2] - trapezoid_height * 0.07), 
                 color = :red, fontsize = font_size_small)
        end
        
        # No Title - removed to simplify display
    end
    
    # Set up animation
    on(time) do t
        draw_scene(t)
    end
    
    # Initial drawing
    draw_scene(time[])
    
    # Animation loop
    @async begin
        while true
            # Update time value and loop animation
            new_time = time[] + 0.005
            if new_time > 1.0
                new_time = 0.0
            end
            time[] = new_time
            sleep(1/60)  # aim for 60fps
        end
    end
    
    return fig
end

# Run the animation
fig = animate_mems(width=900)
display(fig)
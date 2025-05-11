using Plots
using LaTeXStrings
# Function to draw a trapezoid - wider at bottom
function draw_trapezoid!(plt, x_center, width_top, width_bottom, height, color)
    # Calculate the left offset to center the trapezoid properly
    left_offset = (width_bottom - width_top) / 2
    
    # Calculate the four corners of the trapezoid
    top_left_x = x_center - width_top/2
    top_right_x = x_center + width_top/2
    bottom_right_x = x_center + width_bottom/2
    bottom_left_x = x_center - width_bottom/2
    
    top_y = -height/2
    bottom_y = height/2
    
    # Define the points in the correct order (closing the polygon)
    xs = [top_left_x, top_right_x, bottom_right_x, bottom_left_x, top_left_x]
    ys = [top_y, top_y, bottom_y, bottom_y, top_y]
    
    # Plot the filled trapezoid
    plot!(plt, xs, ys, seriestype=:shape, fillcolor=color, linecolor=:black, linewidth=2)
end

# Function to create the electrode drawing
function create_electrode_diagram()
    # Create a new plot with the desired settings
    plt = plot(
        aspect_ratio=:equal,
        xlim=(-150, 150),
        ylim=(-60, 60),
        legend=false,
        grid=false,
        showaxis=false,
        ticks=false,
        background=:white,
        size=(800, 400)
    )
    
    # Define electrode parameters
    height = 60          # Height of electrodes
    width_top = 30       # Width at top (narrower)
    width_bottom = 50    # Width at bottom (wider)
    spacing = 70         # Spacing between electrode centers
    
    # Draw the fixed electrodes (gray)
    draw_trapezoid!(plt, -spacing, width_top, width_bottom, height, :lightgray)
    draw_trapezoid!(plt, spacing, width_top, width_bottom, height, :lightgray)
    
    # Draw the mobile electrode (blue)
    draw_trapezoid!(plt, 0, width_top, width_bottom, height, :steelblue)
    
    return plt
end

# Create and display the electrode diagram
p = create_electrode_diagram()
display(p)


using Plots
using LaTeXStrings

function draw_trapezoid(x, y, width_top, width_bottom, height, color; filled=true)
    # For a trapezoid with the wider side at the top
    xs = [x, x + width_top, x + width_bottom, x, x]
    ys = [y + height, y + height, y, y, y + height]
    
    if filled
        plot!(Shape(xs, ys), fillcolor=color, linecolor=:black, linewidth=1.5)
    else
        plot!(xs, ys, linecolor=:black, linewidth=1.5)
    end
end

function draw_electrode_set()
    # Plot settings
    plot(size=(800, 600), 
         legend=false, 
         framestyle=:none,
         xlims=(-60, 60), 
         ylims=(-25, 25), 
         aspect_ratio=:equal,
         background=:white,
         margin=0,
         grid=false, 
         showaxis=false, 
         ticks=false)
    
    # Dimensions
    trap_height = 30
    trap_width_top = 20
    trap_width_bottom = 12
    gap = 5
    
    # Draw fixed electrodes (gray)
    draw_trapezoid(-50, -trap_height/2, trap_width_top, trap_width_bottom, trap_height, :lightgray)
    draw_trapezoid(30, -trap_height/2, trap_width_top, trap_width_bottom, trap_height, :lightgray)
    
    # Draw mobile electrode (blue)
    draw_trapezoid(-10, -trap_height/2, trap_width_top, trap_width_bottom, trap_height, :steelblue)
    
    return current()
end

# Create and save the plot
p = draw_electrode_set()
display(p)
savefig(p, "trapezoidal_electrodes.png")
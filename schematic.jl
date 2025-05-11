using Plots
using LaTeXStrings

function draw_trapezoid(x, y, width_top, width_bottom, height, color; filled=true)
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
    plot(size=(800, 600), legend=false, framestyle=:box,
         xlims=(-60, 60), ylims=(-40, 40), aspect_ratio=:equal,
         title="Trapezoidal Comb-Drive Electrodes",
         fontfamily="Computer Modern", dpi=300,
         grid=false, showaxis=false, ticks=false)
    
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
    
    # Add labels and dimensions
    annotate!(0, 25, text(L"L_f", 12, :black))
    annotate!(36, 25, text(L"C_{air\_r}", 12, :black))
    annotate!(-15, 25, text(L"W_f", 12, :black))
    annotate!(-30, 0, text(L"g_0 + z_s", 12, :black))
    annotate!(20, 0, text(L"g_0 - z_s", 12, :black))
    annotate!(6, -35, text(L"Direction of motion \rightarrow", 12, :black))
    
    # Add arrows for dimensions
    plot!([0, 0], [15, 35], arrow=true, color=:black, linewidth=1.5)
    plot!([-15, -15], [10, 25], arrow=true, color=:black, linewidth=1.5)
    plot!([-30, -15], [0, 0], arrow=true, color=:black, linewidth=1.5)
    plot!([5, 20], [0, 0], arrow=true, color=:black, linewidth=1.5)
    
    # Add capacitance blocks
    rectangle(x, y, w, h) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
    plot!(rectangle(-18, 5, 8, 6), fillcolor=:gold, linewidth=1.5, linecolor=:black)
    annotate!(-14, 8, text(L"C_{r}", 10, :black))
    
    plot!(rectangle(12, 5, 8, 6), fillcolor=:gold, linewidth=1.5, linecolor=:black)
    annotate!(16, 8, text(L"C_{l}", 10, :black))
    
    # Angle indicator
    θ = LinRange(π/4, π/2, 20)
    r = 10
    plot!([-10 .+ r*cos.(θ)], [-trap_height/2 .+ r*sin.(θ)], color=:black, linewidth=1.5)
    annotate!([-10 + r*cos(π/3)], [-trap_height/2 + r*sin(π/3)], text(L"\alpha", 12, :black))
    
    return current()
end

# Create and save the plot
p = draw_electrode_set()
display(p)
savefig(p, "trapezoidal_electrodes.png")
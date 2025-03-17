"""
    plotbranchingprocess(tree, p=plot())

Plot the simulated trajectories of a branching process stored in a tree of type `SimTree` on a plot `p`. For a multivariate process, set the index of the variable to plot with the keyword argument `col` (default value 1 for a univariate process); the index corresponds to the column of the `x` field of the tree.
"""
function plotbranchingprocess(tree; p=plot(), col=1)
    # set a default colorscheme based on the height (number of generations) of the tree
    ncolr = treeheight(tree)+1
    if ncolr<3
        ncolr=3
    elseif ncolr>11
        ncolr=11
    end
    
    colors = colorschemes[Symbol(@sprintf("RdYlBu_%d", ncolr))].colors
    plotsimtree!(p,tree,col,colors)
    # set the axis labels
    push!(p, Guide.xlabel("time"))
    push!(p, Guide.ylabel("x"))
    return p
end

"""
    plotsimtree!(p, tree, col, colors)

Recursively add the simulated trajectories of variable `col` in a branching process stored in a tree of type `SimTree` to the plot `p`, giving each generation a separate color from the array `colors`.
""" 
function plotsimtree!(p, tree::SimTree, col, colors)
    # each generation has a different color
    coli = min(treeheight(tree)+1, length(colors))
    color = colors[coli]
    # add the root node trajectory to the plot
    push!(p,layer(x=tree.t,y=tree.x[:,col],color=[color],Geom.line))
    for node in tree.children
        plotsimtree!(p,node,col,colors)
    end
end

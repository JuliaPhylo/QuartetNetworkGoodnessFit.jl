"""
    ultrametrize!(net::HybridNetwork, verbose::Bool)

Assign values to missing branch lengths in `net` to make the
network time-consistent (all paths from the root to a given hybrid node
have the same length) and ultrametric (all paths from the root to
the tips have the same length), if possible.

Warnings are given if it's not possible and if `verbose` is true.

The major tree is used to assign node ages. If an edge has a missing
length, this length is changed to the following:
- 0 if the edge is internal and non-minor,
- the value needed to make the network time-consistent if the edge is minor,
- the smallest value possible to make the network ultrametric
  if the edge is external.
"""
function ultrametrize!(net::HybridNetwork, verbose::Bool)
    directEdges!(net)
    preorder!(net) # creates / modified net.nodes_changed
    # first: calculate the distance from the root to each node, in dr
    dr = Vector{Float64}(undef, length(net.nodes_changed))
    dr[1] = 0.0 # root
    for i in 2:length(net.nodes_changed) # pre-order
        @inbounds n = net.nodes_changed[i]
        e = PhyloNetworks.getMajorParentEdge(n)
        p = PhyloNetworks.getParent(e)
        pi = findfirst(x -> x===p, net.nodes_changed)
        dr[i] = dr[pi] + max(e.length, 0.0) # interpret missing as zero
        if !n.leaf && e.length == -1.0
            e.length = 0.0
        end
        if n.hybrid
            emin = PhyloNetworks.getMinorParentEdge(n)
            pmin = PhyloNetworks.getParent(emin)
            pi = findfirst(x -> x===pmin, net.nodes_changed)
            if emin.length == -1.0
                emin.length = dr[i] - dr[pi]
                if emin.length < 0.0
                    verbose && @warn("could not make network time consistent at hybrid $(n.name)")
                    emin.length = 0.0
                end
            else
                dmin = dr[pi] + emin.length
                !verbose || isapprox(dmin, dr[i]) ||
                    @warn("hybrid $(n.name) at distance $(dr[i]) from the root via its major parent, but $dmin via its minor parent")
            end
        end
    end
    height = maximum(dr) # max distance from root to tip: will be tree height
    # adjust external branch lengths
    for i in 2:length(net.nodes_changed)
        @inbounds n = net.nodes_changed[i]
        n.leaf || continue
        e = PhyloNetworks.getMajorParentEdge(n)
        if e.length == -1.0
            e.length = height - dr[i]
        elseif !isapprox(dr[i], height)
            verbose && @warn("could not make network ultrametric at leaf $(n.name)")
        end
    end
    return net
end

"""
    reroot!(net, refnet)

Reroot `net` to minimize the hardwired cluster distance between
the `net` (with the new root position) and the reference network `refnet`.
Candidate root positions are limited to internal nodes (excluding leaves)
that are compatible with the direction of hybrid edges.
"""
function reroot!(net, refnet)
    root_saved = net.root # fall back in case none are admissible
    nnodes = length(net.node)
    bestdissimilarity = typemax(Int)
    besti = nothing
    for i in nnodes:-1:1
        net.node[i].leaf && continue
        net.root = i
        try
            directEdges!(net)
        catch e
            isa(e, PhyloNetworks.RootMismatch) || rethrow(e)
            continue
        end
        # now, i is admissible: internal and compatible with direction
        diss = hardwiredClusterDistance(net, refnet, true) # rooted = true now
        if diss < bestdissimilarity
            bestdissimilarity = diss
            besti = i
        end
        bestdissimilarity == 0 && break # cannot do better than 0!
    end
    if isnothing(besti) # all internal nodes conflicted with edge directions!
        net.root = root_saved
        try directEdges!(net); catch; end
    elseif bestdissimilarity > 0
        net.root = besti
        directEdges!(net)
    end
    return bestdissimilarity
end

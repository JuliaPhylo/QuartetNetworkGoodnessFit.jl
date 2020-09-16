"""
    ultrametrize!(net::HybridNetwork, verbose::Bool)

Assign values to missing branch lengths in `net` to make the
network time-consistent (all paths from the root to a given hybrid node
have the same length) and ultrametric (all paths from the root to
the tips have the same length), if possible.
Warnings are given if it's not possible and if `verbose` is true.

Output: true if the modified network is ultrametric, false otherwise.

The major tree is used to calculate the distance from nodes to the root.
If a tree edge has a missing length, this length is changed to the following:
- 0 if the edge is internal,
- the smallest value possible to make the network ultrametric
  if the edge is external.

It is assumed that hybrid nodes are not leaves, such that external edges
are necessarily tree edges.
If a hybrid edge has a missing length, this length is changed as follows:
- If both partner hybrid edges lack a length: the shortest lengths are assigned
  to make the network time-consistent at the hybrid node. In particular,
  either the major edge or the minor edge is assigned length 0.0.
- Otherwise: the value needed to make the network time-consistent considering
  based on the partner edge's length if this value is non-negative,
  and 0 if the ideal value is negative.

"""
function ultrametrize!(net::HybridNetwork, verbose::Bool)
    ultrametric = true
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
        drmaj = dr[pi]    # distance to root of major parent
        if !n.hybrid # n is a tree node. only 1 parent edge: e
            dr[i] = drmaj + max(e.length, 0.0) # interpret missing as zero
            if !n.leaf && e.length == -1.0
                e.length = 0.0 # leave external edge lengths unchanged
            end
            continue
        end
        # now: n is a hybrid node, with parent edges e and emin
        emin = PhyloNetworks.getMinorParentEdge(n)
        pmin = PhyloNetworks.getParent(emin)
        pmini = findfirst(x -> x===pmin, net.nodes_changed)
        drmin = dr[pmini] # distance to root of minor parent
        majmissing = e.length == -1.0
        minmissing = emin.length == -1.0
        if majmissing && minmissing # both parent edges lack a length
            if drmin > drmaj
                emin.length = 0.0
                e.length = drmin - drmaj
                dr[i] = drmin
            else
                e.length = 0.0
                emin.length = drmaj - drmin
                dr[i] = drmaj
            end
        elseif majmissing # major length missing, but minor parent has length
            dr[i] = drmin + emin.length # ideal situation
            e.length = dr[i] - drmaj
            if e.length < 0.0
                if !isapprox(dr[i], drmaj)
                  verbose && @warn("could not adjust length of major edge $(emin.number) to make network time consistent at hybrid $(n.name)")
                  ultrametric = false
                end
                e.length = 0.0
                dr[i] = drmaj
                # line above: to calculate dr[i] (distance to root) via major tree
                #   delete if instead we wanted to calculate dr[i]
                #   based on edge lengths that were originally non-missing.
            end
        elseif minmissing # the major parent has a length, the minor doesn't
            dr[i] = drmaj + e.length
            emin.length = dr[i] - drmin
            if emin.length < 0.0
                if !isapprox(dr[i], drmin)
                  verbose && @warn("could not adjust length of minor edge $(emin.number) to make network time consistent at hybrid $(n.name)")
                  ultrametric = false
                end
                emin.length = 0.0
            end
        else # both parent edges have a length
            dr[i] = drmaj + e.length
            dalt  = drmin + emin.length
            if !isapprox(dalt, dr[i])
              verbose && @warn("hybrid $(n.name) at distance $(dr[i]) from the root via its major parent, but $dalt via its minor parent")
              ultrametric = false
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
            ultrametric = false
        end
    end
    return ultrametric
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

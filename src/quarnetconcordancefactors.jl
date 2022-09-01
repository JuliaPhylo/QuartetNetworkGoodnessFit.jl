"""
fixit: use i below, index for the quartet to use
"""
function network_expectedCF(net::HybridNetwork; showprogressbar=true)
    net.node[net.root].leaf && error("The root can't be a leaf.")
    taxa = sort!(tipLabels(net))
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    nCk = PhyloNetworks.nchoose1234(ntax) # matrix to rank 4-taxon sets
    qtype = MVector{3,Float64} # 3 floats: CF12_34, CF13_24, CF14_23; initialized at 0.0
    numq = nCk[ntax+1,4]
    quartet = Vector{PhyloNetworks.QuartetT{qtype}}(undef, numq)
    ts = [1,2,3,4]
    for qi in 1:numq
        quartet[qi] = PhyloNetworks.QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.))
        # next: find the 4-taxon set with the next rank,
        #       faster than using the direct mapping function
        ind = findfirst(x -> x>1, diff(ts))
        if ind === nothing ind = 4; end
        ts[ind] += 1
        for j in 1:(ind-1)
            ts[j] = j
        end
    end
    if showprogressbar
        nstars = (numq < 50 ? numq : 50)
        nquarnets_perstar = (numq/nstars)
        println("Calculation quartet CFs for $numq quartets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end
    for qi in 1:numq
        network_expectedCF!(quartet[qi], net, taxa, taxonnumber)
        if showprogressbar && qi >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * nquarnets_perstar))
        end
    end
    showprogressbar && print("\n")
    return quartet, taxa
end
"""
fixit: use i below, index for the quartet to use
"""
function network_expectedCF!(quartet, net, taxa, taxonnumber)
    # delete all taxa except for the 4 in the quartet
    net = deepcopy(net)
    for taxon in taxa
        taxonnumber[taxon] in quartet.taxonnumber && continue
        deleteleaf!(net, taxon)
    end
    quartet.data .= network_expectedCF_4taxa!(net, taxa[quartet.taxonnumber])
end
# below: net is assumed to have 4 taxa only, quartet ordering to follow ordering of taxon names in `fourtaxa`
function network_expectedCF_4taxa!(net, fourtaxa)
    # find and delete degree-2 blobs along external edges
    bcc = biconnectedComponents(net, true) # true: ignore trivial blobs
    entry = PhyloNetworks.biconnectedcomponent_entrynodes(net, bcc)
    entryindex = indexin(entry, net.nodes_changed)
    exitnodes = PhyloNetworks.biconnectedcomponent_exitnodes(net, bcc, false) # don't redo the preordering
    bloborder = sortperm(entryindex) # pre-ordering for blobs in their own blob tree
    function isexternal(ib) # is bcc[ib] of degree 2 and adjacent to an external edge?
        # yes if: 1 single exit adjacent to a leaf
        length(exitnodes[ib]) != 1 && return false
        ch = PhyloNetworks.getChildren(exitnodes[ib][1])
        return length(ch) == 1 && ch[1].leaf
    end
    for ib in reverse(bloborder)
        isexternal(ib) || continue # keep bcc[ib] if not external of degree 2
        for he in bcc[ib]
            he.isMajor && continue
            # deletion of a hybrid can hide the deletion of another: check that he is still in net
            any(e -> e===he, net.edge) || continue
            PhyloNetworks.deletehybridedge!(net,he)
        end
    end
    deleteaboveLSA!(net)
    PhyloNetworks.removedegree2nodes!(net) # just to make sure
    ndes = 4 # number of taxa descendant from lowest hybrid node
    if net.numHybrids > 0
        preorder!(net)
        # find a lowest hybrid node and # of taxa below it
        lowesth = net.nodes_changed[findlast(n -> n.hybrid, net.nodes_changed)]
        childedge = PhyloNetworks.getChildEdge(lowesth)
        ndes = length(PhyloNetworks.descendants(childedge))
    end
    if ndes > 2 # simple formula for qCF: find cut edge and its length
        cutpool = nothing
        if net.numHybrids == 0
            cutpool = net.edge
        else
            n2 = PhyloNetworks.getChild(childedge)
            cutpool = [e for e in n2.edge if e !== childedge]
        end
        sistertofirst = 2    # arbitrarily correct if polytomy
        internallength = 0.0 # correct if polytomy
        for e in cutpool # skip external edges
            PhyloNetworks.getChild(e).leaf && continue
            internallength += e.length
            hwc = hardwiredCluster(e, fourtaxa)
            sistertofirst = findnext(x -> x == hwc[1], hwc, 2)
        end
        minorcf = exp(-internallength)/3
        majorcf = 1.0 - 2 * minorcf
        qCF = (sistertofirst == 2 ? (majorcf,minorcf,minorcf) :
              (sistertofirst == 3 ? (minorcf,majorcf,minorcf) :
                                    (minorcf,minorcf,majorcf)   ))
        return qCF
    elseif ndes == 1 # weighted qCFs average of the 2 displayed networks
        qCF = (0.0,0.0,0.0)
        parenthedge = [e for e in lowesth.edge if PhyloNetworks.getChild(e) === lowesth]
        all(h.hybrid for h in parenthedge) || error("hybrid $(lowesth.number) has a parent edge that's a tree edge")
        parenthindex = indexin(parenthedge, net.edge)
        nhe = length(parenthedge)
        for i in 1:nhe # keep parenthedge[i], remove all others
            gamma = parenthedge[i].gamma
            if i < nhe
                net = deepcopy(net) # to save memory allocation
            end
            for j in 1:nhe
                j == i && continue # don't delete hybrid edge i!
                PhyloNetworks.deletehybridedge!(net, net.edge[parenthindex[j]])
            end
            qCF .+= gamma .* network_expectedCF_4taxa!(net, fourtaxa)
        end
        return qCF
    elseif ndes == 2 # weighted qCFs average of 3 networks: 2 displayed, 1 "parental"
        # if ab below h: (1-e^-t) for ab|cd + e^-t * (γ1^2 CF(N11) + γ2^2 CF(N22) + 2γ1γ2 symmCF(N12))
        # to do
    end
    @show net
end

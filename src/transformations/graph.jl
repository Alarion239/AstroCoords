using StaticArrays
using LinearAlgebra: I

# —————————————————————————————————————————————————————————————————————————————
#   TABLES: store the graph's atomic edges, adjacency, and compiled paths
# —————————————————————————————————————————————————————————————————————————————

# Maps (SourceFrame,TargetFrame) → edge kind symbol
const edge_kind = IdDict{Tuple{Type,Type},Symbol}()

# Maps (SourceFrame,TargetFrame) → atomic payload
#   • For :static_matrix → SMatrix{3,3}
#   • For :static_affine → Tuple{SMatrix{3,3}, SVector{3}}
#   • For :dynamic_* & :callable → a Callable
const edge_payload = IdDict{Tuple{Type,Type},Any}()

# Maps (SourceFrame,TargetFrame) → priority (lower = higher priority)
const edge_priority = IdDict{Tuple{Type,Type},Float64}()

# Adjacency list: for each frame type, which neighbors it points to
const adjacency_list = IdDict{Type,Vector{Type}}()

# Tracks which (Source,Target) methods have already been emitted
const emitted_paths = Set{Tuple{Type,Type}}()

# Flag to indicate whether we've done the one-shot build
const graph_built = Ref(false)

# —————————————————————————————————————————————————————————————————————————————
#   HELPER FUNCTIONS
# —————————————————————————————————————————————————————————————————————————————

# Utility to append to adjacency list
function _append!(dict::IdDict{Type,Vector{Type}}, key::Type, value::Type)
    if haskey(dict, key)
        push!(dict[key], value)
    else
        dict[key] = [value]
    end
end

# Check if two transformation kinds can be fused
function _can_fuse(kind1::Symbol, kind2::Symbol)
    return kind1 in (:static_matrix, :static_affine) && 
           kind2 in (:static_matrix, :static_affine)
end

# Fuse two transformations
function _fuse(kind1::Symbol, R1, Δ1, kind2::Symbol, payload2)
    if kind2 == :static_matrix
        R2 = payload2
        # Affine: new_vec = R2 * (R1 * vec + Δ1) = (R2*R1)*vec + R2*Δ1
        return :static_affine, R2 * R1, R2 * Δ1
    else  # :static_affine
        R2, Δ2 = payload2
        # new_vec = R2 * (R1 * vec + Δ1) + Δ2 = (R2*R1)*vec + (R2*Δ1 + Δ2)
        return :static_affine, R2 * R1, R2 * Δ1 .+ Δ2
    end
end

# Initialize a new segment from a single transformation
function _initialize_segment(kind::Symbol, payload)
    if kind == :static_matrix
        return :static_matrix, payload, SVector(0.0, 0.0, 0.0)
    elseif kind == :static_affine
        R, Δ = payload
        return :static_affine, R, Δ
    else
        # Dynamic/callable transformations can't be part of a fused segment
        return kind, nothing, nothing
    end
end

# —————————————————————————————————————————————————————————————————————————————
#   ATOMIC-EDGE MACRO
# —————————————————————————————————————————————————————————————————————————————

"""
    @transformation KIND [priority=1.0] Source => Target begin
        …kernel code…
    end

KIND can be:
  :static_matrix   – returns a 3×3 SMatrix
  :static_affine   – returns (R::SMatrix, Δ::SVector)
  :dynamic_matrix  – returns (src_frame,dst_frame)->SMatrix
  :dynamic_affine  – returns (src_frame,dst_frame)->(SMatrix,SVector)
  :callable        – returns (vec,src_frame,dst_frame)->new_vec
"""
macro transformation(args...)
    # Parse arguments - could be (kind, expr, body) or (kind, priority, expr, body)
    if length(args) == 3
        kind_sym, ex, body = args
        priority = 1.0
    elseif length(args) == 4 && args[2].head == :(=) && args[2].args[1] == :priority
        kind_sym = args[1]
        priority = args[2].args[2]
        ex = args[3]
        body = args[4]
    else
        error("Invalid @transformation syntax")
    end
    
    # Expect syntax `A => B`
    @assert ex.head === :call && ex.args[1] === :(=>)
    
    Source, Target = ex.args[2:3]
    
    quote
        let src = $(esc(Source)), tgt = $(esc(Target))
            # 1. Record the kind, payload, and priority
            edge_kind[(src, tgt)] = $(esc(kind_sym))
            edge_payload[(src, tgt)] = $(esc(body))
            edge_priority[(src, tgt)] = $(esc(priority))
            
            # 2. Update adjacency
            _append!(adjacency_list, src, tgt)
            
            # 3. If we've already built once, patch incrementally
            if graph_built[]
                update_paths!(src, tgt, $(esc(kind_sym)))
            end
        end
        nothing
    end
end

# —————————————————————————————————————————————————————————————————————————————
#   PUBLIC GENERIC TRANSFORM
# —————————————————————————————————————————————————————————————————————————————

# Core engine signature: transforms a 3-vector from frame From to frame To
function transform(vec::SVector{3}, ::Type{From}, ::Type{To}) where {
    From<:AbstractFrame,
    To<:AbstractFrame}
    throw(MethodError(transform, (typeof(vec), From, To)))
end

# Convenience method that accepts Cartesian and returns Cartesian
function transform(cart::Cartesian, from::Type{F1}, to::Type{F2}) where {
    F1<:AbstractFrame, F2<:AbstractFrame}
    vec = SVector(cart.x, cart.y, cart.z)
    result = transform(vec, from, to)
    return Cartesian(result[1], result[2], result[3])
end

# —————————————————————————————————————————————————————————————————————————————
#   ONE-SHOT BUILD: build_paths!
# —————————————————————————————————————————————————————————————————————————————

"""
    build_paths!()

Traverse every registered atomic edge and emit fully-fused methods for all 
reachable (Source,Target) pairs. Should be called once (e.g. in `__init__`).
"""
function build_paths!()
    for (src, dst) in keys(edge_kind)
        # Get the actual edge transformation
        kind = edge_kind[(src, dst)]
        payload = edge_payload[(src, dst)]
        
        # Initialize with the actual transformation, not identity
        init_kind, init_R, init_Δ = _initialize_segment(kind, payload)
        _dfs_segment!(src, dst, init_kind, init_R, init_Δ)
    end
    graph_built[] = true
end

# —————————————————————————————————————————————————————————————————————————————
#   SEGMENT-AWARE DFS: _dfs_segment!
# —————————————————————————————————————————————————————————————————————————————

"""
    _dfs_segment!(start, node, seg_kind, R_acc, Δ_acc)

Recursively walk the graph:

1. Fuse current `(seg_kind,R_acc,Δ_acc)` with the edge `(node→next)`.
2. If the edge kind is static_matrix or static_affine:
     – continue fusing (update R_acc, Δ_acc).
   Else (dynamic or callable):
     – emit the *completed* segment as code;
     – start a fresh segment on the barrier node.
3. At each `node == destination`, emit any remaining fused segment.
"""
function _dfs_segment!(start::Type, node::Type, 
                       seg_kind::Symbol,
                       R_acc, Δ_acc)
    # For every outgoing edge from `node`
    for next in get(adjacency_list, node, Type[])
        key = (node, next)
        kind2 = edge_kind[key]
        payload2 = edge_payload[key]
        
        # Determine fusion or barrier
        if _can_fuse(seg_kind, kind2)
            # Fuse into same segment
            new_kind, new_R, new_Δ = _fuse(seg_kind, R_acc, Δ_acc, kind2, payload2)
            _dfs_segment!(start, next, new_kind, new_R, new_Δ)
        else
            # Barrier: emit current segment for start->node, then start new segment
            _maybe_emit(start, node, seg_kind, R_acc, Δ_acc)
            # Next segment seeds with this single hop
            seed_kind, seed_R, seed_Δ = _initialize_segment(kind2, payload2)
            _dfs_segment!(start, next, seed_kind, seed_R, seed_Δ)
        end
    end
    
    # If no further edges, finalize the current segment if we reached somewhere new
    if node != start || seg_kind != :static_matrix || R_acc != I || any(!=(0), Δ_acc)
        _maybe_emit(start, node, seg_kind, R_acc, Δ_acc)
    end
end

# —————————————————————————————————————————————————————————————————————————————
#   EMIT OR SKIP A COMPILED METHOD: _maybe_emit
# —————————————————————————————————————————————————————————————————————————————

"""
    _maybe_emit(src, dst, seg_kind, R, Δ)

If we haven't yet emitted a method for (src,dst), generate one 
`transform(vec, src, dst)` that runs the fused segment.

- For :static_matrix → `return R * vec`
- For :static_affine → `return R * vec .+ Δ`
- Otherwise (dynamic/callable barrier) → use stored callable
"""
function _maybe_emit(src::Type, dst::Type,
                     seg_kind::Symbol,
                     R, Δ)
    key = (src, dst)
    if key ∉ emitted_paths
        if seg_kind == :static_matrix
            # Only emit if it's not identity (optimization)
            if R != I
                # Create a let block to capture R as compile-time constant
                @eval transform(v::SVector{3}, ::Type{$src}, ::Type{$dst}) = $R * v
            elseif src != dst
                # Identity transform but different frames
                @eval transform(v::SVector{3}, ::Type{$src}, ::Type{$dst}) = v
            end
        elseif seg_kind == :static_affine
            # Create a let block to capture R and Δ as compile-time constants
            let Rc = R, Δc = Δ
                @eval transform(v::SVector{3}, ::Type{$src}, ::Type{$dst}) = $Rc * v .+ $Δc
            end
        elseif seg_kind in (:dynamic_matrix, :dynamic_affine, :callable)
            # For dynamic transforms, we need the actual edge payload
            edge_key = (src, dst)
            if haskey(edge_payload, edge_key)
                payload = edge_payload[edge_key]
                if seg_kind == :dynamic_matrix
                    @eval function transform(v::SVector{3}, src::$src, dst::$dst)
                        matrix_fn = $payload
                        R = matrix_fn(src, dst)
                        return R * v
                    end
                elseif seg_kind == :dynamic_affine
                    @eval function transform(v::SVector{3}, src::$src, dst::$dst) 
                        affine_fn = $payload
                        R, Δ = affine_fn(src, dst)
                        return R * v .+ Δ
                    end
                else  # :callable
                    @eval function transform(v::SVector{3}, src::$src, dst::$dst)
                        callable = $payload
                        return callable(v, src, dst)
                    end
                end
            end
        end
        push!(emitted_paths, key)
    end
end

# —————————————————————————————————————————————————————————————————————————————
#   INCREMENTAL PATCH: update_paths!
# —————————————————————————————————————————————————————————————————————————————

"""
    update_paths!(new_src, new_dst, new_kind)

When a new atomic edge is added after the one-shot build, re-explore affected 
paths. This is a simplified version that just re-runs DFS from the new edge.
"""
function update_paths!(new_src::Type, new_dst::Type, new_kind::Symbol)
    # Start a new DFS from this edge
    seed_kind, seed_R, seed_Δ = _initialize_segment(new_kind, edge_payload[(new_src, new_dst)])
    _dfs_segment!(new_src, new_dst, seed_kind, seed_R, seed_Δ)
end

# —————————————————————————————————————————————————————————————————————————————
#   EXPORT
# —————————————————————————————————————————————————————————————————————————————


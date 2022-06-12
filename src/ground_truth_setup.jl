
######### constants for simulated data #########
const dat_t = Float32;          # data_matrix_type

# concentration parameter for Dirichlet that generate Categorical Distributions
const motif_alpha_pc = 0.1;             
const DNA_char_cap = ['A','C','G','T'];
const DNA_char_nocap = ['a','c','g','t'];
const mixture_weights_alpha = 50;               # α for generating mixture weight using Dirichlet
                                                # assume for now mixture weights are roughly the same            
                                                # Note that modes are vector of unit ranges
const percent_complement = .5;                  # probability that a motif is represented in its reverse complement direction
const DNA_complement = Dict('A'=>'T','C'=>'G','G'=>'C','T'=>'A');                     

const background = Categorical(0.25*ones(4));   # simple i.i.d flat background for now
const bg_array = [.25, .25, .25, .25];
const sample_near_center = 0;                   #= increment and decrement the avail. position 
                                                    from the start and the end, resp. =#

const v_f64 = Vector{Float64};
const cat_d = Categorical{Float64, v_f64};
const u_conc = Union{Float64, v_f64, Nothing};

# tests
const num_seq = 800;
const num_expr_each = 1;

# dummy reference # might change the type later for more flexibility
const dummy = Dict('A'=>Array{dat_t}([1, 0, 0, 0]), 
                   'C'=>Array{dat_t}([0, 1, 0, 0]),
                   'G'=>Array{dat_t}([0, 0, 1, 0]), 
                   'T'=>Array{dat_t}([0, 0, 0, 1]));


################################################



#=
e.g. max([1:5, 2:9]) returns 9
=#
function max_v(v::Vector{UnitRange{Int}})
    _max = -Inf;
    for vec in v
        η = maximum(vec);
        _max = η > _max ? η : _max;
    end
    return _max
end


# product categorical 
struct product_categorical
    pc::Vector{cat_d}      # product categorical 
    len::Int              # length of the product
    dim::Int              # dimension for each categorical in the product
    alpha::u_conc           # concentration parameter; 
             
    function product_categorical(L::Int, 
                                 dim::Int, 
                                 alpha::Union{Real, v_f64})             
        if typeof(alpha) == v_f64 
            length(alpha) != dim && error("alpha must be of length same as the dimension");                            
        end
        D = typeof(alpha) == v_f64 ? Dirichlet(alpha) : Dirichlet(dim, alpha);
        new([Categorical(rand(D)) for _ in 1:L], L, dim, alpha)
    end

    function product_categorical(mat::Matrix{T}) where T <: Real
        L, width = size(mat);
        @assert width == 4 "the input matrix must have width 4"
        new([Categorical(mat[i,:]) for i = 1:L], L, 4, nothing)
    end
end


# view(v::product_categorical) = [round(v.pc[j].p[i], digits=2) for j = 1:v.len, i = 1:v.dim];
Base.rand(v::product_categorical) = [rand(v.pc[i]) for i = 1:v.len];


# motifs definitions
struct single_block_motif
    P::product_categorical # P stands for probability measure
    len::Int64
    function single_block_motif(L::T) where T <: Integer
        new(product_categorical(L, 4, motif_alpha_pc), L)
    end
    function single_block_motif(mat::Matrix{T}) where T <: Real
        new(product_categorical(mat), size(mat,1));
    end
end

# example ----------------------------------
# a = single_block_motif(5);
# view(a)
# ------------------------------------------

# an intermediate type for gap motifs
struct k_block_motif
    P::Vector{single_block_motif}
    K::Int64
    """
    K: number of blocks in the motif
    L_vec: vector of lengths L's for each product categorical
    """
    function k_block_motif(length_vec::Vector{Int})
        K = length(length_vec);
        new([single_block_motif(length_vec[i]) for i = 1:K], K)
    end
end

#= 
motif with k blocks and gaps in between every adjacent blocks
Input:
    L_vec:      Length of each of the blocks of the motif
    gap_vec:    Length of each of the gap in between blocks --
                gap_vecᵢ specifies the gap in between ith and 
                (i+1)th part        
    gap_dist:   Multinomial that generates the alphabets in an i.i.d fashion
=# 

struct gapped_k_block_motif
    P_motifs::k_block_motif
    gap_len::Vector{Int}    # the maximal gap length that can occur in between each blocks
    P_gap::Categorical{Float64, Vector{Float64}}
    K::Int
    function gapped_k_block_motif(L_vec::Vector{T}, gap_vec::Vector{T}, gap_dist=bg_array) where T <: Integer       
        K = length(L_vec);
        @assert length(gap_vec) == K-1 "Length of the gap_vec (2nd argument) must be length(L_vec)-1"
        new(k_block_motif(L_vec), gap_vec, Categorical(gap_dist), K)
    end
end


struct mixture_gapped_k_block_motif
    modes::Vector{UnitRange{Int}} # inclusive
    mixture_weights::cat_d
    motif::gapped_k_block_motif
    num_modes::Int
    function mixture_gapped_k_block_motif(
        L_vec::Vector{Int}, 
        gap_vec::Vector{Int}, 
        modes::Vector{UnitRange{Int64}},
        gap_dist=bg_array
    )   
        num_modes = length(modes);
        new(modes, 
            Categorical(rand(Dirichlet(num_modes, mixture_weights_alpha))),
            gapped_k_block_motif(L_vec, gap_vec, gap_dist),         
            num_modes
        )
    end
end

# a simpler mixture motif than the one above

struct mixture_k_block_motif
    modes::Vector{UnitRange{Int}} # inclusive
    mixture_weights::cat_d
    motif::single_block_motif
    num_modes::Int
    function mixture_k_block_motif(modes::Vector{UnitRange{Int64}})   
        K = max_v(modes);
        num_modes = length(modes);
        new(modes, 
            Categorical(rand(Dirichlet(num_modes, mixture_weights_alpha))),
            single_block_motif(K),            
            num_modes
        )
    end
end

# motif type
const motif_type = Union{single_block_motif, 
                         gapped_k_block_motif, 
                         mixture_k_block_motif,
                         mixture_gapped_k_block_motif
};

# type for a (simulated) single dna background string with a motif in it 
const sim_dna_str_w_motif = NamedTuple{(:str, :motif_where, :mode, :complement), 
                    Tuple{String, UnitRange{Int64}, Int64, Bool}
                    }; # mode = 0 => no motif
const public_dna_str_w_motif = NamedTuple{(:str, :motif_where, :mode), 
                    Tuple{String, UnitRange{Int64}, Int64}
                    };

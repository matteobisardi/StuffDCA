export fasta2matrix, map_gaps, read_par_BM, stupid, stupid2

stupid() = 8
stupid2() = 5

"""
  function fasta2matrix(filename::AbstractString; max_gap_fraction = 1)


  Reads a fastafile and returnsa matrix MxN.
  Might be cool to be able to choose the alphabet
  and to avoid having a different function for a file
  already containing numbers. 

"""

function fasta2matrix(filename::AbstractString; max_gap_fraction = 1)
    f = FastaReader(filename)

    max_gap_fraction = Float64(max_gap_fraction)

    # pass 1

    seqs = Int[]
    inds = Int[]
    fseqlen = 0

    for (name, seq) in f
        ngaps = 0
        if f.num_parsed == 1
            ls = length(seq)
            resize!(inds, ls)
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    fseqlen += 1
                    inds[fseqlen] = i
                    c == '-' && (ngaps += 1)
                end
            end
        else
            ls = length(seq)
            ls == length(inds) || error("inputs are not aligned")
            tstfseqlen = 0
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    tstfseqlen += 1
                    inds[tstfseqlen] == i || error("inconsistent inputs")
                    c == '-' && (ngaps += 1)
                end
            end
            tstfseqlen == fseqlen || error("inconsistent inputs")
        end
        ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
    end

    length(seqs) > 0 || error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

    # pass 2
    Z = Array{Int16}(undef, fseqlen, length(seqs))
    seqid = 1
    for (name, seq) in f
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue
        for i = 1:fseqlen
            c = seq[inds[i]]
            Z[i, seqid] = letter2num(c)
        end
        seqid += 1
    end
    @assert seqid == length(seqs) + 1

    close(f)

    return Z'
end


let alphabet = [ 1, 21, 2, 3, 4, 5, 6, 7, 8, 21, 9, 10, 11, 12, 21, 13, 14, 15, 16, 17, 21, 18, 19, 21, 20]
               # A, B,  C, D, E, F, G, H, I, J,  K, L,  M,  N,  O,  P,  Q,  R,  S,  T,  U,  V,  W,  X,  Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
     end
end



"""
        function map_gaps(matrix)
        
        Returns a matrix with 21's instead of 0's,
        to change between different mappings.        

"""

function map_gaps(matrix)
        for i in 1:size(matrix)[1]
		for j in 1:size(matrix)[2]
			if (matrix[i, j] == 0) 
				matrix[i, j] = 21
			end
		end
	end
	return matrix
end



export read_par_BM

"""
	read_par_BM(path::AbstractString, q, N)

	Reads the parameters of a Potts model in format
	
	J i j a b
	...
	h i a 

	and returns them in tensor format.
	J[a, b, i, j] is such that  1 == "A", ... "21" == "-",
	and the same for  h[a, i].

"""

function read_par_BM(path::AbstractString, q, N)
	data = readdlm(path,' ', use_mmap = true)[:, 2:6]
	J = Array{Float64}(undef, q, q, N, N)
	h = Array{Float64}(undef, q, N)
	n_J = Int(q*q*N*(N-1)/2)
	n_h = q*N

	for k in 1:n_J
		i, j, a, b, par_j = data[k, :]
		i += 1
		j += 1
		a == 0 && (a = 21)
		b == 0 && (b = 21)
		J[a, b, i, j] = par_j
	end

	for l in (n_J + 1): n_h + n_J
		i, a, par_h = data[l, :]
		i += 1
		a == 0 && (a = 21)
		h[a, i] = par_h
	end

	return h, J
end


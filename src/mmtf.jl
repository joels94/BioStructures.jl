export
    MMTFDict

"""
A MMTF dictionary.
Keys are field names as a `String` and values are various types.
"""
struct MMTFDict
    dict::Dict{String, Any}
end

MMTFDict() = MMTFDict(Dict())

Base.getindex(mmtf_dict::MMTFDict, field::AbstractString) = mmtf_dict.dict[field]

function Base.setindex!(mmtf_dict::MMTFDict,
                    val,
                    field::AbstractString)
    mmtf_dict.dict[field] = val
    return mmtf_dict
end

Base.keys(mmtf_dict::MMTFDict) = keys(mmtf_dict.dict)
Base.values(mmtf_dict::MMTFDict) = values(mmtf_dict.dict)
Base.haskey(mmtf_dict::MMTFDict, key) = haskey(mmtf_dict.dict, key)

function Base.show(io::IO, mmtf_dict::MMTFDict)
    print(io, "MMTF dictionary with $(length(keys(mmtf_dict))) fields")
end


function Base.read(input::IO,
            ::Type{MMTF};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true)
    d = parsemmtf(input)
    struc = ProteinStructure(structure_name)
    model_i = 0
    chain_i = 0
    group_i = 0
    atom_i = 0
    for modelchaincount in d["chainsPerModel"]
        model_i += 1
        for ci in 1:modelchaincount
            chain_i += 1
            for gi in 1:d["groupsPerChain"][chain_i]
                group_i += 1
                group = d["groupList"][d["groupTypeList"][group_i]]
                for ai in 1:group["atomNameList"]["length"]
                    atom_i += 1
                    unsafe_addatomtomodel!(struc[model_i], AtomRecord(
                        het,
                        id,
                        name,
                        altlocid,
                        resname,
                        chain,
                        resnum,
                        d["insCodeList"][group_i],
                        [
                            d["xCoordList"][atom_i],
                            d["yCoordList"][atom_i],
                            d["zCoordList"][atom_i],
                        ],
                        d["occupancyList"][atom_i],
                        d["bFactorList"][atom_i],
                        group["elementList"][ai],
                        group["formalChargeList"][ai]
                    ))
                end
            end
        end
    end
    fixlists!(struc)
    return struc
end

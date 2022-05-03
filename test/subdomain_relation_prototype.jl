using ModelingToolkit, Symbolics, DomainSets

struct VarDomainSubdomainRelation
    A::Symbolics.VarDomainPairing
    B::Symbolics.VarDomainPairing

    function VarDomainSubdomainRelation(A::Symbolics.VarDomainPairing, B::Symbolics.VarDomainPairing)
        if A ⊆ B
            new(A, B)
        else
            throw(ArgumentError("A.domain must be a subset of B.domain"))
        end
    end
end

function subset(A::Symbolics.VarDomainPairing, B::Symbolics.VarDomainPairing)::VarDomainSubdomainRelation
    VarDomainSubdomainRelation(A, B)
end

function Base.issubset(A::Symbolics.VarDomainPairing, B::Symbolics.VarDomainPairing)::Bool
    Base.issubset(A.domain, B.domain)
end

⊂(A::Symbolics.VarDomainPairing, B::Symbolics.VarDomainPairing) = subset(A, B)

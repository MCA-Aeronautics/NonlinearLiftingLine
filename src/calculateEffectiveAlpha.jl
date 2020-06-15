function calculateEffectiveAlpha(freestream,inducedVelocity,anglesOfAttack)

    inducedAlpha = zeros(length(freestream[:,1]))
    effectiveAlphas = similar(inducedAlpha)
    for i = 1:length(freestream[:,1])
       
        inducedAlpha[i] = atan(inducedVelocity[i],freestream[i,1]) # Just x-component of freestream

        effectiveAlphas[i] = anglesOfAttack[i] .+ inducedAlpha[i]

    end

    return effectiveAlphas

end
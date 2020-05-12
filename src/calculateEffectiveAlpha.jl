function calculateEffectiveAlpha(freestream,inducedVelocity,angle)

    inducedAlpha = zeros(length(freestream))
    for i = 1:length(freestream)
       
        inducedAlpha[i] = atan(inducedVelocity[i],freestream[i])

        #println("Induced Angle at panel ",i," = ",inducedAlpha[i]*180/pi)

    end

    effectiveAlpha = angle .+ inducedAlpha

    return effectiveAlpha

end
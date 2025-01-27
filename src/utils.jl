"""
    normalizecols!(A,p=2)

Normalize the columns of a matrix `A` in place such that each column has p-norm equal to one.
"""
function normalizecols!(A,p=2)
    for col in axes(A,2)
        A[:,col] = normalize(A[:,col],p)
    end
    return A
end

"""
    normalizecols(A,p=2)

Normalize the columns of a matrix `A` such that each column has p-norm equal to one.
"""
function normalizecols(A,p=2)
    B = copy(A)
    for col in axes(A,2)
        B[:,col] = normalize(A[:,col],p)
    end
    return B
end

"""
    scalecols!(A,scalingfactors)

Scale each column of a matrix `A` in place by the corresponding entry in the vector `scalingfactors`.
"""
function scalecols!(A,scalingfactors)
    for col in axes(A,2)
        A[:,col] = A[:,col] ./ scalingfactors[col]
    end
    return A
end

"""
    invgeomsum(y,K)

Returns the value `x` for which ``y = 1 + x + x^2 + ... + x^K``.
"""
function invgeomsum(y,K)
    # we only want to do this for positive arguments
    if y < 1
        return NaN
    end
    # for K+1 we have a simple solution
    if y == K+1
        return 1.
    end
    # set the bracketing interval
    if y < K+1
        xspan = (0., 1.)
    else
        # in this case h(x)=x^{K+1} is a lower bound on geomsum(x,K) and its inverse an upper bound for x
        xspan = (1., y^(1/(K+1)))
    end
    # define the nonlinear root finding problem
    f(u,p) = geomsum(u,p) - y
    prob = IntervalNonlinearProblem(f, xspan, K)
    # solve the problem and return the result
    sol = solve(prob)
    return sol.u
end

# function invgeomsum_pol(y,K)
#     # check if we are in the singular case y=K+1
#     if y == K+1
#         return 1.
#     end
#     # coefficients of the polynomial  p(z) = z^{K+1} - y*z + y - 1, in ascending order, from the lowest to the highest
#     p = [0. for _ in 0:K+1]
#     p[1] = y - 1
#     p[2] = -y
#     p[end] = 1.
#     # find the roots of the polynomial
#     r = roots(p)
#     # 1 is always a root, but has been dealt with above
#     # find the root that is not 1 and returns the correct geometric sum
#     gs = geomsum.(r, K)
    
#     return r
# end

"""
    geomsum(x,K)

Returns the sum of the geometric series ``1 + x + x^2 + ... + x^K``.
"""
function geomsum(x,K)
    if x==1
        return K+1
    else
        return (1-x^(K+1))/(1-x)
    end 
end

"""
    varconfint(s², n, α)

Compute the confidence interval for the variance of a normal distribution with estimated variance `s²` based on a sample of size `n` at significance level `ϵ`.
"""
varconfint(s², n, ϵ) = (n-1)*s²./quantile(Chisq(n-1), [1-ϵ/2, ϵ/2])

"""
    stdconfint(s², n, ϵ)

Compute the confidence interval for the standard deviation of a normal distribution with estimated standard deviation `s` based on a sample of size `n` at significance level `ϵ`.
"""
stdconfint(s, n, ϵ) = sqrt.(varconfint(s^2, n, ϵ))


"""
    cleanroc(gtpval,scores,pcut)

Compute ROC curve and rank correlation for a given set of `scores` and ground truth p-values `gtpval`, considering all points with p-values below `pcut` as true. NaNs are removed from the ground truth p-values and scores before computing the ROC curve. The ROC curve is computed with `n` points.
"""
function cleanroc(gtpval,scores,pcut=0.05,n=100)
    # remove NaNs
    tf = .!isnan.(gtpval) .& .!isnan.(scores)
    gtpval = gtpval[tf]
    scores = scores[tf]
    # compute the rank correlation between the scores and the ground truth p-values, where low p-values should correspond to high scores
    c = corspearman(1. .- gtpval,scores)
    # define the binary ground truth
    gt = gtpval .< pcut
    # compute the ROC curve
    if n>0
        thres = range(0.0,maximum(scores),n)
    else
        thres = sort(unique(scores))
    end
    rc = roc(gt,scores,thres)

    return rc, c
end
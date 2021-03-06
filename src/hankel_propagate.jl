#=

\[U(r,\phi) = \sum_{p=-\infty}^\infty
c_p(-j)^p
\exp(jp\phi)U_p(r)\]
where
\[c_p = \frac{1}{2\pi}
\int\limits_0^{2\pi} d\theta\;
\exp(-jp\theta)g_\Theta(\theta)\]
and
\[U_p(r) = -j\frac{k}{2\pi z}
\exp(jkz)\exp\left(j\frac{k}{2z}r^2\right)
\mathcal{H}_p\{g_R(\rho)\}(kr/2\pi z)\]
=#

using Hankel
using DSP

# This only implements the radial transform for now. If U is an array,
# it is assumed to be evaluated at the correct radial points.
function hankel_propagate(U::Union{Function,AbstractArray},
                          ρmax::Real, N,
                          k = 1.0, z = 1.0,
                          p = 0)
    guizar = Guizar(p, ρmax, N)
    r, H = guizar(U)
    μ = k/(2π*z)
    UU = -im*μ*exp(im*k*z)*exp(im*k*r.^2/2z).*H

    r /= μ
    r, UU
end

function naive_hankel(U::AbstractVector,
                      ρ::AbstractVector,
                      p)
    Δρ = ρ[2]-ρ[1]
    r = linspace(0,1,length(ρ))*1.0/Δρ
    r, 2π * map(eachindex(r)) do i
        sum(ρ.*besselj(p, 2π*r[i]*ρ).*U)*Δρ
    end
end

function naive_hankel_propagate(U::AbstractVector,
                                ρ::AbstractVector,
                                k = 1.0, z = 1.0,
                                p = 0)
    assert(all(ρ .>= 0))
    r, H = naive_hankel(U, ρ, p)
    μ = k/(2π*z)
    UU = -im*μ*exp(im*k*z)*exp(im*k*r.^2/2z).*H

    r /= μ
    r, UU
end

export hankel_propagate

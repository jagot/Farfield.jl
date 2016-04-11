#= We follow the notation of Goodman, Introduction to Fourier optics, 1996:


\[ U(x,y) = -j \frac{k}{2\pi z}
\exp(jkz)\exp\left[j\frac{k}{2z}(x^2+y^2)\right]
\mathcal{F}\{U(\xi,\eta)\}(kx/2\pi z,ky/2\pi z) \]
=#

using DSP

function fft_propagate(U::AbstractVector,
                       ξ::AbstractVector,
                       k = 1.0, z = 1.0)
    Δξ = ξ[2]-ξ[1]
    x = fftfreq(length(ξ), 1.0/Δξ)
    μ = k/(2π*z)
    UU = -im*μ*exp(im*k*z)*exp(im*k*x.^2/2z).*fft(U) * Δξ
    x /= μ
    x,UU
end

function fft_propagate(U::AbstractMatrix,
                       ξ::AbstractVector,
                       η::AbstractVector,
                       k = 1.0, z = 1.0)
    Δξ = ξ[2]-ξ[1]
    Δη = η[2]-η[1]
    x = fftfreq(length(ξ), 1.0/Δξ)
    y = fftfreq(length(η), 1.0/Δη)
    μ = k/(2π*z)
    UU = -im*μ*exp(im*k*z)*exp(im*k*(x.^2 + y.^2)/2z).*fft(U) * Δξ * Δη
    x /= μ
    y /= μ
    x,y,UU
end

export fft_propagate

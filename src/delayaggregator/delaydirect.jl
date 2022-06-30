mutable struct DelayDirectJumpAggregation{T,S,F1,F2,RNG,IType} <: AbstractDSSAJumpAggregator
    next_jump::Int
    prev_jump::Int
    next_jump_time::T
    end_time::T
    cur_rates::Vector{T} # cur_rates here is the cumsum of rates
    sum_rate::T
    ma_jumps::S
    rates::F1
    affects!::F2
    save_positions::Tuple{Bool,Bool}
    rng::RNG
    next_delay::Vector{Int}
    num_next_delay::Vector{Int}
    time_to_next_jump::T
    next_delay_time::T
    shadow_integrator::IType
    copied::Bool
end
mutable struct ShadowIntegrator{uType,chanType,T}
    u::uType
    de_chan::chanType
    delayjumpsets::DelayJumpSet
    cur_rates::Vector{T}
end

function DelayDirectJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; u0, kwargs...) where {T,S,F1,F2,RNG}
    ttnj = zero(et)
    ndt = zero(et)
    shadow_integrator = ShadowIntegrator{typeof(u0),Vector{Vector{T}}, T}(copy(u0), [Vector{T}()], DelayJumpSet(Dict(),Dict(),Dict()), copy(crs))
    nd = Int64[]
    nnd = Int64[] # in the Direct Method the number of next delay equals always 1
    DelayDirectJumpAggregation{T,S,F1,F2,RNG,typeof(shadow_integrator)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, nd, nnd, ttnj, ndt, shadow_integrator, false)
end

function aggregate(aggregator::DelayDirect, u, p, t, end_time, constant_jumps, ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using function wrappers
    rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)
    build_jump_aggregation(DelayDirectJumpAggregation, u, p, t, end_time, ma_jumps, rates, affects!, save_positions, rng; u0=copy(u), kwargs...)
end

function initialize!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    p.shadow_integrator.delayjumpsets = integrator.delayjumpsets
    generate_jumps!(p, integrator, u, params, t)
    nothing
end

@inline function execute_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    update_state_delay_direct!(p, integrator, integrator.u, t)
    nothing
end

function generate_jumps!(p::DelayDirectJumpAggregation, integrator, u, params, t)
    generate_time_to_next_jump!(p, integrator, params, t)
    @fastmath p.next_jump_time = t + p.time_to_next_jump
    @inbounds p.next_jump = searchsortedfirst(p.cur_rates, rand(p.rng) * p.sum_rate) # 
    nothing
end
"""
    Create time_to_next_jump based on the shawdow integrator, the goal is to use p.copied to avoid the case when no changes happened wrt integrator.u such that the allocation can be optimised
"""
@inline function generate_time_to_next_jump!(p::DelayDirectJumpAggregation, integrator, params, t)

    fill_cum_rates_and_sum!(p, integrator.u, params, t)
    r1 = rand(p.rng)
    if isempty(reduce(vcat, integrator.de_chan))
        ttnj = -log(r1) / p.cur_rates[end]
        ttnj_last = ttnj
        p.copied = false
    else
        find_next_delay_num!(p, integrator.de_chan)
        aₜ = p.cur_rates[end] * p.next_delay_time
        F = one(t) - exp(-aₜ)
        aₜ_ = zero(aₜ)
        prev_T1 = zero(t)
        cur_T1 = zero(t)
        cur_T1 += p.next_delay_time
        if F < r1 # if enter the loop
            shadow_integrator = p.shadow_integrator
            shadow_integrator.u = copy(integrator.u) 
            shadow_integrator.de_chan = deepcopy(integrator.de_chan)
            p.copied = true
        else
            p.copied = false
        end
        while F < r1
            shift_delay_channel!(shadow_integrator.de_chan, p.next_delay_time)
            update_delay_channel!(shadow_integrator.de_chan)
            update_delay_complete!(p, shadow_integrator)

            calculate_sum_rate!(p, shadow_integrator, shadow_integrator.u, params, t + cur_T1)

            find_next_delay_num!(p, shadow_integrator.de_chan)
            prev_T1 = cur_T1 # to avoid cur_T1 = Inf 
            cur_T1 += p.next_delay_time

            aₜ_ = aₜ # backup aₜ
            aₜ += shadow_integrator.cur_rates[end] * (p.next_delay_time)
            F = one(t) - exp(-aₜ)
        end
        sum_ = p.copied ? shadow_integrator.cur_rates[end] : p.sum_rate
        ttnj_last = (-log(one(t) - r1) - aₜ_) / sum_
        ttnj = prev_T1 + ttnj_last
    end
    
    if p.copied
        shift_delay_channel!(shadow_integrator.de_chan, ttnj_last)
        update_delay_channel!(shadow_integrator.de_chan)
        fill_cum_rates_and_sum!(p, shadow_integrator.u, params, t + ttnj)
    end
    p.time_to_next_jump = ttnj
    nothing
end



"""
    update_delay_at_tstop_test!

this function will change 
- p.next_delay, p.num_next_delay, p.next_delay_time 
- integrator.u, integrator.de_chan
"""
function update_delay_at_tstop_test!(p, integrator, params, t, tgap)
    # direct_algo!(p, integrator, params, t; tgap = tgap)
    cur_T1 = zero(t)
    prev_T1 = zero(t)
    find_next_delay_num!(p, integrator.de_chan)
    cur_T1 += p.next_delay_time
    while cur_T1 <= tgap
        shift_delay_channel!(integrator.de_chan, p.next_delay_time)
        update_delay_channel!(integrator.de_chan)
        update_delay_complete!(p, integrator)
        find_next_delay_num!(p, integrator.de_chan) 
        prev_T1 = cur_T1 # to avoid cur_T1 = Inf 
        cur_T1 += p.next_delay_time    
    end
    ttnj_last = tgap - prev_T1
    shift_delay_channel!(integrator.de_chan, ttnj_last)
    update_delay_channel!(integrator.de_chan)
    nothing
end

"""
    function fill_cum_rates_and_sum!(p::DelayDirectJumpAggregation, u, params, t)

This function modifies `p.cur_rates` and `p.sum_rate`
"""
function fill_cum_rates_and_sum!(p::DelayDirectJumpAggregation, u, params, t)
    prev_rate = zero(t)
    new_rate = zero(t)
    cur_rates = p.cur_rates
    # mass action rates
    majumps = p.ma_jumps
    idx = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate = cur_rates[i]
    end
    # constant jump rates
    idx += 1
    rates = p.rates
    @inbounds for i in eachindex(p.rates)
        new_rate = rates[i](u, params, t)
        cur_rates[idx] = new_rate + prev_rate
        prev_rate = cur_rates[idx]
        idx += 1
    end

    @inbounds sum_rate = cur_rates[end]
    p.sum_rate = sum_rate
    nothing
end
"""
    Only changes s.sum_rate
"""
function calculate_sum_rate!(p, s::ShadowIntegrator, u, params, t)
    prev_rate = zero(t)
    new_rate = zero(t)
    cur_rates = s.cur_rates
    majumps = p.ma_jumps
    idx = get_num_majumps(majumps)
    @inbounds for i in 1:idx
        new_rate = evalrxrate(u, i, majumps)
        cur_rates[i] = new_rate + prev_rate
        prev_rate = cur_rates[i]
    end
    # constant jump rates
    idx += 1
    rates = p.rates
    @inbounds for i in eachindex(rates)
        new_rate = rates[i](u, params, t)
        cur_rates[idx] = new_rate + prev_rate
        prev_rate = cur_rates[idx]
        idx += 1
    end
end

@inline function update_state_delay_direct!(p::DelayDirectJumpAggregation, integrator, u, t)
    @unpack ma_jumps, next_jump, time_to_next_jump = p
    @unpack delay_trigger_set, delay_interrupt_set = integrator.delayjumpsets
    if p.copied
        integrator.u = copy(p.shadow_integrator.u)
        integrator.de_chan = deepcopy(p.shadow_integrator.de_chan)
    else
        # update_delay_at_tstop_test!(p, integrator, p, t, time_to_next_jump)
        shift_delay_channel!(integrator.de_chan, time_to_next_jump)
        update_delay_channel!(integrator.de_chan)
    end

    num_ma_rates = get_num_majumps(ma_jumps)
    if next_jump <= num_ma_rates # if the next jump is a mass action jump
        if u isa SVector
            integrator.u =  executerx(integrator.u, next_jump, ma_jumps)
        else
            @inbounds executerx!(integrator.u, next_jump, ma_jumps)
        end
    else
        idx = next_jump - num_ma_rates
        @inbounds p.affects![idx](integrator)
    end

    if next_jump in delay_interrupt_set
        update_delay_interrupt!(p, integrator)
    elseif next_jump in delay_trigger_set
        update_delay_trigger!(p, integrator)
    end
    # save jump that was just executed
    p.prev_jump = next_jump
    nothing
end


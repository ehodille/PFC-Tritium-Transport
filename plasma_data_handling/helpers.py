"""
Helper functions for plasma data handling and pulse processing.

These functions are used to create time-dependent boundary conditions
and process pulse profiles for FESTIM simulations.
"""

from scenario import Pulse
from scenario_cplx import Pulse as Pulse_c

def periodic_pulse_function(current_time: float, pulse: Pulse, value, value_off=343.0):
    """Creates bake function with ramp up rate and ramp down rate.

    Args:
        current_time (float): time within the pulse 
        pulse (Pulse): pulse of HISP Pulse class
        value (float): steady-state value 
        value_off (float): value at t=0 and t=final time. 
    """
    
    if current_time == pulse.total_duration:
        return value_off
    elif current_time % pulse.total_duration < pulse.ramp_up:  # ramp up 
        return (value - value_off) / (pulse.ramp_up) * current_time + value_off  # y = mx + b, slope is temp/ramp up time
    elif current_time % pulse.total_duration < pulse.ramp_up + pulse.steady_state:  # steady state
        return value
    else:  # ramp down, waiting
        lower_value = value - (value - value_off)/pulse.ramp_down * (current_time - (pulse.ramp_up + pulse.steady_state))  # y = mx + b, slope is temp/ramp down time
        if lower_value >= value_off: 
            return lower_value
        else: 
            return value_off

def periodic_pulse_cplx_function(current_time: float, pulse: Pulse_c, value, value_off=0):
    '''
    Creates a function with multiple steps with ramp up/down between them

    Args:
        current_time (float): time within the pulse
        pulse (Pulse_cplx): pulse of cplx class
        value (float): steady_state value
        value_off (float): value at t_final and t=0
    '''
    time_prev = 0
    value_prev = value_off
    for r in range(len(pulse.steady_STATE)):
        if time_prev <= current_time < time_prev + pulse.transition[r]:
            return (value * pulse.fraction_ss[r] - value_prev) / pulse.transition[r] * (current_time - time_prev) + value_prev
        elif time_prev + pulse.transition[r] <= current_time < time_prev + pulse.transition[r] + pulse.steady_STATE[r]:
            return value * pulse.fraction_ss[r]
        value_prev = value * pulse.fraction_ss[r]
        time_prev += pulse.transition[r] + pulse.steady_STATE[r]
    if time_prev <= current_time <= time_prev + pulse.transition[-1]:
        return (value_off - value_prev) / pulse.transition[-1] * (current_time - time_prev) + value_prev
    if current_time >= pulse.duration_no_waiting:
        return value_off
    if current_time == pulse.total_duration:
        return value_off




def periodic_step_function(x, period_on, period_total, value, value_off=0.0):
    """
    Creates a periodic step function with two periods.
    
    Args:
        x: time or position value
        period_on: duration when the function is "on"
        period_total: total period duration
        value: value when "on"
        value_off: value when "off"
    
    Returns:
        value or value_off depending on position in period
    """
    if x % period_total < period_on:
        return value
    else:
        return value_off

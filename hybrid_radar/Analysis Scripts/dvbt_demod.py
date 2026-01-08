
import pdb
import scipy 
import numpy as np
import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# import 

def find_symbol_starts(input_sig, S, D, dt):
    """
    Looks for start points of the OFDM symbols
    
    Parameters:
    -----------
    input_sig : numpy array
        Input signal
    S : int
        Symbol parameter
    D : int
        Duration parameter
    dt : float
        Time step
    
    Returns:
    --------
    symbols_starts : numpy array
        Array of symbol start indices
    symbol_phases : numpy array
        Array of symbol phases at start points
    """
    
    # Initialize P1 array
    P1 = np.zeros(len(input_sig), dtype=complex)
    
    # Calculate correlation metric
    for i in range(len(input_sig) - S - D - 1):
        P1[i] = np.sum(np.conj(input_sig[i:i+D+1]) * input_sig[i+S:i+S+D+1])
    
    # Uncomment to plot (optional)
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(np.abs(P1))
    # plt.title('Magnitude of P1')
    # plt.show()
    # 
    # plt.figure()
    # plt.plot(np.angle(P1))
    # plt.title('Phase of P1')
    # plt.show()
    
    P1_angle = np.angle(P1)
    P1 = np.abs(P1)
    
    # Commented out section from MATLAB (kept as reference)
    # t = np.linspace(0, dt * len(P1), len(P1))
    # corr = np.exp(1j * 0.150455 / S / dt * t)
    # 
    # temp = input_sig * corr
    # 
    # P1_temp = np.zeros(len(input_sig), dtype=complex)
    # for i in range(len(input_sig) - S - D - 1):
    #     P1_temp[i] = np.sum(np.conj(temp[i:i+D+1]) * temp[i+S:i+S+D+1])
    # 
    # plt.figure()
    # plt.plot(np.angle(P1_temp))
    # plt.show()
    
    symbols_starts = []
    
    # Find first symbol start
    m = np.max(P1[0:2*S])
    I_new = np.argmax(P1[0:2*S])
    symbols_starts.append(I_new)
    
    # Find subsequent symbol starts
    while True:
        if I_new + D + S + 200 > len(P1):
            break
        
        I_prev = I_new
        
        # Define search window
        search_start = max(0, I_prev + D + S - 1000)
        search_end = min(len(P1), I_prev + D + S + 1001)
        
        # Find maximum in search window
        search_window = P1[search_start:search_end]
        m_new = np.max(search_window)
        I_new = np.argmax(search_window)
        I_new = I_new + search_start
        
        # Check if peak is valid
        if abs(m_new - m) / m < 0.5:
            symbols_starts.append(I_new)
        else:
            print('something has gone wrong with tracking beginning of frame')
    
    # Convert to numpy array
    symbols_starts = np.array(symbols_starts)
    
    # Extract symbol phases
    symbol_phases = np.zeros(len(symbols_starts))
    for i in range(len(symbols_starts)):
        symbol_phases[i] = P1_angle[symbols_starts[i]]
    
    return symbols_starts, symbol_phases

def find_closest_neighbors(array, n):
    """
    Find closest values on either side of n
    
    Parameters:
    -----------
    array : numpy array or list
        Sorted array of numbers (ascending order)
    n : float or int
        Target number
    
    Returns:
    --------
    lower : int
        Index of closest value less than or equal to n
    upper : int
        Index of closest value greater than or equal to n
    """
    
    
    # Convert to numpy array if needed
    array = np.array(array)
    
    # Validate inputs
    if not np.all(array[:-1] <= array[1:]):
        raise ValueError('Array must be sorted in ascending order')
    
    if n < array[0] or n > array[-1]:
        raise ValueError('n must be between the smallest and largest values in the array')
    
    # Find the index where n would be inserted to maintain sorted order
    idx = np.where(array >= n)[0]
    
    # Handle edge cases
    if len(idx) == 0:
        # n is larger than all elements (shouldn't happen due to validation)
        lower = len(array) - 1
        upper = len(array) - 1
    elif array[idx[0]] == n:
        # Exact match found
        lower = idx[0]
        upper = idx[0]
    elif idx[0] == 0:
        # n is smaller than first element (shouldn't happen due to validation)
        lower = 0
        upper = 0
    else:
        # Normal case: n is between two elements
        lower = idx[0] - 1
        upper = idx[0]
    
    return lower, upper

def channel_compensation(XF,CS, carrier_locations, numb_carriers):
    temp = np.abs(CS)
    dont_use_carrier = []
    
    # Find carriers with anomalous magnitudes
    for i in range(1, len(CS) - 1):
        fl_avg = 0.5 * (temp[i-1] + temp[i+1])
        if np.abs(temp[i] - fl_avg) / fl_avg > 0.1:
            dont_use_carrier.append(i)
    
    # Create good carriers list by removing bad ones
    good_carriers = carrier_locations.copy()
    CS_good = CS.copy()
    
    # Remove bad carriers (in reverse order to avoid index shifting)
    for i in sorted(dont_use_carrier, reverse=True):
        if i < len(good_carriers):
            good_carriers = np.delete(good_carriers, i)
            CS_good = np.delete(CS_good, i)
    
    XF_comp = XF.copy()
    Comps = np.zeros(len(XF_comp), dtype=complex)
    
    for i in range(numb_carriers):
        if np.any(carrier_locations == i):
            # Find closest carrier location
            II = np.argmin(np.abs(i - carrier_locations))
            XF_comp[i] = XF[i] / CS[II]
            Comps[i] = CS[II]
        else:
            # Find closest neighbors and interpolate
            lower, upper = find_closest_neighbors(good_carriers, i)
            
            # Linear interpolation for complex values
            CS_interp = interp1([good_carriers[lower], good_carriers[upper]], 
                               [CS_good[lower], CS_good[upper]], i)
            
            XF_comp[i] = XF[i] / CS_interp
            Comps[i] = CS_interp
    
    return XF_comp, Comps

def generate_prbs(n, init_state=None):
    """
    Generate PRBS sequence based on polynomial X^11 + X^2 + 1
    
    Parameters:
    -----------
    n : int
        Number of bits to generate
    init_state : numpy array or list, optional
        Initial 11-bit state as vector (default: all ones)
    
    Returns:
    --------
    prbs_seq : numpy array
        PRBS sequence of length n (column vector)
    
    Examples:
    ---------
    seq = generate_prbs(100)  # Generate 100 bits
    seq = generate_prbs(50, [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1])  # Custom init state
    """
    
    # Default initial state (all ones) if not provided
    if init_state is None:
        init_state = np.ones(11, dtype=int)
    else:
        init_state = np.array(init_state, dtype=int)
    
    # Ensure init_state is the correct length
    if len(init_state) != 11:
        raise ValueError('Initial state must be 11 bits long')
    
    # Initialize shift register (ensure it's a column vector)
    shift_reg = init_state.flatten()
    
    # Preallocate output
    prbs_seq = np.zeros(n, dtype=int)
    
    # Generate PRBS sequence
    # Polynomial: X^11 + X^2 + 1
    # Feedback taps at positions 11 and 2
    for i in range(n):
        # Output is the last bit in the shift register
        prbs_seq[i] = shift_reg[-1]
        
        # Calculate feedback: XOR of tap positions 11 and 2
        # Note: Using indices 10 and 8 for 0-based indexing (positions 11 and 9 in 1-based)
        feedback = shift_reg[10] ^ shift_reg[8]
        
        # Shift register right and insert feedback at the beginning
        shift_reg = np.concatenate(([feedback], shift_reg[:-1]))
    
    return prbs_seq

def get_actual_symbols(XF, num_carriers, sypo, carrier_locations, tps_carriers):
    """
    Get actual symbols from received signal
    
    Parameters:
    -----------
    XF : numpy array
        Frequency domain received signal
    num_carriers : int
        Number of carriers
    sypo : numpy array or list
        Symbol positions (constellation points)
    carrier_locations : numpy array or list
        Locations of pilot carriers
    tps_carriers : numpy array or list
        TPS (Transmission Parameter Signaling) carrier locations
    
    Returns:
    --------
    i_symbol : numpy array
        In-phase symbol indices
    q_symbol : numpy array
        Quadrature symbol indices
    """
    
    # Convert inputs to numpy arrays if needed
    sypo = np.array(sypo)
    carrier_locations = np.array(carrier_locations)
    tps_carriers = np.array(tps_carriers)
    
    # Initialize output arrays
    i_symbol = np.zeros(num_carriers, dtype=int)
    q_symbol = np.zeros(num_carriers, dtype=int)
    
    # Process each carrier
    for i in range(num_carriers):
        # Skip if this is a pilot or TPS carrier
        if np.any(carrier_locations == i):  # or np.any(tps_carriers == i)
            continue
        else:
            # Find closest symbol position for real part
            i_symbol[i] = np.argmin(np.abs(np.real(XF[i]) - sypo))
            
            # Find closest symbol position for imaginary part
            q_symbol[i] = np.argmin(np.abs(np.imag(XF[i]) - sypo))
    
    return i_symbol, q_symbol

def get_channel_state(X_pilots_int,prbs_seq,carrier_locations):

    pilot_locations=(4/3)*2*(0.5-prbs_seq)+1j*0

    CS=X_pilots_int./pilot_locations(carrier_locations)'

def get_all_carrier_locations(offset, existing_carrier_locations, numb_carriers):

    """
    Get all carrier locations based on offset and number of carriers
    
    Parameters:
    -----------
    offset : int
        Offset value (1, 2, 3, or 4)
    existing_carrier_locations : list or numpy array
        Existing carrier locations to merge with new ones
    numb_carriers : int
        Number of carriers
    
    Returns:
    --------
    all_carriers : numpy array
        Unique sorted array of all carrier locations
    """
        
    # Determine start value based on offset
    if offset == 1:
        start = 10
    elif offset == 2:
        start = 13
    elif offset == 3:
        start = 16
    elif offset == 4:
        start = 19
    else:
        start = 10
        print('this should not happen')
    
    # Generate carrier locations with step of 12
    temp = np.arange(start, numb_carriers + 1, 12)
    
    # Combine with existing locations and get unique values
    all_carriers = np.unique(np.concatenate([temp, existing_carrier_locations]))
    
    return all_carriers


def get_carrier_offset(XF_int, guess_offset, offset_probability, offset_order, offset_spacing):
    """
    Determine carrier offset based on XF values and probability metrics
    
    Parameters:
    -----------
    XF_int : numpy array
        Input array (complex or real)
    guess_offset : int
        Initial guess for offset (1-4)
    offset_probability : float
        Probability threshold for offset determination
    offset_order : list or numpy array
        Order of offsets to check (length 4)
    offset_spacing : int
        Spacing between offset positions
    
    Returns:
    --------
    output_offset : int
        Determined carrier offset (1-4)
    """
    
    # Get absolute values
    XF = np.abs(XF_int)
    
    # Initialize temp array
    temp = np.zeros(4)
    
    # Calculate sums for each offset (converting to 0-based indexing)
    temp[0] = (XF[9] + XF[offset_order[0] - 1 + offset_spacing] + 
               XF[offset_order[0] - 1 + 2*offset_spacing] + 
               XF[offset_order[0] - 1 + 3*offset_spacing])
    
    temp[1] = (XF[12] + XF[offset_order[1] - 1 + offset_spacing] + 
               XF[offset_order[1] - 1 + 2*offset_spacing] + 
               XF[offset_order[1] - 1 + 3*offset_spacing])
    
    temp[2] = (XF[15] + XF[offset_order[2] - 1 + offset_spacing] + 
               XF[offset_order[2] - 1 + 2*offset_spacing] + 
               XF[offset_order[2] - 1 + 3*offset_spacing])
    
    temp[3] = (XF[18] + XF[offset_order[3] - 1 + offset_spacing] + 
               XF[offset_order[3] - 1 + 2*offset_spacing] + 
               XF[offset_order[3] - 1 + 3*offset_spacing])
    
    # Find maximum value and its index
    I = np.argmax(temp)
    B = temp[I]
    
    # Get second maximum
    temp2 = temp.copy()
    temp2[I] = -np.inf  # Set max to very small value to find second max
    second = np.max(temp2)
    
    # Calculate guess metric
    myguess = (abs(B - second) / second) * 100
    
    # Determine output offset (convert back to 1-based indexing for output)
    if guess_offset == (I + 1):
        output_offset = I + 1
    elif offset_probability > myguess:
        output_offset = guess_offset
        print("metrics disagree on what is carrier offset")
    else:
        output_offset = I + 1
        print("metrics disagree on what is carrier offset")
    
    return output_offset

def map_symbols_to_new_locations(i_symbol, q_symbol, numb_carriers, prbs_seq, sypo, carrier_locations, tps_carriers):
    """
    Map symbols to new carrier locations
    
    Parameters:
    -----------
    i_symbol : numpy array or list
        In-phase symbol indices
    q_symbol : numpy array or list
        Quadrature symbol indices
    numb_carriers : int
        Number of carriers
    prbs_seq : numpy array or list
        Pseudo-random binary sequence
    sypo : numpy array or list
        Symbol position array
    carrier_locations : numpy array or list
        Locations of carriers
    tps_carriers : numpy array or list
        TPS (Transmission Parameter Signaling) carrier locations
    
    Returns:
    --------
    perfect_symbols : numpy array
        Complex array of mapped symbols
    """
    
    # Initialize complex array
    perfect_symbols = np.zeros(numb_carriers, dtype=complex)
    
    # Convert inputs to numpy arrays for easier checking
    carrier_locations = np.asarray(carrier_locations)
    tps_carriers = np.asarray(tps_carriers)
    
    # Loop through all carriers (convert to 0-based indexing)
    for i in range(numb_carriers):
        # Check if this is a carrier location (MATLAB uses 1-based, Python uses 0-based)
        if np.any(carrier_locations == (i + 1)):
            # Map to PRBS-based symbol
            perfect_symbols[i] = (4/3) * 2 * (0.5 - prbs_seq[i]) + 1j * 0
        
        # Check if this is a TPS carrier
        elif np.any(tps_carriers == (i + 1)):
            if sypo[i_symbol[i]] < 0:
                perfect_symbols[i] = -1
            else:
                perfect_symbols[i] = 1
        
        # Normal symbol mapping
        else:
            perfect_symbols[i] = sypo[i_symbol[i]] + 1j * sypo[q_symbol[i]]
    
    return perfect_symbols

def update_offset(old_offset):
    """
    Update offset value, cycling from 1 to 4
    
    Parameters:
    -----------
    old_offset : int
        Current offset value (1, 2, 3, or 4)
    
    Returns:
    --------
    guess_offset : int
        Updated offset value
    """
    if old_offset == 1 or old_offset == 2 or old_offset == 3:
        guess_offset = old_offset + 1
    else:
        guess_offset = 1
    
    return guess_offset

def produce_IQ_data(perfect_symbols, numb_carriers, sample_length, CS):
    """
    Produce IQ data from perfect symbols using IFFT
    
    Parameters:
    -----------
    perfect_symbols : numpy array
        Complex array of symbols
    numb_carriers : int
        Number of carriers
    sample_length : int
        Length of output sample
    CS : int or float
        (Parameter not used in function body)
    
    Returns:
    --------
    res2_signal_per : numpy array
        Complex IQ signal data
    """
    
    # Apply ifftshift
    temp4_per = np.fft.ifftshift(perfect_symbols)
    
    # Split and zero-pad the signal
    # Convert to 0-based indexing: MATLAB's 1:3409 becomes 0:3409, and 3410:end becomes 3409:
    temp5_per = (sample_length / numb_carriers) * np.concatenate([
        temp4_per[0:3409],
        np.zeros(sample_length - numb_carriers),
        temp4_per[3409:]
    ])
    
    # Apply inverse FFT
    res2_signal_per = np.fft.ifft(temp5_per)
    
    return res2_signal_per



def produce_ideal_ofdm_symbol(input_signal, S, D, numb_carriers, permanant_carriers, tps_carriers, dt):
    """
    Takes in reference (or surveillance signal) and demodulates it and
    recreates the signal based on the actual bits, produces a cleaned
    signal which can be used as reference
    
    Parameters:
    -----------
    input_signal : numpy array
        Complex input signal
    S : int
        Symbol length parameter
    D : int
        Guard interval/delay parameter
    numb_carriers : int
        Number of carriers
    permanant_carriers : numpy array
        Permanent carrier locations
    tps_carriers : numpy array
        TPS carrier locations
    dt : float
        Time step
    
    Returns:
    --------
    output_signal : numpy array
        Complex output signal
    total_offsets : numpy array
        Array of offsets for each symbol
    """
    
    # Initialize output signal
    output_signal = np.zeros(len(input_signal), dtype=complex)
    
    # Initialize parameters
    offset_order = np.array([10, 13, 16, 19])
    offset_spacing = 12
    offset_probability = 0
    offset = 4
    
    # Generate PRBS sequence
    prbs_seq = generate_prbs(numb_carriers)
    
    # Scale factor and symbol positions
    scale_factor = 6.3636
    sypo = np.array([-7, -5, -3, -1, 1, 3, 5, 7]) / scale_factor
    
    # Find symbol starts
    symbol_starts, symbol_phases = find_symbol_starts(input_signal, S, D, dt)
    
    # Save symbol starts data (uncomment if needed)
    # np.savez("symbol_starts_data_1.npz", symbol_starts=symbol_starts, symbol_phases=symbol_phases)
    
    # Initialize total offsets
    total_offsets = np.zeros(len(symbol_starts))
    
    # Process each symbol
    for i in range(len(symbol_starts)):
        # Extract signal segment (convert to 0-based indexing)
        start_idx = symbol_starts[i] + D
        end_idx = symbol_starts[i] + D + S + 1
        msignal = input_signal[start_idx:end_idx]
        check_me = msignal.copy()
        
        # Create time vector and correction
        t = np.linspace(0, dt * len(msignal), len(msignal))
        corr = np.exp(-1j * symbol_phases[i] / S / dt * t)
        msignal = msignal * corr
        
        # Interpolation
        XX = np.linspace(-1 - 5 + 7, 6819 - 5 + 4, len(msignal))
        XV = np.linspace(1, numb_carriers, numb_carriers)
        
        # Use scipy's interp1d with cubic spline
        interp_func = interp1d(XX, msignal, kind='cubic', fill_value='extrapolate')
        new_signal = interp_func(XV)
        
        # FFT and shift
        XF_int = np.fft.fft(new_signal, numb_carriers)
        XF_int = np.fft.fftshift(XF_int)
        
        # Update offset guess
        guess_offset = offset
        guess_offset = update_offset(guess_offset)
        
        if i > 0:
            if abs((abs(symbol_starts[i] - symbol_starts[i-1]) / 13440) - 1) < 0.2:
                offset_probability = 100
            else:
                offset_probability = 0
                print("this should not really happen")
        
        # Get carrier offset
        offset = get_carrier_offset(XF_int[0:60], guess_offset, offset_probability, 
                                    offset_order, offset_spacing)
        
        total_offsets[i] = offset
        
        # Get carrier locations
        carrier_locations = get_all_carrier_locations(offset, permanant_carriers, numb_carriers)
        
        # Extract pilot symbols
        X_pilots_int = XF_int[carrier_locations - 1]  # Adjust for 0-based indexing
        
        # Get channel state
        CS = get_channel_state(X_pilots_int, prbs_seq, carrier_locations)
        
        # Channel compensation
        XF_comp, Comps = channel_compensation(XF_int, CS, carrier_locations, numb_carriers)
        
        # Get actual symbols
        i_symbol, q_symbol = get_actual_symbols(XF_comp, numb_carriers, sypo, 
                                                carrier_locations, tps_carriers)
        
        # Map symbols to new locations
        XF_ideal = map_symbols_to_new_locations(i_symbol, q_symbol, numb_carriers, 
                                                prbs_seq, sypo, carrier_locations, tps_carriers)
        
        # Apply compensation
        XF_comp = XF_comp * Comps
        
        # Produce IQ data
        output_IQ = produce_IQ_data(XF_comp, numb_carriers, len(msignal), CS)
        
        # Apply phase correction
        corr = np.exp(1j * symbol_phases[i] / S / dt * t)
        output_IQ = output_IQ * corr
        
        # Scale output
        mmm = np.max(np.abs(check_me))
        m_small = np.max(np.abs(output_IQ))
        output_IQ = mmm / m_small * output_IQ
        
        # Assign to output signal (convert to 0-based indexing)
        output_signal[start_idx:end_idx] = output_IQ
        
        # Handle guard interval
        guard_start = symbol_starts[i] - 1
        guard_end = symbol_starts[i] + D
        output_signal[guard_start:guard_end] = output_IQ[S - D:]
    
    return output_signal, total_offsets


# Note: The following functions need to be defined separately:
# - generate_prbs(numb_carriers)
# - find_symbol_starts(input_signal, S, D, dt)
# - update_offset(guess_offset)
# - get_carrier_offset(XF_int, guess_offset, offset_probability, offset_order, offset_spacing)
# - get_all_carrier_locations(offset, permanant_carriers, numb_carriers)
# - get_channel_state(X_pilots_int, prbs_seq, carrier_locations)
# - channel_compensation(XF_int, CS, carrier_locations, numb_carriers)
# - get_actual_symbols(XF_comp, numb_carriers, sypo, carrier_locations, tps_carriers)
# - map_symbols_to_new_locations(i_symbol, q_symbol, numb_carriers, prbs_seq, sypo, carrier_locations, tps_carriers)
# - produce_IQ_data(XF_comp, numb_carriers, length, CS)





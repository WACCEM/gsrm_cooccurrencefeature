import time
import calendar
import numpy as np
import pandas as pd

#-------------------------------------------------------------------
def convert_cftime_to_standard_calendar(cftime_times):
    """
    Convert cftime datetime objects (non-standard calendars) to numpy datetime64 (standard calendar).
    
    This function handles conversion from calendars like DatetimeNoLeap, Datetime360Day, etc.
    to standard proleptic_gregorian calendar (numpy.datetime64).
    
    Args:
        cftime_times: array-like of cftime datetime objects
            Timestamps with non-standard calendar (e.g., DatetimeNoLeap, Datetime360Day)
            
    Returns:
        numpy.ndarray: Array of numpy.datetime64 objects with standard calendar
    """
    
    # Check if input is a single timestamp
    is_single_object = not hasattr(cftime_times, '__iter__')
    
    # Convert to list for uniform processing
    times_list = [cftime_times] if is_single_object else cftime_times
    
    # Convert each cftime timestamp to numpy datetime64
    converted_times = []
    for t in times_list:
        # Check if it's already a standard datetime type
        if isinstance(t, (np.datetime64, pd.Timestamp)):
            converted_times.append(np.datetime64(t, 'ns'))
        # If it's a cftime object, extract components and create datetime64
        elif hasattr(t, 'year'):
            # Create a pandas Timestamp from the cftime components
            # Note: This may shift dates for non-standard calendars that have different
            # day counts (e.g., Feb 30 in 360_day calendar doesn't exist in standard)
            try:
                pd_time = pd.Timestamp(
                    year=t.year, month=t.month, day=t.day,
                    hour=t.hour, minute=t.minute, second=t.second
                )
                converted_times.append(pd_time.to_datetime64())
            except ValueError as e:
                # Handle invalid dates (e.g., Feb 30)
                # For simplicity, we'll skip invalid dates or adjust them
                print(f"Warning: Could not convert {t} to standard calendar: {e}")
                # Try to adjust the day to the last valid day of the month
                last_day = calendar.monthrange(t.year, t.month)[1]
                adjusted_day = min(t.day, last_day)
                pd_time = pd.Timestamp(
                    year=t.year, month=t.month, day=adjusted_day,
                    hour=t.hour, minute=t.minute, second=t.second
                )
                converted_times.append(pd_time.to_datetime64())
        else:
            raise TypeError(f"Unsupported time type: {type(t)}")
    
    # Return a single object or array based on input type
    if is_single_object:
        return converted_times[0]
    else:
        return np.array(converted_times, dtype='datetime64[ns]')
    
#-------------------------------------------------------------------
def convert_to_matching_calendar(std_times, target_calendar):
    """
    Convert standard calendar (proleptic_gregorian) timestamps to match a target calendar.
    
    Args:
        std_times: array of numpy.datetime64, pandas.DatetimeIndex or pandas.Timestamp
            Timestamps with standard (proleptic_gregorian) calendar
        target_calendar: str
            Target calendar to convert to ('365_day', '360_day', 'noleap', etc.)
            
    Returns:
        cftime.datetime objects using the target calendar
    """

    
    # Initialize the appropriate cftime date type based on target calendar
    calendar_types = {
        '365_day': cftime.DatetimeNoLeap,
        'noleap': cftime.DatetimeNoLeap,
        '360_day': cftime.Datetime360Day,
        'all_leap': cftime.DatetimeAllLeap,
        'julian': cftime.DatetimeJulian,
        # Add other calendars as needed
    }
    
    if target_calendar in ['proleptic_gregorian', 'gregorian', 'standard']:
        # No conversion needed
        return std_times
    
    if target_calendar not in calendar_types:
        raise ValueError(f"Unsupported calendar: {target_calendar}")
        
    datetime_type = calendar_types[target_calendar]
    
    # Check if input is a single timestamp
    is_single_object = not hasattr(std_times, '__iter__') or isinstance(std_times, pd.Timestamp)
    
    # Convert to list for uniform processing
    times_list = [std_times] if is_single_object else std_times
    
    # Convert each timestamp to the target calendar
    converted_times = []
    for t in times_list:
        # Convert numpy.datetime64 to pandas.Timestamp which has the necessary attributes
        if isinstance(t, np.datetime64):
            ts = pd.Timestamp(t)
            converted_times.append(datetime_type(
                ts.year, ts.month, ts.day, 
                ts.hour, ts.minute, ts.second
            ))
        else:
            # For pandas.Timestamp or datetime objects that already have year, month attributes
            converted_times.append(datetime_type(
                t.year, t.month, t.day, 
                t.hour, t.minute, t.second
            ))
    
    # Return a single object or a list based on input type
    if is_single_object:
        return converted_times[0]
    else:
        return converted_times
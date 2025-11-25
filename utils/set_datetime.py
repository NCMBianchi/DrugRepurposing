"""
SET_DATETIME UTILITY MODULE: set the date
Created on February 21st 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]
"""

import datetime

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def set_date(date_tog, date_or = datetime.date(2024,7,31)):
    """
    Sets the date to today or a custom date depending on a toggle variable.

    :param date_toggle: 1 = override, 0 = today.
    :param date_override: which date to use as an override, default = July 31st 2024.
    :return: ___.
    """
    if date_tog == 0:
        today = datetime.date.today()
    elif date_tog == 1:
        today = date_or

    return today
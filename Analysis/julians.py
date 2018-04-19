import datetime
import jdcal
import math as m
""" 
Simulation start time: Oct 1st, 2015 
"""
yy = 0
mm = 0
dd = 0

date_factor = 1
if (100 * yy + mm - 190002.5) < 0:
    date_factor = -1

# simulation start time in JD (Julian Day)
jd_start = 367 * yy - m.trunc(7 * (yy + m.trunc((mm + 9) / 12)) / 4) + m.trunc(
    (275 * mm) / 9) + dd + 1721013.5 - 0.5 * date_factor

# dt = datetime.date(yy, mm, dd)
# jd_start2 = int(sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)))

# Test both method equivalent
# print(jd_start == jd_start2)


print(jd_start)
# dat = datetime.datetime.strptime('16234', '%y%j').date()



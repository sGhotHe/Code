import numpy as np
import sys

def month_days(year):
    """
    Returns a list of days of for each month in this year
    """
    result = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    if divmod(year, 4)[1] == 0:
        if divmod(year, 100)[1] == 0 and divmod(year, 400)[1] > 0:
            result[1] = 28
        else:
            result[1] = 29

    return result

def date2jul(year, month, day, **args):
    """
    Used to get the julian day 
    ---Keywords------------
       year   :      year
       month  :      month
       day    :      day
    ---Arguments-----------
       hour   :      hour
       minute :      minute
       second :      second
       microsecond : microsecond
    ---Return-----
       jul    :      julian day
    """
    if "hour" in args:
    	hour = args["hour"]
    else:
    	hour = 0
    if "minute" in args:
    	minute = args["minute"]
    else:
    	minute = 0
    if "second" in args:
    	second = args["second"]
    else:
    	second = 0
    if "microsecond" in args:
    	microsecond = args["microsecond"]
    else:
    	microsecond = 0
    mds = month_days(year)
    if month < 2:
        jul = day
    else:
        jul = mds[0 : month - 1].sum() + day
    jul = jul + hour / 24.0
    jul = jul + minute / 24.0 / 60.0
    jul = jul + second / 24.0 / 60.0 / 60.0
    jul = jul + microsecond / 1000000.0 / 24.0 / 60.0 / 60.0
    return jul
    
def jul2date(year, jul, **args):
    """ 
    Used to transform the julian day to date
    ---keywords-----------
        year   :  year
        jul    :  julian day
    ---Arguments----------------
        form   : which form to return
    ---Return-------------------
       form=1: return yyyy-mm-dd : HH:MM:SS
       form=2: return yyyy-mm-dd
       if keyword form is not defined then return [year,month,day,hour,minute,second] 
    """
    '''
    if jul > 365:
        jul -= 365
    '''
    # deleted by Hebs 21/3/7/15:43
    mds = month_days(year)
    start = 0
    nd = int(jul)
    hour = int((jul - nd) * 24.0)
    minute = int((jul - nd - hour / 24.0) * 24.0 * 60)
    second = int((jul - nd - hour / 24.0 - minute / (24.0 * 60.0)) * (24.0 * 60 * 60))
    for m in np.arange(12):
        start = start + mds[m]
        if start >= nd:
            month = m + 1
            day = nd - (start - mds[m])
            break
        else:
            if m == 11:
                print(("The jul:", jul, ",is not valid in year ", year))
                return -1
    if "form" in args:
        form = args["form"]

        if args["form"] == False:
            return [year, month, day]
        else:
            if month == -1:
                return month
            else:
                if form == 1:
                    return str(
                        "%04d-%02d-%02d : %02d:%02d:%02d"
                        % (year, month, day, hour, minute, second)
                    )
                else:
                    return str("%04d-%02d-%02d" % (year, month, day))
    else:
        return [year, month, day, hour, minute, second]
    
def jul2strDate(year, jul, **args):
    """
    used to transform the julian day to date of str type of date

    ------input------
    year  : year
    jul   : julian day of the year

    ------args-------
    form  : whichn form to return.   considering the date is 2016.9.11 12:20:50
            if form=='M.D'   : then return 9.11
            if form=='Y.M.D' : then return 2016.9.11
            if form=='h:m:s' : then return 12:20:50
            if form=='M/Y/D h.m.s' : then return 9/2016/11 12.20.50
    ------return------
    sdate : str of date.
    """
    if "form" in args:
        form = args["form"]
    else:
        form = "M.D"
    [year, month, day, hour, minute, second] = jul2date(year, jul)
    year = year % 100
    if "Y" in form:
        form = form.replace("Y", str(year))
    if "M" in form:
        form = form.replace("M", "%02d" % month)
    if "D" in form:
        form = form.replace("D", "%02d" % day)
    if "h" in form:
        form = form.replace("h", "%02d" % hour)
    if "m" in form:
        form = form.replace("m", "%02d" % minute)
    if "s" in form:
        form = form.replace("s", "%02d" % second)
    sdate = form
    return sdate

def strdate2jul(strdate, **args):
    """
    Used to Convert string of a date to julian day
    ---Keywords----------------------------
        strdate  :   string of the date, form: 'yyyy-mm-dd' or 'yyyymmdd'
    ---Arguments---------------------------
        form     :   The string form of date, form=1: 'yyyy-mm-dd', form=2: 'yyyymmdd'
                     default: form=1
    ---Return------------------------------
        year     :   The year
        jul      :   The julian day
    """
    if "form" in args:
        form = args["form"]
    else:
        form = 1
    if form == 1:
        try:
            date = [float(x) for x in strdate.split("-")]
            year = date[0]
            month = date[1]
            day = date[2]
        except:
            print("The string form is not right")
            return False
    else:
        try:
            year = int(strdate[0:4])
            month = int(strdate[4:6])
            day = int(strdate[6:])
        except:
            print("The string form is not right")
            return False, False
    jul = date2jul(year, month, day)
    return year, jul

def iday2month(iday):
	'''
	This function is to use the day in a standard year iday to calculate which month it in
	input:
		iday    : day in a standard year, int, 1 to 365
	output:
		month   : month, int
	'''
	if iday<1 or iday>365:
		print('Illegal iday. Please check.')
		sys.exit()
	month_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	for i in range(len(month_day)):
		if sum(month_day[:i+1])>=iday:
			return i+1

if __name__ == '__main__':
	'''
	print(jul2strDate(2020, 366.72916667, form="Y:M:D:h:m:s"))
	print(date2jul(2021, 1, 1, hour=18))
	'''
	print(iday2month(365))
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
